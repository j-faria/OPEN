program main

	use params
	use nestwrapper
	use like

	use strings
	use gputils

	implicit none

	integer i
	integer :: iou, ierr
	integer :: n
	character(len=40) line
	character(len=1) delim
	character(len=5), dimension(10) :: args
    character(len=100) :: fmt
	integer :: nargs, k, ind1, ind2
	integer, dimension(:), allocatable :: n_each_observ

	real(kind=8), parameter :: kmax = 2129d0 ! m/s
	real(kind=8), dimension(6) :: a

    namelist /NEST_parameters/ sdim, &
                               nest_IS, nest_updInt, nest_resume, nest_maxIter, nest_fb, &
                               nest_root, nest_context
    ! read configuration values from namelist
    iou = 8
    open(unit=iou, file="./OPEN/multinest/namelist1", status='old', action='read', delim='quote', iostat=ierr)
    if (ierr /= 0) then
      stop 'Failed to open namelist file'
    else
      read(iou, nml=NEST_parameters)
      close(iou)
    end if

    !! check compatibility between dimensionality and context
    if ((nest_context == 111 .and. sdim < 6) .or. &
        (nest_context == 112 .and. sdim < 7) .or. &
        (nest_context == 121 .and. sdim < 11) .or. &
        (nest_context == 122 .and. sdim < 12) ) then
		stop 'Conflict between "sdim" and "nest_context"'
    end if

    if (mod(nest_context, 10) == 2) then
    	if (nest_context / 100 == 2) then
#ifdef MPI
    		call MPI_INIT(ierr)
    		call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
    		if (i==0) print *, '==> GP and jitter are incompatible, nest_context cannot be 2.y.2'
    		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    		call MPI_ABORT(MPI_COMM_WORLD,1)
#else
			print *, '==> GP and jitter are incompatible, nest_context cannot be 2.y.2'
			STOP 1
#endif
    	end if
    	using_jitter = .true.
    	nextra = 1
    else
    	using_jitter = .false.
    	nextra = 0
    	if (nest_context / 100 == 2) then
    		using_gp = .true.
    		nextra = 2  ! hyperparameters of GP
    	end if
    end if

	!doing debug?
	doing_debug = .false.


	!! load data
	open(unit=15, file='input.rv', status="old")
	read(15, *) line ! header

	read(15, '(a)') line
	line = trim(line(index(line, "#")+1:))
	read(line, *) nobserv ! number of observatories
    if (nobserv < 1) stop 'Error parsing input file'

	allocate(n_each_observ(nobserv)) ! number of measurements in each observatory
	read(15, '(a)') line
	line = trim(line(index(line, "#")+1:)) ! remove the #
	delim = ' '
	call parse(line, delim, args, nargs) ! split whitespaces
	if (nargs /= nobserv) stop 'Error parsing input file'
	do k=1, nargs
		read( args(k), '(i5)' ) n_each_observ(k) ! store
	end do

	read(15, '(a)') line
	line = trim(line(index(line, "#")+1:))
	read(line, *) n  ! total number of measurements
	if (n /= sum(n_each_observ)) stop 'Error parsing input file'


    !! some parameters and allocations depend on the ones set on the namelist / input file
    nplanets = mod(nest_context,100) / 10 ! this is good for now
    sdim = 5*nplanets + nobserv + nextra
    allocate(spriorran(sdim, 2))
    allocate(nest_pWrap(sdim))
    nest_nPar=sdim


	!initialize data size dependant variables
	call likelihood_init(n)
	do i = 1, n
        read(15, *) times(i), rvs(i), errors(i)
    end do
    
    if (nobserv == 1) then ! if only one observatory, observ is always 1
    	observ = 1
    else ! else it takes a different value for each observatory
    	observ(1:n_each_observ(1)) = 1
    	do i = 1, nobserv-1
    		ind1 = sum(n_each_observ(:i))
    		ind2 = sum(n_each_observ(:i+1))
    		observ(ind1+1:ind2) = i+1
    	end do
    end if

    ! extra parameters are systematic velocities for each observatory 
    ! plus, if present, the jitter or the hyperparameters
    nextra = nextra + nobserv


    if (using_gp) then
        !write(*,*) 'Using Gaussian Process model'
    	!k3 = DiagonalKernel((/1.d0/))
    	!kernel_to_pass => k3
    	!k5 = ExpSquaredKernel((/ 1.d0, 1.d0 /))
    	!kernel_to_pass => k5
        k6 = ExpSineSquaredKernel((/ 1.d0, 25.d0 /))
        kernel_to_pass => k6

    	gp1 = GP(k6%evaluate_kernel(times, times), kernel_to_pass)
    	gp1%mean_fun => mean_fun_keplerian
        !print *, gp1%gp_kernel%pars
		gp_n_planets = 1
		gp_n_planet_pars = 6
    end if

	!no parameters to wrap around
	nest_pWrap = 0
	
	! here we set prior limits, 
	! the mathematical form is only used when rescaling
	do i = 1, sdim-nextra, 5
		!! Period, Jeffreys, 0.2d - 365000d
		spriorran(i,1)= 0.2d0 !0.2d0
		spriorran(i,2)= 365000d0 !365000.d0

		!! semi amplitude K, Mod. Jeffreys
		spriorran(i+1,1)=1d0
		! since the upper limit depends on e and P, it can only be set
		! when rescaling. We just initialize it here to a big number
		spriorran(i+1,2)=500d0

		!! eccentricity, Uniform, 0-1
		spriorran(i+2,1)=0d0
		spriorran(i+2,2)=1d0

		!! long. periastron, Uniform, 0-2pi rad
		spriorran(i+3,1)=0d0
		spriorran(i+3,2)=twopi		

		!! chi, Uniform, 
		spriorran(i+4,1)= minval(times)
		spriorran(i+4,2)= spriorran(5,1) + spriorran(1,2)
	end do

	if (using_jitter) then
    ! parameter array organization in this case:
    ! P1, K1, ecc1, omega1, t01, [P2, K2, ecc2, omega2, t02], jitter, vsys_obs1, [vsys_obs2]
		i = sdim-nextra+1 ! index of jitter parameter
		spriorran(i,1)= 1d0
		spriorran(i,2)= kmax

		!! systematic velocity(ies), Uniform, -kmax - kmax
		i = sdim-nextra+2 
		spriorran(i:,1)= -kmax
		spriorran(i:,2)= kmax

    else if (using_gp) then
    ! parameter array organization in this case:
    ! P1, K1, ecc1, omega1, t01, [P2, K2, ecc2, omega2, t02], vsys_obs1, [vsys_obs2], hyperpar1, hyperpar2
        !! systematic velocity(ies), Uniform, -kmax - kmax
        i = sdim-nextra+1
        spriorran(i:,1)= -kmax
        spriorran(i:,2)= kmax

        !! hyperparameters
        i = sdim-nextra+nobserv+1
        spriorran(i,1)= 0.9d0
        spriorran(i,2)= 1.1d0

        spriorran(i+1:,1)= 24.d0
        spriorran(i+1:,2)= 26.d0


	else
    ! parameter array organization in this case:
    ! P1, K1, ecc1, omega1, t01, [P2, K2, ecc2, omega2, t02], vsys_obs1, [vsys_obs2]
		!! systematic velocity(ies), Uniform, -kmax - kmax
		i = sdim-nextra+1
		spriorran(i:,1)= -kmax
		spriorran(i:,2)= kmax

	end if	

!    call MPI_INIT(ierr)
!    call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
!    if (i==0) then
!        write(fmt,'(a,i2,a)')  '(',sdim,'f13.4)'
!        write(*, fmt) spriorran(:,1)
!        write(*, fmt) spriorran(:,2)
!    end if
!    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!
! 	stop "we have to stop here"
   	call nest_Sample

   	! deallocate memory
   	!call nullify(kernel_to_pass)
   	deallocate(n_each_observ)
   	call likelihood_finish

contains

!**********************************************************************
!   Find number of uncommented lines of ASCII file
!**********************************************************************
  integer function n_uncomment_lines(filename, comment)

    character(len=*),intent(in) :: filename ! the ascii file 
    character(len=*),intent(in) :: comment  ! the comment character
    character(len=80) :: line  
    integer :: count_total    ! how many lines have been read so far
    integer :: count_uncomment
    integer :: un
    integer :: ioerr    ! used for I/O errors

    un = 18
    open(unit=un, file = filename, status='old')
    do
        read(un, '(a)',iostat=ioerr) line
        if(ioerr.ne.0) exit
        if(line(1:1) == comment) then
            count_total = count_total + 1
        else
            count_uncomment = count_uncomment + 1
        end if      
    end do
    close(un)
    n_uncomment_lines = count_uncomment

  end function n_uncomment_lines
end program


