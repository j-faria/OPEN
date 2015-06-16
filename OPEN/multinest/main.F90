program main

    use params
    use nestwrapper
    use like

    use strings
    use gputils
    use array_utils, only: linspace, variance

    implicit none

    integer :: i, iou, ierr
    integer :: PGOPEN
    integer :: n

    character(len=100) cmdline
    character(len=40) line
    character(len=1) delim
    character(len=5), dimension(10) :: args
    character(len=100) :: fmt
    integer :: nargs, k, ind1, ind2
    integer, dimension(:), allocatable :: n_each_observ

    real(kind=8), parameter :: kmax = 2129d0 ! m/s

    namelist /NEST_parameters/ sdim, nest_nlive, &
                               nest_IS, nest_updInt, nest_resume, nest_maxIter, nest_fb, nest_MAPfb, nest_liveplot, &
                               nest_root, nest_context, &
                               training, trained_parameters, trained_std, &
                               trend, trend_degree, &
                               lin_dep, n_lin_dep

    !MPI initializations
#ifdef MPI
    call MPI_INIT(ierr)
    if (ierr/=MPI_SUCCESS) then
            write(*,*)'Error starting MPI. Terminating.'
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
    end if
#endif

    !! read configuration values from namelist
    ! namelist file from first cmd line argument
    call get_command_argument(1, cmdline, status=ierr)
    if (ierr /= 0) error stop 'Failed to retrieve namelist filename from cmd line'

    iou = 8
    open(unit=iou, file=trim(cmdline), status='old', action='read', delim='quote', iostat=ierr)
    if (ierr /= 0) then
      error stop 'Failed to open namelist file'
    else
      read(iou, nml=NEST_parameters)
      close(iou)
    end if


#ifdef PLOT 
    ! initialize graphics device
    if (nest_liveplot) then
#ifdef MPI
        call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
        if (i==0) then
            IF (PGOPEN('/XSERVE') .LE. 0) STOP 'Error creating plot window'
            CALL PGPAGE
            CALL PGBBUF
        end if
#else
    IF (PGOPEN('/XSERVE') .LE. 0) STOP 'Error creating plot window'
    CALL PGPAGE
    CALL PGBBUF   
#endif
    end if
#endif

!    !! check compatibility between dimensionality and context
!    if ((nest_context == 111 .and. sdim < 6) .or. &
!        (nest_context == 112 .and. sdim < 7) .or. &
!        (nest_context == 121 .and. sdim < 11) .or. &
!        (nest_context == 122 .and. sdim < 12) ) then
!        stop 'Conflict between "sdim" and "nest_context"'
!    end if

    trend = .false.
    if (nest_context/1000 > 0) then
        trend = .true.
        trend_degree = nest_context/1000
        nest_context = nest_context - (nest_context / 1000)*1000
    end if

    if (mod(nest_context, 10) == 2) then
        if (nest_context / 100 == 2) then
#ifdef MPI
            call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
            if (i==0) then
                print *, '==> GP and jitter are incompatible, nest_context cannot be 2.y.2'
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_ABORT(MPI_COMM_WORLD,1)
            end if
#else
            stop '==> GP and jitter are incompatible, nest_context cannot be 2.y.2'
#endif
        end if
        using_jitter = .true.
        nextra = 1
    else
        using_jitter = .false.
        nextra = 0
        if (nest_context / 100 == 2 .or. nest_context / 100 == 3) then
            ! if 2, using GP plus planets; if 3, using GP only, without planets
            using_gp = .true.
            nextra = 5  ! hyperparameters of GP
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
    if (nobserv < 1) stop 'Error parsing input file (1)'

    allocate(n_each_observ(nobserv)) ! number of measurements in each observatory
    read(15, '(a)') line
    line = trim(line(index(line, "#")+1:)) ! remove the #
    delim = ' '
    call parse(line, delim, args, nargs) ! split whitespaces
    if (nargs /= nobserv) stop 'Error parsing input file (2)'
    do k=1, nargs
        read( args(k), '(i5)' ) n_each_observ(k) ! store
    end do

    read(15, '(a)') line
    line = trim(line(index(line, "#")+1:))
    read(line, *) n  ! total number of measurements
    if (n /= sum(n_each_observ)) stop 'Error parsing input file (3)'


    !! some parameters and allocations depend on the ones set on the namelist / input file
    nextra = nextra + trend_degree 
    nplanets = mod(nest_context,100) / 10 ! this is good for now
    sdim = 5*nplanets + nobserv + nextra
    allocate(spriorran(sdim, 2))
    allocate(nest_pWrap(sdim))
    nest_nPar=sdim

    !! the number of parameters to cluster is usually set to 3
    !! but it cannot be greater than sdim
    if (sdim < 4) then
        nest_nClsPar=1
    else
        nest_nClsPar=5
    end if


    !! initialize data size dependant variables
    call likelihood_init(n)


    !include linear dependences in the model
    if (lin_dep) then
        !nextra = nextra + n_lin_dep
        if (n_lin_dep == 1) then
            allocate(linvar1(n))
        else if (n_lin_dep == 2) then
            allocate(linvar1(n))
            allocate(linvar2(n))
        else if (n_lin_dep == 3) then
            allocate(linvar1(n))
            allocate(linvar2(n))
            allocate(linvar3(n))
        end if
    end if
    
    do i = 1, n
        if (lin_dep) then
            if (n_lin_dep == 1) read(15, *) times(i), rvs(i), errors(i), linvar1(i)
            if (n_lin_dep == 2) read(15, *) times(i), rvs(i), errors(i), linvar1(i), linvar2(i)
            if (n_lin_dep == 3) read(15, *) times(i), rvs(i), errors(i), linvar1(i), linvar2(i), linvar3(i)
        else
            read(15, *) times(i), rvs(i), errors(i)
        end if
    end do

    if (size(times) <= sdim) then
#ifdef MPI
        call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
        if (i==0) print *, 'Number of observations is less than dimensionality. Aborting!'
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (i==0) call MPI_ABORT(MPI_COMM_WORLD,1)
#else
        stop 'Number of observations is less than dimensionality. Aborting!'
#endif
    end if

    ! calculate the mean time to be used in the trend
    mean_time = sum(times) / size(times)


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

#ifdef PLOT 
    !! build oversampled time array
    !print *, minval(times), maxval(times), 5*n
    times_oversampled = linspace(minval(times), maxval(times), 5*n)
#endif

    !! initialize the GP "object"
    if (using_gp) then
        
        !call init_HODLR(N)

        ! this is the jitter, a white noise homoskedastic kernel
        k3 = DiagonalKernel((/1.d0/))

        ! the quasi periodic terms
        k5 = ExpSquaredKernel((/ 1.d0 /)) ! pars(1) --> lengthscale
        k6 = ExpSineSquaredKernel((/ 25.d0, 1.d0 /)) ! pars(1) --> period; pars(2) --> correlation scale

        ! the product of the two kernels
        k8 = ProductKernels((/1.d0, 1.d0 /)) ! multiplying constants for kernel1 and kernel2
        k8%kernel1 => k5
        k8%kernel2 => k6

        ! the sum of the jitter and quasi-periodic kernels
        k7 = SumKernels((/1d0, 1d0/)) ! multiplying constants for kernel1 and kernel2
        k7%kernel1 => k8 ! quasi-periodic
        k7%kernel2 => k3 ! jitter

        ! final kernel
        kernel_to_pass => k7
        ! we also want the GP object to have information about its subkernels
        subk1 => k8
        subk2 => k3

        gp1 = GP(k7%evaluate_kernel(times, times), kernel_to_pass, subk1, subk2)
!         if (nest_context / 100 == 2) then
!             gp1%mean_fun => mean_fun_keplerian
!         else
            gp1%mean_fun => mean_fun_constant
!         endif
        gp_n_planets = nplanets
        gp_n_planet_pars = 5*nplanets + nobserv
    end if

    !! extra parameters are systematic velocities for each observatory 
    !! plus, if present, the jitter, trends or the hyperparameters
    nextra = nextra + nobserv


    !no parameters to wrap around
    nest_pWrap = 0
    
    ! here we set prior limits, 
    ! the mathematical form is only used when rescaling
    do i = 1, sdim-nextra, 5
        !! Period, Jeffreys, 0.2d - 365000d
        spriorran(i,1)= 1d0 !0.2d0
!         spriorran(i,2)= 365000d0
       spriorran(i,2)= 1.d0*(maxval(times) - minval(times)) ! don't look for periods longer than timespan

        !! semi amplitude K, Mod. Jeffreys (or Uniform)
!        spriorran(i+1,1)=1.d0
        spriorran(i+1,1)= sum(errors) / n
!         spriorran(i+1,2)= maxval(rvs) - minval(rvs)
        ! the true upper limit depends on e and P, and it will only be set when rescaling.
!        spriorran(i+1,2)=kmax 
        spriorran(i+1,2)=10.d0

        !! eccentricity, Uniform, 0-1
!         spriorran(i+2,1)=0d0
!         spriorran(i+2,2)=1d0
        !! eccentricity, Beta(0.867, 3.03), based on Kipping (2013)
        spriorran(i+2,1)=0.867d0
        spriorran(i+2,2)=3.03d0

        !! long. periastron, Uniform, 0-2pi rad
        spriorran(i+3,1)=0d0
        spriorran(i+3,2)=twopi
        nest_pWrap(i+3) = 1  ! wrap around this parameter

        !! chi, Uniform, 
        spriorran(i+4,1)= minval(times)
        spriorran(i+4,2)= spriorran(5,1) + spriorran(1,2)
!         spriorran(i+4,1)=2453988.24856
!         spriorran(i+4,2)=2453988.24896
    end do

    if (using_jitter) then
    ! parameter array organization in this case:
    ! P1, K1, ecc1, omega1, t01, [P2, K2, ecc2, omega2, t02], jitter, vsys_obs1, [vsys_obs2]
        i = sdim-nextra+1 ! index of jitter parameter
        !spriorran(i,1)= 1d0
        spriorran(i,1)= sum(errors) / n
        spriorran(i,2)= 10d0 !kmax
!         spriorran(i,2)= maxval(rvs) - minval(rvs)
        ! when data is in km/s
!         spriorran(i,1)= 0.001d0
!         spriorran(i,2)= kmax * 1d-3
        

        !! systematic velocity(ies), Uniform, -kmax - kmax
        i = sdim-nextra+2 
        spriorran(i,1)= minval(rvs) ! -kmax
        spriorran(i,2)= maxval(rvs) ! kmax

        !! trend 
        if (trend_degree > 0) then
            i = sdim-nextra+3
            spriorran(i,1)= 0d0
            spriorran(i,2)= 0.1d0
        end if

    else if (using_gp) then
    ! parameter array organization in this case:
    ! P1, K1, ecc1, omega1, t01, [P2, K2, ecc2, omega2, t02], vsys_obs1, [vsys_obs2], hyperpar1, hyperpar2
        !! systematic velocity(ies), Uniform, -kmax - kmax
        i = sdim-nextra+1
        spriorran(i,1)= -kmax
        spriorran(i,2)= kmax

        !! hyperparameters
        if (training) then 
            ! the prior limits in this case are just the values of the 
            ! trained parameters, only the GP amplitude is left to vary
            
            ! amplitude of GP
            i = sdim-nextra+nobserv+1
            spriorran(i,1)= 0.d0
            spriorran(i,2)= maxval(rvs) - minval(rvs)
            ! jitter 
            spriorran(i+1,1)= 0.d0
            spriorran(i+1,2)= 0.d0 !maxval(rvs) - minval(rvs) !maxval(errors)
            ! timescale for growth / decay of active regions (d)
            spriorran(i+2,1)= trained_parameters(3)
            spriorran(i+2,2)= trained_std(3)
            ! periodicity (recurrence) timescale -> rotation period of the star
            spriorran(i+3,1)= trained_parameters(4)
            spriorran(i+3,2)= trained_std(4)
            ! smoothing parameter (?)
            spriorran(i+4,1)= trained_parameters(5)
            spriorran(i+4,2)= trained_std(5)
        else
            ! amplitude of GP
            i = sdim-nextra+nobserv+1
            spriorran(i,1)= 0.d0
            spriorran(i,2)= maxval(rvs) - minval(rvs)
            ! jitter 
            spriorran(i+1,1)= 0.d0
            spriorran(i+1,2)= 0.d0
            ! timescale for growth / decay of active regions (d)
            spriorran(i+2,1)= 1d0
            spriorran(i+2,2)= 1000.d0
            ! periodicity (recurrence) timescale -> rotation period of the star
            spriorran(i+3,1)= 10.d0
            spriorran(i+3,2)= 20.d0
            ! smoothing parameter (?)
            spriorran(i+4,1)= 0.1d0
            spriorran(i+4,2)= 10d0
        end if
        
    else
    ! parameter array organization in this case:
    ! P1, K1, ecc1, omega1, t01, [P2, K2, ecc2, omega2, t02], vsys_obs1, [vsys_obs2]
        !! systematic velocity(ies), Uniform, -kmax - kmax
        i = sdim-nextra+1
        spriorran(i,1)= minval(rvs) ! -kmax
        spriorran(i,2)= maxval(rvs) ! kmax

        !! trend 
        if (trend_degree > 0) then
            i = sdim-nextra+2
            spriorran(i,1)= 0d0
            spriorran(i,2)= 0.1d0
        end if

    end if    

!     call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
!     if (i==0) then
!         print *, 'Prior limits:'
!         write(fmt,'(a,i2,a)')  '(',sdim,'f13.4)'
!         write(*, fmt) spriorran(:,1)
!         write(*, fmt) spriorran(:,2)
!        !write(*,*) observ
!     end if
!     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    print *, training, trained_parameters
!     stop "we have to stop here"

    call nest_Sample

    ! deallocate memory
    !call nullify(kernel_to_pass)
    deallocate(n_each_observ)
    call likelihood_finish

#ifdef PLOT 
    if (nest_liveplot) then
    ! finish graphics devide
#ifdef MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
    if (i==0) then
        CALL PGEBUF
        CALL PGCLOS
    end if
#else
    CALL PGEBUF
    CALL PGCLOS
#endif
    end if
#endif

#ifdef MPI
    call MPI_FINALIZE(ierr) 
#endif

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


