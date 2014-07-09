program main

	use params
	use nestwrapper
	use like
      
	implicit none
	
	integer i
	integer :: iou, ierr
	integer :: n
	character(len=10) line
	integer, parameter :: nobserv=1

	real(kind=8), parameter :: kmax = 2129d0 ! m/s

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
    if ((nest_context == 11 .and. sdim /= 6) .or. &
        (nest_context == 12 .and. sdim /= 7) .or. &
        (nest_context == 21 .and. sdim /= 11) ) then
		stop 'Conflict between "sdim" and "nest_context"'
    end if


    if (mod(nest_context, 10) == 2) then
    	using_jitter = .true.
    else
    	using_jitter = .false.
    end if

    !! some parameters depend on the ones set on the namelist
    nplanets = sdim / 5
    allocate(spriorran(sdim, 2))
    allocate(nest_pWrap(sdim))
    nest_nPar=sdim

	!doing debug?
	doing_debug = .false.


	!load data
	open(unit=15, file='input.rv', status="old")
	read(15, *) line ! header
	read(15, '(a)') line ! number of values
	line = trim(line(index(line, "#")+1:))
	read(line, *) n 
	
	!initialize data size dependant variables
	call likelihood_init(n)
	do i = 1, n
        read(15, *) times(i), rvs(i), errors(i)
    end do


	!no parameters to wrap around
	nest_pWrap = 0
	
	! here we set prior limits, 
	! the mathematical form is only used when rescaling
	if (nest_context == 11) then
		!! Period, Jeffreys, 0.2d - 365000d
		spriorran(1,1)= 500d0 !0.2d0
		spriorran(1,2)= 1000d0 !365000.d0
		!! semi amplitude K, Mod. Jeffreys
		spriorran(2,1)=0d0
		! since the upper limit depends on e and P, it can only be set
		! when rescaling. We just initialize it here to a big number
		spriorran(2,2)=500d0
		!! eccentricity, Uniform, 0-1
		spriorran(3,1)=0d0
		spriorran(3,2)=1d0		
		!! long. periastron, Uniform, 0-2pi rad
		spriorran(4,1)=0d0
		spriorran(4,2)=twopi		
		!! chi, Uniform, 0-1
		spriorran(5,1)= minval(times)
		spriorran(5,2)= spriorran(5,1) + spriorran(1,2)
		!! systematic velocity, Uniform
		!! Vmin = -Kmax, Vmax = Kmax
		spriorran(6,1)= -kmax
		spriorran(6,2)= kmax

	else if (nest_context == 21) then
		!!-> planet 1
		!! Period, Jeffreys, 0.2d - 365000d
		spriorran(1,1)= 500d0 !0.2d0
		spriorran(1,2)= 1000d0 !365000.d0
		!! semi amplitude K, Mod. Jeffreys
		spriorran(2,1)=0d0
		! since the upper limit depends on e and P, it can only be set
		! when rescaling. We just initialize it here to a big number
		spriorran(2,2)=500d0
		!! eccentricity, Uniform, 0-1
		spriorran(3,1)=0d0
		spriorran(3,2)=1d0		
		!! long. periastron, Uniform, 0-2pi rad
		spriorran(4,1)=0d0
		spriorran(4,2)=twopi		
		!! chi, Uniform, 0-1
		spriorran(5,1)= minval(times)
		spriorran(5,2)= spriorran(5,1) + spriorran(1,2)

		!!-> planet 2
		!! Period, Jeffreys, 0.2d - 365000d
		spriorran(6,1)= 50d0 !0.2d0
		spriorran(6,2)= 200d0 !365000.d0
		!! semi amplitude K, Mod. Jeffreys
		spriorran(7,1)=0d0
		! since the upper limit depends on e and P, it can only be set
		! when rescaling. We just initialize it here to a big number
		spriorran(7,2)=50d0
		!! eccentricity, Uniform, 0-1
		spriorran(8,1)=0d0
		spriorran(8,2)=1d0		
		!! long. periastron, Uniform, 0-2pi rad
		spriorran(9,1)=0d0
		spriorran(9,2)=twopi		
		!! chi, Uniform, 0-1
		spriorran(10,1)= minval(times)
		spriorran(10,2)= spriorran(5,1) + spriorran(1,2)

		!! systematic velocity, Uniform
		!! Vmin = -Kmax, Vmax = Kmax
		spriorran(11,1)= -kmax
		spriorran(11,2)= kmax

	end if	


! 	stop "we have to stop here"
   	call nest_Sample

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


