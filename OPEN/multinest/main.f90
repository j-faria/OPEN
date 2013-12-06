program main

	use params
	use nestwrapper
	use like
      
	implicit none
	
	integer i
	integer, parameter :: n=119
	integer, parameter :: nobserv=1

	real(kind=8), parameter :: kmax = 2129d0 ! m/s
	real(kind=8) :: slhood

	!doing debug?
	doing_debug = .false.

	!initialize data size dependant variables
	call likelihood_init(n)

	!load data
	open(unit=15, file='input.rv', status="old")
	do i = 1, n
        read(15, *) times(i), rvs(i), errors(i)
    end do


	!no parameters to wrap around
	nest_pWrap = 0
	
	! here we set prior limits, 
	! the mathematical form is only used when rescaling

	!! Period, Jeffreys, 0.2d - 365000d
	spriorran(1,1)= 1000.d0 !0.2d0
	spriorran(1,2)= 2000.d0 !365000.d0
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
	spriorran(5,1)=2.45128604d6
	spriorran(5,2)=2.45148604d6
	!! systematic velocity, Uniform
	!! Vmin = -Kmax, Vmax = Kmax
	spriorran(6,1)= -kmax
	spriorran(6,2)= kmax		


! 	stop "we have to stop here"
   	call nest_Sample

end program
