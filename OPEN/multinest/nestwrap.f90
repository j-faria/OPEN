module nestwrapper

! Nested sampling includes

use Nested
use params
use like
use priors
use array_utils, only: linspace, get_diagonal, sort
use gputils
   
implicit none

contains


	subroutine nest_Sample
	! 'main' function that actually calls MultiNest
		implicit none
		
	   	integer nclusters ! total number of clusters found
		integer context
		integer maxNode ! variables used by the posterior routine
	   
	    context = nest_context  ! why do I have to do this here?

	   	! calling MultiNest
	   	call nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol, &
	   		         nest_efr,sdim,nest_nPar, nest_nClsPar,nest_maxModes, &
	   		         nest_updInt,nest_Ztol,nest_root,nest_rseed,nest_pWrap, &
	   		         nest_fb,nest_resume,nest_outfile,nest_initMPI, &
	   		         nest_logZero,nest_maxIter,getLogLike,dumper,context)

	end subroutine nest_Sample


	subroutine getLogLike(Cube,n_dim,nPar,lnew,context)
	! Wrapper around Likelihood Function which rescales parameters
	! and then calls the actual likelihood calculations
		implicit none
		
		! Input arguments
		integer n_dim ! total number of free parameters
		integer nPar ! total number of free plus derived parameters

		!Input/Output arguments
		! on entry has the ndim parameters in unit-hypercube
		! on exit has the physical parameters plus a copy of any
		! derived parameters you want to store
		double precision Cube(nPar)
		real(kind=8) :: P, K, ecc, omega, t0
		real(kind=8) :: P1, K1, ecc1, omega1, t01
		real(kind=8) :: P2, K2, ecc2, omega2, t02
		real(kind=8) :: Vsys
		 						
		! Output arguments
		double precision lnew ! loglikelihood
		integer context ! additional information user wants to pass

		! Local variables
		double precision kmax
		integer i, j
		

		! Transform parameters to physical space using assigned priors
		! Cube(1:nPar) = P, K, ecc, omega, t0  for each planet     |
		!				 jitter at the next to last position	   | - if jitter model				
		!                vsys at the last position                 |
		!
		! Cube(1:nPar) = P, K, ecc, omega, t0  for each planet						|
		!				 vsys for each observatory at the next to last positions	| - if GP model				
		!                hyperparameters at the last positions						|
		
		if (using_jitter) then
			! prior for jitter
			i = nPar-nextra+1 ! index of jitter parameter
			Cube(i) = ModJeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))

			! prior for systematic velocity
			i = npar-nextra+2 ! index of systematic velocity
			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2))			

		else if (using_gp) then
			! prior for systematic velocity
			i = nPar-nextra+1 
			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2))

			! prior for hyperparameters
			i = nPar-nextra+nobserv+1
			do j = i, i+3
				Cube(j) = UniformPrior(Cube(j), spriorran(j,1), spriorran(j,2))
			end do

		else
			! prior for systematic velocity
			i = npar-nextra+1 ! index of systematic velocity
			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2))
			!Cube(i+1) = UniformPrior(Cube(i+1), spriorran(i+1,1), spriorran(i+1,2))
		end if


		! prior for period(s)
		do i = 1, nPar-nextra, 5
			Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))
		end do

		! priors for ecc, omega, t0
		do i = 3, nPar-nextra, 5
			!Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2)) ! ecc
			Cube(i) = BetaPrior(Cube(i), spriorran(i,1), spriorran(i,2)) ! ecc
			Cube(i+1) = UniformPrior(Cube(i+1), spriorran(i+1,1), spriorran(i+1,2)) ! omega
			Cube(i+2) = UniformPrior(Cube(i+2), spriorran(i+2,1), spriorran(i+2,2)) ! t0
		end do

		! prior for semi-amplitude(s)
		! this has to be set after eccentricity has been rescaled. If the ecc prior is uniform, it 
		! does not matter but if not, we need here the rescaled value.
		! The upper limit is set to Kmax (Pmin / P)^1/3  * (1/ sqrt(1 − e^2))
		! see Eq. (14) of Gregory 2007 [2007MNRAS.374.1321G]
		do i = 2, nPar-nextra, 5
			!Cube(i+1) is this planet's eccentricity
			!Cube(i-1) is this planet's period
			!kmax = spriorran(i,2)
			kmax = spriorran(i,2) * (spriorran(i-1,1) / Cube(i-1))**(1/3.d0) * 1.d0/(sqrt(1-Cube(i+1)**2))
			!print *, kmax, Cube(i+1)
			Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), kmax)
			!print *, Cube(i)
		end do


		!call loglike function here 
		call slikelihood(Cube,lnew,context)

	end subroutine getLogLike


	subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, INSlogZ, logZerr, context, ending)
	! dumper routine, called after every updInt*10 iterations
	! and at the end of the sampling.
		implicit none
		integer nSamples ! number of samples in posterior array
		integer nlive ! number of live points
		integer nPar ! number of parameters saved (physical plus derived)
		! array containing the last set of live points:
		double precision, pointer :: physLive(:,:)	
		! array with the posterior distribution
		double precision, pointer :: posterior(:,:)	
		! array with mean, sigmas, maxlike & MAP parameters:
		double precision, pointer :: paramConstr(:)	
		double precision maxLogLike ! max loglikelihood value
		double precision logZ, INSlogZ ! log evidence
		double precision logZerr ! error on log evidence
		integer context ! any additional information user wants to pass
		logical ending ! is this the final call?
		! local variables
		double precision, dimension(size(times) + 1000) :: t, mu, std
		double precision, dimension(size(times) + 1000, size(times) + 1000) :: cov
		character(len=100) gppredictfile
		integer i, map_index
		character(len=100) :: fmt

		!write(fmt,'(a,i2,a)')  '(',nPar,'f13.4)'
		!write(*,fmt) paramConstr(nPar*3+1:nPar*4)

		write(*,*) ' '
		if (nplanets == 1 .and. using_gp) then  ! 1 planet + 4 hyper
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',nPar - 4,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*4 - 4)
			
			write(*,'(4a13)') (/"t1", "t2", "t3", "t4" /)
			write(fmt,'(a,i2,a)')  '(', 4, 'f13.4)'
			write(*,fmt) paramConstr(nPar*4 - 3:)

		else if (nplanets == 1) then  ! 1 planet
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',nPar,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*4)
		
		else if (nplanets == 2 .and. using_gp) then  ! 2 planets + 4 hyper
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*3+5)
			write(fmt,'(a,i2,a)')  '(',6,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+6:nPar*4 - 4)

			write(*,'(4a13)') (/"t1", "t2", "t3", "t4" /)
			write(fmt,'(a,i2,a)')  '(', 4, 'f13.4)'
			write(*,fmt) paramConstr(nPar*4 - 3:)

		else if (nplanets == 2) then  ! 2 planets
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*3+5)
			write(fmt,'(a,i2,a)')  '(',6,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+6:nPar*4)

		else if (nplanets == 3 .and. using_gp) then  ! 3 planets + 4 hyper
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*3+10)
			write(fmt,'(a,i2,a)')  '(',6,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+11:nPar*4 - 4)

			write(*,'(4a13)') (/"t1", "t2", "t3", "t4" /)
			write(fmt,'(a,i2,a)')  '(', 4, 'f13.4)'
			write(*,fmt) paramConstr(nPar*4 - 3:)

		else if (nplanets == 3) then  ! 3 planets
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*3+10)
			write(fmt,'(a,i2,a)')  '(',6,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+11:nPar*4)

		else if (nplanets == 0 .and. using_gp) then  ! 4 hyper (plus one systematic velocity)
			write(*,'(a13)') (/" vsys" /)
			write(fmt,'(a,i2,a)')  '(', 1, 'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1)
			write(*,'(4a13)') (/"t1", "t2", "t3", "t4" /)
			write(fmt,'(a,i2,a)')  '(', 4, 'f13.4)'
			write(*,fmt) paramConstr(nPar*3+2:)

		else if (nplanets == 0) then  ! one systematic velocity
			write(*,'(a13)') (/" vsys" /)
			write(fmt,'(a,i2,a)')  '(', 1, 'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1)

		end if
		write(*,*) ' '

		if (ending .and. using_gp) then
			gppredictfile = TRIM(nest_root)//'gp.dat'

			! prediction times are lineary spaced
			t(1:1000) = linspace(minval(times), maxval(times), 1000)
			! and put together with observed times
			t(1001:) = times
			call sort(t, size(t))  ! sort it

			!! set hyperparameters to their MAP values
			map_index = nPar*3 ! this is where the MAP parameters start in the array
			gp1%gp_kernel%pars(1) = paramConstr(map_index+gp_n_planet_pars+1)
			call gp1%gp_kernel%set_kernel_pars(1, (/paramConstr(map_index+gp_n_planet_pars+2)/) )
			call gp1%gp_kernel%set_kernel_pars(2, paramConstr(map_index+gp_n_planet_pars+3:))

			!print *, paramConstr(map_index+gp_n_planet_pars+1)
			!print *, paramConstr(map_index+gp_n_planet_pars+2)
			!print *, paramConstr(map_index+gp_n_planet_pars+3:)
			!print *, paramConstr(nPar*3+1:nPar*4-4)
			!print *, ''
			!write(*, '(10f13.6)') gp1%gp_kernel%evaluate_kernel(times, times)
			!print *, ''

			write(*,*) 'Calculating predictions...'
			call gp1%predict(times, rvs, paramConstr(nPar*3+1:nPar*4-4), t, mu, cov, yerr=errors)
			std = sqrt(get_diagonal(cov))

			write(*,*) 'Writing file ', gppredictfile
			open(unit=59, file=gppredictfile, form='formatted', status='replace')
			write(59, '(3f16.4)') (t(i), mu(i), std(i), i=1, size(times) + 1000)
			close(59)

		end if

	end subroutine dumper


end module nestwrapper
