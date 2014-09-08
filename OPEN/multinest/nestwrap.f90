module nestwrapper

! Nested sampling includes

use Nested
use params
use like
use priors
use array_utils, only: linspace, get_diagonal
   
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
		end if

		! prior for period(s)
		do i = 1, nPar-nextra, 5
			Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))
		end do

		! prior for semi-amplitude(s)
		do i = 2, nPar-nextra, 5
			Cube(i) = ModJeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))
		end do

		! priors for ecc, omega, t0
		do i = 3, nPar-nextra, 5
			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2)) ! ecc
			Cube(i+1) = UniformPrior(Cube(i+1), spriorran(i+1,1), spriorran(i+1,2)) ! omega
			Cube(i+2) = UniformPrior(Cube(i+2), spriorran(i+2,1), spriorran(i+2,2)) ! t0
		end do


!		if (nest_context == 11) then  ! 1 planet
!			P = UniformPrior(Cube(1), spriorran(1,1), spriorran(1,2))
!			ecc = UniformPrior(Cube(3), spriorran(3,1), spriorran(3,2))
!			omega = UniformPrior(Cube(4), spriorran(4,1), spriorran(4,2))
!			t0 = UniformPrior(Cube(5), spriorran(5,1), spriorran(5,2))
!			K = UniformPrior(Cube(2), spriorran(2,1), spriorran(2,2)) ! for now
!
!			Vsys = UniformPrior(Cube(6), spriorran(6,1), spriorran(6,2))
!
!			Cube(1:nPar) = (/P, K, ecc, omega, t0, Vsys/)
!
!		else if (nest_context == 21) then
!			P1 = UniformPrior(Cube(1), spriorran(1,1), spriorran(1,2))
!			ecc1 = UniformPrior(Cube(3), spriorran(3,1), spriorran(3,2))
!			omega1 = UniformPrior(Cube(4), spriorran(4,1), spriorran(4,2))
!			t01 = UniformPrior(Cube(5), spriorran(5,1), spriorran(5,2))
!			K1 = UniformPrior(Cube(2), spriorran(2,1), spriorran(2,2)) ! for now !!!!!
!
!			P2 = UniformPrior(Cube(6), spriorran(6,1), spriorran(6,2))
!			ecc2 = UniformPrior(Cube(8), spriorran(8,1), spriorran(8,2))
!			omega2 = UniformPrior(Cube(9), spriorran(9,1), spriorran(9,2))
!			t02 = UniformPrior(Cube(10), spriorran(10,1), spriorran(10,2))
!			K2 = UniformPrior(Cube(7), spriorran(7,1), spriorran(7,2)) ! for now !!!!!		
!
!			Vsys = UniformPrior(Cube(11), spriorran(11,1), spriorran(11,2))
!
!			Cube(1:nPar) = (/P1, K1, ecc1, omega1, t01, &
!			                 P2, K2, ecc2, omega2, t02, &
!			                 Vsys/)
!
!		end if
		!write(unit=*, fmt=*) '+++', Cube(1:nPar)

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
		double precision, dimension(1000) :: t, mu, std
		double precision, dimension(1000, 1000) :: cov
		character(len=100) gppredictfile
		integer i, map_index
		character(len=100) :: fmt

		! now do something
		!if (doing_debug) write(*,*) paramConstr(nPar*3+1:nPar*4-2)
		write(fmt,'(a,i2,a)')  '(',nPar,'f13.4)'
		write(*,fmt) paramConstr(nPar*3+1:nPar*4)
		!print *, ending, using_gp

		if (ending .and. using_gp) then
			gppredictfile = TRIM(nest_root)//'gp.dat'
			write(*,*) 'Writing file ', gppredictfile

			t = linspace(minval(times), maxval(times), 1000)

			!! set hyperparameters to their MAP values
			map_index = nPar*3 ! this is where the MAP parameters start in the array
			gp1%gp_kernel%pars(1) = paramConstr(map_index+gp_n_planet_pars+1)
			call gp1%gp_kernel%set_kernel_pars(1, (/paramConstr(map_index+gp_n_planet_pars+2)/) )
			call gp1%gp_kernel%set_kernel_pars(2, paramConstr(map_index+gp_n_planet_pars+3:))

			!print *, paramConstr(nPar*3+1:nPar*4-2)
			call gp1%predict(times, rvs, paramConstr(nPar*3+1:nPar*4-4), t, mu, cov, yerr=errors)
			std = sqrt(get_diagonal(cov))

			open(unit=59, file=gppredictfile, form='formatted', status='replace')
			write(59, '(3f16.4)') (t(i), mu(i), std(i), i=1,1000)
			close(59)

		end if

	end subroutine dumper


end module nestwrapper
