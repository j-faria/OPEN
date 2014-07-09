module nestwrapper

! Nested sampling includes

use Nested
use params
use like
use priors
   
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
		

		! Transform parameters to physical space using assigned priors
		! Cube(1:nPar) = P, K, ecc, omega, t0  for each planet
		!                vsys at the last position 
		if (nest_context == 11) then  ! 1 planet
			P = UniformPrior(Cube(1), spriorran(1,1), spriorran(1,2))
			ecc = UniformPrior(Cube(3), spriorran(3,1), spriorran(3,2))
			omega = UniformPrior(Cube(4), spriorran(4,1), spriorran(4,2))
			t0 = UniformPrior(Cube(5), spriorran(5,1), spriorran(5,2))
			K = UniformPrior(Cube(2), spriorran(2,1), spriorran(2,2)) ! for now

			Vsys = UniformPrior(Cube(6), spriorran(6,1), spriorran(6,2))

			Cube(1:nPar) = (/P, K, ecc, omega, t0, Vsys/)

		else if (nest_context == 21) then
			P1 = UniformPrior(Cube(1), spriorran(1,1), spriorran(1,2))
			ecc1 = UniformPrior(Cube(3), spriorran(3,1), spriorran(3,2))
			omega1 = UniformPrior(Cube(4), spriorran(4,1), spriorran(4,2))
			t01 = UniformPrior(Cube(5), spriorran(5,1), spriorran(5,2))
			K1 = UniformPrior(Cube(2), spriorran(2,1), spriorran(2,2)) ! for now !!!!!

			P2 = UniformPrior(Cube(6), spriorran(6,1), spriorran(6,2))
			ecc2 = UniformPrior(Cube(8), spriorran(8,1), spriorran(8,2))
			omega2 = UniformPrior(Cube(9), spriorran(9,1), spriorran(9,2))
			t02 = UniformPrior(Cube(10), spriorran(10,1), spriorran(10,2))
			K2 = UniformPrior(Cube(7), spriorran(7,1), spriorran(7,2)) ! for now !!!!!		

			Vsys = UniformPrior(Cube(11), spriorran(11,1), spriorran(11,2))

			Cube(1:nPar) = (/P1, K1, ecc1, omega1, t01, &
			                 P2, K2, ecc2, omega2, t02, &
			                 Vsys/)

		end if
		!write(unit=*, fmt=*) '+++', Cube(1:nPar)

		!call loglike function here 
		call slikelihood(Cube,lnew,context)

	end subroutine getLogLike


	subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, INSlogZ, logZerr, context)
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
		
		! now do something
		!if (doing_debug) write(*,*) paramConstr(:)
		write(*,*) paramConstr(1:nPar)

	end subroutine dumper

	! these should be in priors.f90 but I'm afraid of recompiling it :)
	!=======================================================================
	! Uniform[0:1]  ->  Exponential[0:inf]
	function ExpPrior(r)
      	implicit none
      	double precision r,ExpPrior
      	
      	ExpPrior=-log(1.d0-r)

	end function ExpPrior

	!=======================================================================
	! Uniform[0:1]  ->  Jeffreys[x1:x2]
	function JeffreysPrior(r, x1, x2)
		implicit none
		double precision r,x1,x2,JeffreysPrior

		JeffreysPrior=x1*(x2/x1)**r
		
	end function JeffreysPrior


end module nestwrapper
