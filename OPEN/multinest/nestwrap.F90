module nestwrapper

! Nested sampling includes

use Nested
use params
use like
use priors
use array_utils, only: linspace, get_diagonal, sort, swap
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
	   		         nest_fb,nest_MAPfb,nest_liveplot,nest_resume,nest_outfile,nest_initMPI, &
	   		         nest_logZero,nest_maxIter,getLogLike,dumper,live_plot,context)

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
		 						
		! Output arguments
		double precision lnew ! loglikelihood
		integer context ! additional information user wants to pass

		! Local variables
		double precision kmax
		integer i, j
		
		! Transform parameters to physical space using assigned priors
		! Cube(1:nPar) = P, K, ecc, omega, t0  for each planet             |
		!				 jitter at the next to last position	           | - if jitter model				
		!                vsys for each observatory at the last positions   |
		!
		! Cube(1:nPar) = P, K, ecc, omega, t0  for each planet						|
		!				 vsys for each observatory at the next to last positions	| - if GP model				
		!                hyperparameters at the last positions						|

		if (using_jitter) then
			! prior for jitter
			i = nPar-nextra+1 ! index of jitter parameter
			Cube(i) = ModJeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))
			!Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))

			! prior for systematic velocity
			i = npar-nextra+2 ! index of (first) systematic velocity
			Cube(i:) = UniformPrior(Cube(i:), spriorran(i:,1), spriorran(i:,2))			

		else if (using_gp) then
			! prior for systematic velocity
			i = nPar-nextra+1 
			!!!! this needs to be fixed for more than 1 observatory !!!!
			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2))

			! prior for hyperparameters
			i = nPar-nextra+nobserv+1
			do j = i, i+3
				Cube(j) = UniformPrior(Cube(j), spriorran(j,1), spriorran(j,2))
			end do

		else
			! prior for systematic velocity
			i = npar-nextra+1 ! index of (first) systematic velocity
			Cube(i:) = UniformPrior(Cube(i:), spriorran(i:,1), spriorran(i:,2))
		end if


		! prior for period(s)
		do i = 1, nPar-nextra, 5
			! uncomment the following line for Jeffreys prior
			Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))
			! uncomment the following two lines for individual prior (2 planets case)
! 			if (i==1) Cube(i) = GaussianPrior(Cube(i), 0.85359165d0, 5.6d-7)
! 			if (i==6) Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), spriorran(i,2))

! 			if (i==1) Cube(i) = GaussianPrior(Cube(i), 2.2185733d0, 1.9d-6)
		end do

		! priors for ecc, omega, t0
		do i = 3, nPar-nextra, 5
			! uncomment the following line for uniform prior for eccentricity
! 			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2)) ! ecc
			! uncomment the following line for beta prior for eccentricity (based on Kipping)
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
			! set kmax here
! 			kmax = spriorran(i,2)
! 			kmax = spriorran(i,2) * (spriorran(i-1,1) / Cube(i-1))**(1/3.d0) * 1.d0/(sqrt(1-Cube(i+1)**2))
			! uncomment the following line for Modified Jeffreys prior
! 			Cube(i) = ModJeffreysPrior(Cube(i), spriorran(i,1), kmax)
			! uncomment the following line for Jeffreys prior
			!Cube(i) = JeffreysPrior(Cube(i), spriorran(i,1), kmax)
			! uncomment the following line for a uniform prior			
			Cube(i) = UniformPrior(Cube(i), spriorran(i,1), spriorran(i,2))
		end do

		! For e=0, where pericentre is undefined, ω=0 can be chosen
		! such that t0 gives the time of nodal passage
		do i=3,nPar-1,5
			if (Cube(i) == 0) Cube(i+1) = 0.d0
		end do


! 		! this is experimental !!!!!!!!!!!
! 		! and it doesn't seem to solve the problem...
		if (nplanets == 2) then
			if (Cube(1) > Cube(6)) then
				call swap(Cube(1:5), Cube(6:10))
			end if
		end if

! 		print *, Cube
		!call loglike function here 
		call slikelihood(Cube,lnew,context)
! 		stop

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
		print *, maxLogLike

		write(*,*) ' '
		if (nplanets == 1 .and. using_gp) then  ! 1 planet + 4 hyper
			write(*,'(6a13)') (/"    P", "    K", "  ecc", "omega", "   t0", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',nPar - 4,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*4 - 4)
			
			write(*,'(4a13)') (/"t1", "t2", "t3", "t4" /)
			write(fmt,'(a,i2,a)')  '(', 4, 'f13.4)'
			write(*,fmt) paramConstr(nPar*4 - 3:)

		else if (nplanets == 1 .and. using_jitter) then  ! 1 planet + jitter
			write(*,'(7a13)') (/"    P", "    K", "  ecc", "omega", "   t0", "    s", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',nPar,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*4)

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

		else if (nplanets == 2 .and. using_jitter) then  ! 2 planets + jitter
			write(*,'(7a13)') (/"    P", "    K", "  ecc", "omega", "   t0", "    s", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*3+5)
			write(fmt,'(a,i2,a)')  '(',nPar-5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+6:nPar*4)

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

		else if (nplanets == 3 .and. using_jitter) then  ! 3 planets + jitter
			write(*,'(7a13)') (/"    P", "    K", "  ecc", "omega", "   t0", "    s", " vsys" /)
			write(fmt,'(a,i2,a)')  '(',5,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:nPar*3+10)
			write(fmt,'(a,i2,a)')  '(',7,'f13.4)'
			write(*,fmt) paramConstr(nPar*3+11:nPar*4)

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

		else if (nplanets == 0 .and. using_jitter) then  ! jitter + vsys
			write(*,'(2a13)') (/"    s", " vsys" /)
			write(fmt,'(a,i2,a)')  '(', 1+nobserv, 'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:)

		else if (nplanets == 0) then  ! systematic velocity
			write(*,'(a13)') (/" vsys" /)
			write(fmt,'(a,i2,a)')  '(', nobserv, 'f13.4)'
			write(*,fmt) paramConstr(nPar*3+1:)

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


	subroutine live_plot(nPar, nLpt, phyP, l)
		implicit none
		integer nPar !total no. of parameters to save
		integer nLpt !no. of live points
		double precision phyP(nPar,nLpt), l(nLpt)
#ifdef PLOT
		integer, parameter :: nsamples = 3
		double precision u(nsamples)
		integer j(nsamples), jmaxlike(1), min_index, max_index, nt, i

		call random_number(u)
		!print *, u
		min_index=1
		max_index=nLpt
		j = min_index + FLOOR((max_index+1-min_index)*u)
		jmaxlike = maxloc(l)
		!print *, j
		!print *, phyP(:,j)
		nt = size(times)

		!**********************      PLOT      ******************
		CALL PGERAS
		! Define the Viewport
	    CALL PGSVP(0.01, 0.99, 0.01, 0.99)
		! Define the Window
		CALL PGSWIN(REAL(minval(times)-10), REAL(maxval(times)+10), REAL(minval(rvs)-2), REAL(maxval(rvs)+2))
		!CALL PGSWIN(2449464.5956, 2452856.4222, -150.0, 150.0)
		!CALL PGSWIN(13000., 16000., -150.0, 150.0)
		! Draw a box
	    CALL PGSCI(4) ! blue
	    CALL PGBOX ('BCTS', 0., 0, 'BCTSV', 0.0, 0)

	    CALL PGSCI (6) ! magenta
	    CALL PGPT (nt, REAL(times), REAL(rvs), 17)

		do i=2,nsamples
			call get_rvN(times_oversampled, &
						 phyP(1:nPar-1:5, j(i)), & ! periods for all planets
						 phyP(2:nPar-1:5, j(i)), & ! K for all planets
						 phyP(3:nPar-1:5, j(i)), & ! ecc for all planets
						 phyP(4:nPar-1:5, j(i)), & ! omega for all planets
						 phyP(5:nPar-1:5, j(i)), & ! t0 for all planets
						 phyP(nPar, j(i)), & ! systematic velocity
						 vel_oversampled, 5*nt, nplanets)

			!CALL PGSCI (0) ! white
		    !CALL PGLINE(5*nt, REAL(times_oversampled), REAL(last_vel_oversampled(:,i)))
		    !CALL PGPT (nt, REAL(times), REAL(last_vel), 18)
		    CALL PGSCI (1) ! black
		    !CALL PGPT (nt, REAL(times), REAL(vel), 18)
		    CALL PGLINE(5*nt, REAL(times_oversampled), REAL(vel_oversampled))
		    last_vel_oversampled(:, i) = vel_oversampled
		end do

		call get_rvN(times_oversampled, &
					 phyP(1:nPar-1:5, jmaxlike), & ! periods for all planets
					 phyP(2:nPar-1:5, jmaxlike), & ! K for all planets
					 phyP(3:nPar-1:5, jmaxlike), & ! ecc for all planets
					 phyP(4:nPar-1:5, jmaxlike), & ! omega for all planets
					 phyP(5:nPar-1:5, jmaxlike), & ! t0 for all planets
					 phyP(nPar, jmaxlike), & ! systematic velocity
					 vel_oversampled, 5*nt, nplanets)

		CALL PGSCI (3)
		CALL PGSLW (3)  ! line width
	    CALL PGLINE(5*nt, REAL(times_oversampled), REAL(last_vel_oversampled(:,1)))
	    !CALL PGPT (nt, REAL(times), REAL(last_vel), 18)
	    !CALL PGSCI (1) ! black
	    !CALL PGPT (nt, REAL(times), REAL(vel), 18)
	    !CALL PGLINE(5*nt, REAL(times_oversampled), REAL(vel_oversampled))
	    last_vel_oversampled(:, 1) = vel_oversampled
	

	    CALL PGUPDT
		CALL PGEBUF
#endif

	end subroutine live_plot

end module nestwrapper
