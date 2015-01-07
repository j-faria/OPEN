module like

	use params
	use gputils, only: gp_n_planet_pars
	use utils1, only: logSumExp
	!use Nested, only: MPI_COMM_WORLD

	!use lib_matrix, only: inverse, determinant
	implicit none
  
	real(kind=8), parameter :: pi = 3.1415926535897932384626433832795029d0
	real(kind=8), parameter :: twopi = 2.d0 * pi 
	real(kind=8), parameter :: lnstwopi = log(sqrt(twopi))



contains  

	subroutine likelihood_init(N)
	! Initialize variables that depend on the data dimensions
		integer, intent(in) :: N ! total size of data-set
		allocate(ss(N), alpha(N), tau(N))
		allocate(observ(N))
		allocate(times(N), rvs(N), errors(N))
		allocate(train_var(N))
		allocate(vel(N), dist(N), sigma(N))
		allocate(vel_oversampled(5*N), times_oversampled(5*N))
		allocate(last_vel_oversampled(5*N, 3))
		allocate(covmat(N,N), inv_covmat(N,N))
	end subroutine likelihood_init

	subroutine likelihood_finish()
	! deallocate memory
		deallocate(ss, alpha, tau)
		deallocate(observ)
		deallocate(times, rvs, errors)
		deallocate(vel, dist, sigma)
		deallocate(last_vel_oversampled, vel_oversampled, times_oversampled)
		deallocate(covmat, inv_covmat)
	end subroutine likelihood_finish		

	subroutine slikelihood(Cube,slhood,context)
	! Likelihood subroutine
	! This is called by getLogLike in nestwrap.f90 which does the parameter rescaling. 
	! Cube(1:nPar) has already physical parameters. 
	! The log-likelihood is returned in slhood
		implicit none
	      
		real(kind=8), intent(in) :: Cube(nest_nPar)
		real(kind=8), intent(out) :: slhood
		integer, intent(in) :: context
		real(kind=8) :: lhood, lhood_test(1,1), det, jitter
		integer :: i, n, ierr, iendpar
		character(len=100) :: fmt

		! times, rvs and errors are defined in params and initialized/read in main
		! Cube(1:nest_nPar) = P, K, ecc, omega, t0, Vsys
		!write(*,*) Cube(3)
 		ss = 1.d0
 		alpha = 1.d0
 		tau = 1.d0

 		if (using_gp) then
			!call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
			!if (i==0) then
			!    write(fmt,'(a,i2,a)')  '(',nest_nPar,'f13.4)'
			!    write(*, fmt) Cube(:)
			!    write(*, fmt) Cube(:gp_n_planet_pars)
			!    write(*, fmt) Cube(gp_n_planet_pars+1:)
			!end if
			!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			!STOP 1
			
			!gp1%gp_kernel%pars = (/1.d0, 25.d0/)
			!gp1%gp_kernel%pars = Cube(gp_n_planet_pars+1:)

			!! hyperparameters
			! amplitude of GP
			gp1%gp_kernel%pars(1) = Cube(gp_n_planet_pars+1)
			! timescale for growth / decay of active regions (d)
			call gp1%gp_kernel%set_kernel_pars(1, (/Cube(gp_n_planet_pars+2)/) )
			!gp1%gp_kernel%kernel1%pars = Cube(gp_n_planet_pars+2)
			! periodicity (recurrence) timescale -> rotation period of the star
			! smoothing parameter (?)
			call gp1%gp_kernel%set_kernel_pars(2, Cube(gp_n_planet_pars+3:))
			!gp1%gp_kernel%kernel2%pars = Cube(gp_n_planet_pars+3:)


			slhood = gp1%get_lnlikelihood(times, rvs, Cube(:gp_n_planet_pars), yerr=errors)
			return
 		end if

 		n = size(times)
 		slhood=-huge(1.d0)*epsilon(1.d0)

 		! get the radial velocity model with these parameters (in vel)
 		if (nplanets > 0) then
 			! this works for one observatory using OPEN/ext/get_rvN.f90
			call get_rvN(times, &
						 Cube(1:nest_nPar-1:5), & ! periods for all planets
						 Cube(2:nest_nPar-1:5), & ! K for all planets
						 Cube(3:nest_nPar-1:5), & ! ecc for all planets
						 Cube(4:nest_nPar-1:5), & ! omega for all planets
						 Cube(5:nest_nPar-1:5), & ! t0 for all planets
						 Cube(nest_nPar), & ! systematic velocity
						 vel, n, nplanets)

			! this works from more than one observatory using OPEN/ext/get_rvN_MultiSite.f90
			! well not just yet
		!	iendpar = nest_nPar - nobserv
		!	call get_rvN(times, &
		!				 Cube(1:iendpar:5), & ! periods for all planets
		!				 Cube(2:iendpar:5), & ! K for all planets
		!				 Cube(3:iendpar:5), & ! ecc for all planets
		!				 Cube(4:iendpar:5), & ! omega for all planets
		!				 Cube(5:iendpar:5), & ! t0 for all planets
		!				 Cube(iendpar+1:), & ! systematic velocity
		!				 observ, & ! array with observatory indices
		!				 vel, n, nplanets, nobserv)

			! this doesn't work...
	    	!call get_rvN(times, &
	    	!           Cube(1), Cube(2), Cube(3), Cube(4), Cube(5), &
	    	!           Cube(6), vel, n, nplanets)
		else
			vel = Cube(nest_nPar)
		end if

	    dist = rvs - vel
	    !call get_covmat(times, errors, observ, ss, alpha, tau, covmat, det, inv_covmat)

	    !lhood_test = -0.5d0 * matmul(matmul(reshape(dist, (/1,n/)), covmat), reshape(dist, (/n,1/)))
	    !lhood = lhood_test(1,1) - 0.5d0*log(twopi**n * det)

		if (using_jitter) then
	    	jitter = Cube(nest_nPar-1)
	    	sigma = errors**2 + jitter**2
	    	!lhood = - n*lnstwopi - sum(log(sqrt(sigma)) + 0.5d0 * dist**2 / sigma)
	    	lhood = - 0.5d0*log(twopi**n * product(sqrt(sigma))) -0.5d0 * sum(dist**2 / sigma)
	    else
			lhood = - 0.5d0*log(twopi**n * product(errors)) -0.5d0 * sum(dist**2 / errors**2)
	    end if

		slhood = logSumExp(slhood,lhood)
		!print *, Cube(1:nest_nPar)
		!stop
		if (doing_debug) write(*,'(f8.3)', advance='no') slhood

	end subroutine slikelihood

!	subroutine get_covmat(times, sigma, observ, ss, alpha, tau, covmat, det, inv_covmat)
!	! Calculates the covariance matrix of the observations and returns its 
!	! determinant and inverse.
!		! vectors with times and uncertainties
!		real(kind=8), intent(in), dimension(:) :: times, sigma
!		! vector with flags for points comming from each observatory
!		integer, intent(in), dimension(:) :: observ
!		! nuisance parameters for each observatory
!		real(kind=8), intent(in), dimension(:) :: ss, alpha, tau
!		! on output, the covariance matrix of the observations and its inverse
!		real(kind=8), intent(out), dimension(:,:) :: covmat, inv_covmat
!		real(kind=8), intent(out) :: det ! determinant of covmat
!
!		! local variables
!		integer :: i, j, nt
!
!		covmat = 0.d0; inv_covmat = 0.d0
!		nt = size(sigma)
!		do i=1,nt
!			do j=1,nt
!				! Kronecker delta on the times
!				if (i==j) covmat(i,j) = (sigma(j)/alpha(j))**2
!				! Kronecker delta on the observatories
!				if (observ(i)==observ(j)) then
!					covmat(i,j) = covmat(i,j) + ss(j)**2 * exp(-abs(times(i)-times(j))/tau(j))
!				endif
!			end do
!		end do	
!		det = determinant(covmat)
!		call inverse(covmat, inv_covmat)
!	end subroutine get_covmat     

end module like
