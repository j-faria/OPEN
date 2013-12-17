module like

	use params
	use utils1, only: logSumExp

	use lib_matrix, only: inverse, determinant
	implicit none
  
	real(kind=8), parameter :: pi = 3.1415926535897932384626433832795029d0
	real(kind=8), parameter :: twopi = 2.d0 * pi 


contains  

	subroutine likelihood_init(N)
	! Initialize variables that depend on the data dimensions
		integer, intent(in) :: N ! total size of data-set
		allocate(ss(N), alpha(N), tau(N))
		allocate(observ(N))
		allocate(times(N), rvs(N), errors(N))
		allocate(vel(N), dist(N))
		allocate(covmat(N,N), inv_covmat(N,N))

	end subroutine likelihood_init

	subroutine slikelihood(Cube,slhood)
	! Likelihood subroutine
	! This is called by getLogLike in nestwrap.f90 which does the 
	! parameter rescaling. 
	! Cube(1:nPar) has already physical parameters. 
	! The log-likelihood is returned in slhood
		implicit none
	      
		real(kind=8), intent(in) :: Cube(nest_nPar)
		real(kind=8), intent(out) :: slhood
		real(kind=8) :: lhood, lhood_test(1,1), det
		integer :: i, n

		! times, rvs and errors are defined in params and initialized/read in main
		! Cube(1:nest_nPar) = P, K, ecc, omega, t0, Vsys
 		!write(*,*) Cube(1)
 		observ = 1
 		ss = 1.d0
 		alpha = 1.d0
 		tau = 1.d0

 		n = size(times)
 		slhood=-huge(1.d0)*epsilon(1.d0)

 		! get the radial velocity model with these parameters (in vel)
    	call get_rvN(times, &
    	           Cube(1), Cube(2), Cube(3), Cube(4), Cube(5), &
    	           Cube(6), vel, n, 1)

	    dist = rvs - vel
	    call get_covmat(times, errors, observ, ss, alpha, tau, covmat, det, inv_covmat)

	    lhood_test = -0.5d0 * matmul(matmul(reshape(dist, (/1,n/)), covmat), reshape(dist, (/n,1/)))
	    lhood = lhood_test(1,1) - 0.5d0*log(twopi**n * det)

! 		call likelihood(times, rvs, errors, &
!                       Cube(1), Cube(2), Cube(3), Cube(4), Cube(5), &
!                       0.d0, vel, lhood, 119, 1)

		
		slhood=logSumExp(slhood,lhood)
		if (doing_debug) write(*,'(f8.3)', advance='no') slhood

	end subroutine slikelihood

	subroutine get_covmat(times, sigma, observ, ss, alpha, tau, covmat, det, inv_covmat)
	! Calculates the covariance matrix of the observations and returns its 
	! determinant and inverse.
		! vectors with times and uncertainties
		real(kind=8), intent(in), dimension(:) :: times, sigma
		! vector with flags for points comming from each observatory
		integer, intent(in), dimension(:) :: observ
		! nuisance parameters for each observatory
		real(kind=8), intent(in), dimension(:) :: ss, alpha, tau
		! on output, the covariance matrix of the observations and its inverse
		real(kind=8), intent(out), dimension(:,:) :: covmat, inv_covmat
		real(kind=8), intent(out) :: det ! determinant of covmat

		! local variables
		integer :: i, j, nt

		covmat = 0.d0; inv_covmat = 0.d0
		nt = size(sigma)
		do i=1,nt
			do j=1,nt
				! Kronecker delta on the times
				if (i==j) covmat(i,j) = (sigma(j)/alpha(j))**2
				! Kronecker delta on the observatories
				if (observ(i)==observ(j)) then
					covmat(i,j) = covmat(i,j) + ss(j)**2 * exp(-abs(times(i)-times(j))/tau(j))
				endif
			end do
		end do	
		det = determinant(covmat)
		call inverse(covmat, inv_covmat)
	end subroutine get_covmat     

end module like
