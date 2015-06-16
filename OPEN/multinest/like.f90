module like

    use params
    use gputils, only: gp_n_planet_pars
    use utils1, only: logSumExp
    use Nested, only: MPI_COMM_WORLD

    !use lib_matrix, only: inverse, determinant
    implicit none
  
    real(kind=8), parameter :: pi = 3.1415926535897932384626433832795029d0
    real(kind=8), parameter :: twopi = 2.d0 * pi 
    real(kind=8), parameter :: lntwopi = log(twopi)
    real(kind=8), parameter :: lnstwopi = log(sqrt(twopi))



contains  

    subroutine likelihood_init(N)
    ! Initialize variables that depend on the data dimensions
        integer, intent(in) :: N ! total size of data-set
        allocate(ss(N), alpha(N), tau(N))
        allocate(observ(N))
        allocate(times(N), rvs(N), errors(N))
        allocate(train_var(N))
        allocate(vel(N), r(N), sigma(N))
        !allocate(vel_oversampled(5*N), times_oversampled(5*N))
        !allocate(last_vel_oversampled(5*N, 3))
        !allocate(covmat(N,N), inv_covmat(N,N))
    end subroutine likelihood_init

    subroutine likelihood_finish()
    ! deallocate memory
        deallocate(ss, alpha, tau)
        deallocate(observ)
        deallocate(times, rvs, errors)
        deallocate(vel, r, sigma)
        !deallocate(last_vel_oversampled, vel_oversampled, times_oversampled)
        !deallocate(covmat, inv_covmat)
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
        real(kind=8) :: lhood, lhood_test(1,1), det, jitter, var
        integer :: i, n, ierr, iendpar
        character(len=100) :: fmt

        ! times, rvs and errors are defined in params and initialized/read in main
!         ss = 0.d0
!         alpha = 1.d0
!         tau = 1.d0



        n = size(times)
        slhood=-huge(1.d0)*epsilon(1.d0)

!         write(*,'(f10.3, i3)') (rvs(i), observ(i), i=1,n)
        iendpar = nest_nPar-nextra
        if (using_jitter) iendpar = iendpar+1
        do i=1,nobserv
            where(observ == i) vel = Cube(iendpar+i)
        end do

        ! get the radial velocity model with these parameters (in vel)
        if (nplanets > 0) then
            call get_rvN(times, &
                         Cube(1:5*nplanets:5), & ! periods for all planets
                         Cube(2:5*nplanets:5), & ! K for all planets
                         Cube(3:5*nplanets:5), & ! ecc for all planets
                         Cube(4:5*nplanets:5), & ! omega for all planets
                         Cube(5:5*nplanets:5), & ! t0 for all planets
                         0.d0, & ! systematic velocity
                         vel, n, nplanets)

            ! this works from more than one observatory using OPEN/ext/get_rvN_MultiSite.f90
            ! well not just yet  
!             call get_rvN(times, &
!                          Cube(1:iendpar:5), & ! periods for all planets
!                          Cube(2:iendpar:5), & ! K for all planets
!                          Cube(3:iendpar:5), & ! ecc for all planets
!                          Cube(4:iendpar:5), & ! omega for all planets
!                          Cube(5:iendpar:5), & ! t0 for all planets
!                          Cube(iendpar+1:), & ! systematic velocity(ies)
!                          observ, & ! array with observatory indices
!                          vel, n, nplanets, nobserv)

        end if

        ! residuals: what is left when the planets (or just vsys) is subtracted from the data
        r = rvs - vel  
        
!             call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
!             if (i==0) then
!                 print *, 'r = ', r
!                write(fmt,'(a,i2,a)')  '(',nest_nPar,'f13.4)'
!                write(*, fmt) Cube(:)
! !                write(*, fmt) Cube(:gp_n_planet_pars)
! !                write(*, fmt) Cube(gp_n_planet_pars)
! !                write(*, fmt) Cube(gp_n_planet_pars+2)
! !                write(*, fmt) Cube(gp_n_planet_pars+3)
! !                write(*, fmt) Cube(gp_n_planet_pars+4:gp_n_planet_pars+5)
!             end if
!             call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!             STOP 1


        if (using_gp) then
            ! gp1 is an object of type GP (class is defined in gp_utils.f90)
            ! its hierarchy is quite complex:
            !   gp1%gp_kernel%pars    holds the parameters of the kernel's sum, i.e,
            !                         the multiplying constants behind the 2 kernels
            !    
            !   gp1%sub_kernel1       is the first sub kernel (a product of kernels)
            !   gp1%sub_kernel1%pars  holds the parameters of the kernels' product, i.e.,
            !                         the multiplying constants - do not change these
            !   parameters of the two kernels being multiplied can be accessed and changed by
            !   gp1%sub_kernel1%get_kernel_pars(k)
            !   gp1%sub_kernel1%set_kernel_pars(k, new_pars)
            !   where k=1 corresponds to the ExpSquared kernel (with 1 parameter)
            !         k=2 corresponds to the ExpSineSquared kernel (with 2 parameters)
            !
            !   gp%sub_kernel2        is the second sub kernel (a simple DiagonalKernel)
            !   gp1%sub_kernel1%pars  holds the DiagonalKernel parameter (the jitter value)

            !! update the hyperparameters

            ! amplitude of the quasi-periodic component
            gp1%gp_kernel%pars(1) = Cube(gp_n_planet_pars+1)
            ! jitter
            gp1%sub_kernel2%pars(1) = Cube(gp_n_planet_pars+2)
            ! timescale for growth / decay of active regions (d)
            call gp1%sub_kernel1%set_kernel_pars(1, (/Cube(gp_n_planet_pars+3)/) )
            ! periodicity (recurrence) timescale -> rotation period of the star
            ! correlation scale
            call gp1%sub_kernel1%set_kernel_pars(2, Cube(gp_n_planet_pars+4:))

            !! log likelihood of the residuals 
            ! mean_fun is 0. because vsys is already taken care of
            lhood = gp1%get_lnlikelihood(times, r, (/0.d0/), yerr=errors)


        else if (using_jitter) then
            jitter = Cube(nest_nPar-nextra+1)
            do i=1,n
                var = errors(i)*errors(i) + jitter*jitter
                lhood = lhood - 0.5*lntwopi - 0.5*log(var) - 0.5*(r(i)**2)/var
            end do

!             sigma(:) = errors(:)**2 + jitter**2
            !lhood = - n*lnstwopi - sum(log(sqrt(sigma)) + 0.5d0 * r**2 / sigma)
!             lhood = -0.5d0*n*lntwopi - n*sum(log(sqrt(sigma))) -0.5d0*sum(r**2 / sigma)
!             print *, lhood
            ! lhood = - 0.5d0*log(twopi**n * product(sqrt(sigma))) -0.5d0 * sum(r**2 / sigma)
        else
            do i=1,n
                var = errors(i)*errors(i)
                lhood = lhood - 0.5*lntwopi - 0.5*log(var) - 0.5*(r(i)**2)/var
            end do
            ! lhood = -0.5d0*n*lntwopi - n*sum(log(errors)) -0.5d0*sum(r**2 / errors**2)
            ! lhood = - 0.5d0*log(twopi**n * product(errors)) -0.5d0 * sum(r**2 / errors**2)
        end if

        slhood = logSumExp(slhood,lhood)
!         stop
        if (doing_debug) write(*,'(f8.3)', advance='no') slhood

    end subroutine slikelihood

!    subroutine get_covmat(times, sigma, observ, ss, alpha, tau, covmat, det, inv_covmat)
!    ! Calculates the covariance matrix of the observations and returns its 
!    ! determinant and inverse.
!        ! vectors with times and uncertainties
!        real(kind=8), intent(in), dimension(:) :: times, sigma
!        ! vector with flags for points comming from each observatory
!        integer, intent(in), dimension(:) :: observ
!        ! nuisance parameters for each observatory
!        real(kind=8), intent(in), dimension(:) :: ss, alpha, tau
!        ! on output, the covariance matrix of the observations and its inverse
!        real(kind=8), intent(out), dimension(:,:) :: covmat, inv_covmat
!        real(kind=8), intent(out) :: det ! determinant of covmat
!
!        ! local variables
!        integer :: i, j, nt
!
!        covmat = 0.d0; inv_covmat = 0.d0
!        nt = size(sigma)
!        do i=1,nt
!            do j=1,nt
!                ! Kronecker delta on the times
!                if (i==j) covmat(i,j) = (sigma(j)/alpha(j))**2
!                ! Kronecker delta on the observatories
!                if (observ(i)==observ(j)) then
!                    covmat(i,j) = covmat(i,j) + ss(j)**2 * exp(-abs(times(i)-times(j))/tau(j))
!                endif
!            end do
!        end do    
!        det = determinant(covmat)
!        call inverse(covmat, inv_covmat)
!    end subroutine get_covmat     

end module like
