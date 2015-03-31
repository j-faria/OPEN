subroutine gp_predictor(times, rvs, errors, pred, planetpar, hyperpar, meanf, nplanets, nt)
    use gputils
    implicit none

!f2py intent(in) times(nt)
!f2py intent(in) rvs(nt)
!f2py intent(in) errors(nt)
!f2py intent(out) pred(nt)
!f2py intent(hide) nplanets
!f2py intent(hide) nt

	!! arguments
    integer :: nt
    real(kind=8), intent(in), dimension(nt) :: times, rvs, errors
    real(kind=8), intent(out), dimension(nt) :: pred
    integer, intent(in) :: nplanets
    real(kind=8), intent(in), dimension(5*nplanets+1) :: planetpar
    real(kind=8), intent(in), dimension(5) :: hyperpar
    character(len=*), intent(in) :: meanf


    !! Gaussian process variables
    type(DiagonalKernel), target :: k3
    type(ExpSquaredKernel), target :: k5
    type(ExpSineSquaredKernel), target :: k6
    type(SumKernels), target :: k7
    type(ProductKernels), target :: k8
    type(GP) gp1
    class(Kernel), pointer :: kernel_to_pass, subk1, subk2

    !! local variables
    integer :: nPar
    real(kind=8), dimension(nt, nt) :: cov
    real(kind=8), dimension(nt) :: vel, r

    nPar = size(planetpar)

    k3 = DiagonalKernel((/3d-7/))
    k5 = ExpSquaredKernel((/ sqrt(10.d0) /))
    k6 = ExpSineSquaredKernel((/ 30.d0, sqrt(2.d0) /))

    k8 = ProductKernels((/1.d0, 1.d0 /))
    k8%kernel1 => k5
    k8%kernel2 => k6

    k7 = SumKernels((/10d0, 1d0/))
    k7%kernel1 => k8 ! quasi-periodic
    k7%kernel2 => k3 ! jitter

    kernel_to_pass => k7
    subk1 => k8
    subk2 => k3

    gp1 = GP(k7%evaluate_kernel(times, times), kernel_to_pass, subk1, subk2)
    gp1%mean_fun => mean_fun_constant

    gp_n_planets = nplanets
    gp_n_planet_pars = 5*nplanets + 1

    gp1%gp_kernel%pars(1) = hyperpar(1)
    gp1%sub_kernel2%pars(1) = hyperpar(2)
    call gp1%sub_kernel1%set_kernel_pars(1, (/ hyperpar(3) /) )
    call gp1%sub_kernel1%set_kernel_pars(2, hyperpar(4:))

!     print *, gp1%gp_kernel%pars
!     print *, gp1%sub_kernel2%pars
!     call gp1%sub_kernel1%get_kernel_pars(1)
!     call gp1%sub_kernel1%get_kernel_pars(2)

    ! fill pred with systematic velocities
!     do i=1,nobserv
!         where(observ == i) pred = planetpar(size(planetpar)-)
!     end do
    vel = planetpar(nPar)

    ! get the radial velocity model with these parameters (in vel)
    if (gp_n_planets > 0) then
        call get_rvN(times, &
                     planetpar(1:nPar-1:5), & ! periods for all planets
                     planetpar(2:nPar-1:5), & ! K for all planets
                     planetpar(3:nPar-1:5), & ! ecc for all planets
                     planetpar(4:nPar-1:5), & ! omega for all planets
                     planetpar(5:nPar-1:5), & ! t0 for all planets
                     0.d0, & ! systematic velocity
                     vel, nt, gp_n_planets)
    end if

    ! residuals: what is left when the planets (or just vsys) is subtracted from the data
    r = rvs - vel

    ! calculate predictions at observed times
    call gp1%predict(times, r, (/0.d0/), times, pred, cov, yerr=errors)
    pred = pred + vel


    return

end subroutine gp_predictor
