subroutine gp_predictor(times, rvs, errors, pred, planetpar, hyperpar, nplanets, nt)
    use gputils
    implicit none

!f2py intent(in) times(nt)
!f2py intent(in) rvs(nt)
!f2py intent(in) errors(nt)
!f2py intent(out) pred(nt)
!f2py intent(in) planetpar(nplanets)
!f2py intent(hide) nplanets
!f2py intent(hide) nt

	!! arguments
    integer :: nt
    real(kind=8), intent(in), dimension(nt) :: times, rvs, errors
    real(kind=8), intent(out), dimension(nt) :: pred
    integer, intent(in) :: nplanets
    real(kind=8), intent(in), dimension(5*nplanets+1) :: planetpar
    real(kind=8), intent(in), dimension(4) :: hyperpar

    !! Gaussian process variables
    !type(WhiteKernel), target :: k2
    !type(DiagonalKernel), target :: k3
    !type(ExpKernel), target :: k4
    type(ExpSquaredKernel), target :: k5
    type(ExpSineSquaredKernel), target :: k6
    !type(SumKernels), target :: k7
    type(ProductKernels), target :: k8
    type(GP) gp1
    class(Kernel), pointer :: kernel_to_pass

    !! local variables
    real(kind=8), dimension(nt, nt) :: cov

    k5 = ExpSquaredKernel((/ 1.d0 /))
    k6 = ExpSineSquaredKernel((/ 1.d0, 25.d0 /))
    k8 = ProductKernels((/1.d0, 1.d0 /))
    k8%kernel1 => k5
    k8%kernel2 => k6
    kernel_to_pass => k8

    gp1 = GP(k8%evaluate_kernel(times, times), kernel_to_pass)
    gp1%mean_fun => mean_fun_keplerian

    gp_n_planets = nplanets
    gp_n_planet_pars = 5*nplanets + 1

    gp1%gp_kernel%pars(1) = hyperpar(1)
    call gp1%gp_kernel%set_kernel_pars(1, (/ hyperpar(2) /) )
    call gp1%gp_kernel%set_kernel_pars(2, hyperpar(3:))

    call gp1%predict(times, rvs, planetpar, times, pred, cov, yerr=errors)

    return

end subroutine gp_predictor
