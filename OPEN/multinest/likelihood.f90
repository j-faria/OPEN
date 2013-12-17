!*  Purpose
!*  =======
!*  Likelihood calculation
subroutine likelihood(time, datum, sigma, &
                      P, K, ecc, omega, t0, vsys, vel, lhood, nt, np)
    implicit none

!f2py intent(in) time(nt)
!f2py intent(in) datum(nt)
!f2py intent(in) sigma(nt)
!f2py intent(inout) vel(nt)
!f2py intent(in) P(np)
!f2py intent(in) K(np)
!f2py intent(in) ecc(np)
!f2py intent(in) omega(np)
!f2py intent(in) t0(np)
!f2py intent(in) vsys
!f2py intent(hide) lhood
!f2py intent(hide) nt
!f2py intent(hide) np

! Input arguments
    integer nt, np ! number of observations, planets
    real (kind=8), dimension(nt) :: time, datum, sigma, vel
    real (kind=8), dimension(np) :: P, K, ecc, omega, t0
    real (kind=8) :: vsys, lhood

! Local variables
    integer,parameter :: sp = selected_real_kind(p=6,r=37)
    integer,parameter :: dp = selected_real_kind(p=15,r=307)

    real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
    real(dp), parameter :: twopi = 2.0_dp * pi 

    real(dp) :: lhd     

    ! get the radial velocity model with these parameters (in vel)
    call get_rvN(time, P, K, ecc, omega, t0, vsys, vel, nt, np)

!    lhood = product((1._dp / sqrt(twopi*sigma**2)) * exp(-0.5_dp * ((datum - vel)/sigma)**2))
!     write(*,*) lhood


end subroutine