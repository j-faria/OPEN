!*  Purpose
!*  =======
!*  Calculate the radial velocity curve at times 'time' with given parameters P,
!*  K, ecc, omega and t0 corresponding to period, semi-amplitude, eccentricity,
!*  argument of the periastron and time of periastron. vsys is the systematic
!*  velocity.
!*
!*  This routine works for NP planets but one observatory only!
subroutine get_rvN(time, P, K, ecc, omega, t0, vsys, vel, nt, np)
    implicit none
    
!f2py intent(in) time(nt)
!f2py intent(inout) vel(nt)
!f2py intent(in) P(np)
!f2py intent(in) K(np)
!f2py intent(in) ecc(np)
!f2py intent(in) omega(np)
!f2py intent(in) t0(np)
!f2py intent(in) vsys
!f2py intent(hide) nt
!f2py intent(hide) np

! Input arguments
    integer nt, np ! number of observations, planets
    real (kind=8), dimension(nt) :: time, vel
    real (kind=8), dimension(np) :: p, k, ecc, omega, t0
    real (kind=8) :: vsys

! Local variables
    integer,parameter :: sp = selected_real_kind(p=6,r=37)
    integer,parameter :: dp = selected_real_kind(p=15,r=307)
    integer :: i

    real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
    real(dp), parameter :: twopi = 2.0_dp * pi 

    vel = vsys
    do i=1,np
      vel = vel - rv_curve(time, P(i), K(i), ecc(i), omega(i), t0(i))
    end do 
    

contains
  
  elemental real(dp) function rv_curve(time, P, K, ecc, omega, t0) result(vel)
  ! RV curve for one planet - see Balan & Lahav (2009), MNRAS, 394, 1936
    real(dp), intent(in) :: time
    real(dp), intent(in) :: P, K, ecc, omega, t0
    real(dp) :: M, f 

    M = twopi * (time-t0)/P
    f = true_anom(M, ecc)
    vel = K * (sin(f+omega) + ecc*sin(omega))

  end function rv_curve


  elemental real(dp) function true_anom(M, ecc) result(nu)
    real(dp), intent(in) :: M, ecc
    real(dp) :: E

    E = ecc_anom(M, ecc)
    nu = atan2(sqrt(1._dp-ecc*ecc)*sin(E),cos(E)-ecc)
  
  end function true_anom
  


  elemental real(dp) function ecc_anom(M, ecc) result(E)
    real(dp), intent(in) :: M, ecc
    real(dp), parameter :: derror = 1.0d-3
    real(dp) :: EA_new, EA_old
    
    EA_old = 0.5_dp * pi 
    if (ecc < 0.8d0) then
      EA_new = newton(EA_old, ecc, M)
      do while (abs(EA_old - EA_new) >= derror)
        EA_old = EA_new
        EA_new = newton(EA_old, ecc, M)
      end do
    else
      EA_new = strict_iteration(EA_old, ecc, M)
      do while (abs(EA_old - EA_new) >= derror)
        EA_old = EA_new
        EA_new = strict_iteration(EA_old, ecc, M)
      end do    
    endif  
    E = EA_new

  end function ecc_anom
        
        
  real(dp) function halley(EA, e, M) result(EA_new)
    real(dp), intent(in) :: EA, e, M
    real(dp) :: spsi0, f, fp, fpp

    ! Compute the function and derivatives.
    spsi0 = sin(EA)
    f = EA - e * spsi0 - M
    fp = 1._dp - e * cos(EA)
    fpp = e * spsi0

    ! Take a second order step.
    EA_new = EA - 2.0d0 * f * fp / (2.d0 * fp * fp - f * fpp)

    ! Deal with looping boundaries properly.
    if (EA_new .gt. twopi) then
      EA_new = 0.5 * (EA + twopi)
    elseif (EA_new .lt. 0.0) then
      EA_new = 0.5 * EA 
    endif

  end function halley
        
  elemental real(dp) function strict_iteration(EA, e, M) result(EA_new)
  ! Solve kepler's equation by a simple iteration scheme
    real(dp), intent(in) :: EA, e, M
    EA_new = M+e*sin(EA)
  end function strict_iteration
        
  elemental real(dp) function newton(EA, e, M) result(EA_new)
  ! The following routine is used to solve kepler's equation using
  ! Newton's method. It is very fast and reliable for small values of
  ! e, but can be wildly erratic for e close to 1. See, e.g, the discussion
  ! in Jean Meeus, Astronomical Algorithms, Willmann-Bell, 1991, 181-193.
    real(dp), intent(in) :: EA, e, M
    EA_new = EA + (M + e*sin(EA) - EA)/(1 - e*cos(EA))
  end function newton
  
  
  
end subroutine



        
