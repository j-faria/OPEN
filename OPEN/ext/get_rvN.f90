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
    real (kind=8), dimension(np) :: P, K, ecc, omega, t0
    real (kind=8) :: vsys

! Local variables
    integer :: i
    real(kind=8), parameter :: pi = 3.1415926535897932384626433832795029d0
    real(kind=8), parameter :: twopi = 2.0d0 * pi 

!     vel = vsys
    do i=1,np
      vel = vel + rv_curve(time, P(i), K(i), ecc(i), omega(i), t0(i))
    end do
    vel = vel + vsys
    

contains
  
  elemental real(kind=8) function rv_curve(time, P, K, ecc, omega, t0) result(vel)
  ! RV curve for one planet - see Balan & Lahav (2009), MNRAS, 394, 1936
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in) :: P, K, ecc, omega, t0
    real(kind=8) :: M, f 

    M = twopi * (time-t0)/P
    f = true_anom(M, ecc)
    !vel = K * (sin(f+omega) + ecc*sin(omega))
    vel = K * (cos(omega+f) + ecc*cos(omega))

  end function rv_curve


  elemental real(kind=8) function true_anom(M, ecc) result(nu)
    real(kind=8), intent(in) :: M, ecc
    real(kind=8) :: E, cosE

    E = ecc_anom(M, ecc)

    !! eq 2.6 of Perryman 2011
    !cosE = cos(E)
    !nu = acos((cosE-ecc) / (1._dp-ecc*cosE))

    !! eq 2.7 of Perryman 2011 (seems a tiny bit faster)
    nu = 2.d0 * atan( sqrt((1.d0 + ecc)/(1.d0 - ecc)) * tan(E/2.d0))

    !! old version
    !nu = atan2(sqrt(1._dp-ecc*ecc)*sin(E), cos(E)-ecc)
  
  end function true_anom
  


  elemental real(kind=8) function ecc_anom(M, ecc) result(E)
    real(kind=8), intent(in) :: M, ecc
    real(kind=8), parameter :: derror = 1.0d-7
    integer :: iter
    real(kind=8) :: EA_new, EA_old
    
    ! catch special case
    if (ecc .eq. 0.d0) then
      E = M
      return
    endif

    iter = 0
    EA_old = 0.5d0 * pi 
    if (ecc < 0.8d0) then
      EA_new = newton(EA_old, ecc, M)
      do while (abs(EA_old - EA_new) >= derror .and. iter<200)
        EA_old = EA_new
        EA_new = newton(EA_old, ecc, M)
        iter = iter+1
      end do
    else
      EA_new = strict_iteration(EA_old, ecc, M)
      do while (abs(EA_old - EA_new) >= derror .and. iter<200)
        EA_old = EA_new
        EA_new = strict_iteration(EA_old, ecc, M)
        iter = iter+1
      end do    
    endif  
    E = EA_new

  end function ecc_anom
        
        
  real(kind=8) function halley(EA, e, M) result(EA_new)
    real(kind=8), intent(in) :: EA, e, M
    real(kind=8) :: spsi0, f, fp, fpp

    ! Compute the function and derivatives.
    spsi0 = sin(EA)
    f = EA - e * spsi0 - M
    fp = 1.d0 - e * cos(EA)
    fpp = e * spsi0

    ! Take a second order step.
    EA_new = EA - 2.0d0 * f * fp / (2.d0 * fp * fp - f * fpp)

    ! Deal with looping boundaries properly.
    if (EA_new .gt. twopi) then
      EA_new = 0.5d0 * (EA + twopi)
    elseif (EA_new .lt. 0.0d0) then
      EA_new = 0.5d0 * EA 
    endif

  end function halley
        
  elemental real(kind=8) function strict_iteration(EA, e, M) result(EA_new)
  ! Solve kepler's equation by a simple iteration scheme
    real(kind=8), intent(in) :: EA, e, M
    EA_new = M+e*sin(EA)
  end function strict_iteration
        
  elemental real(kind=8) function newton(EA, e, M) result(EA_new)
  ! The following routine is used to solve kepler's equation using
  ! Newton's method. It is very fast and reliable for small values of
  ! e, but can be wildly erratic for e close to 1. See, e.g, the discussion
  ! in Jean Meeus, Astronomical Algorithms, Willmann-Bell, 1991, 181-193.
    real(kind=8), intent(in) :: EA, e, M
    EA_new = EA + (M + e*sin(EA) - EA)/(1d0 - e*cos(EA))
  end function newton
  
  
  
end subroutine



        
