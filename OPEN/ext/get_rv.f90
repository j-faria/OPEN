!*  Purpose
!*  =======
!*  Calculate the radial velocity curve at times 'time' with given parameters P,
!*  K, ecc, omega and t0 corresponding to period, semi-amplitude, eccentricity,
!*  argument of the periastron and time of periastron.
subroutine get_rv( time, P, K, ecc, omega, t0, vel, nt)
    implicit none
    
!f2py intent(in) time(nt)
!f2py intent(in) vel(nt)
!f2py intent(in) P
!f2py intent(in) K
!f2py intent(in) ecc
!f2py intent(in) omega
!f2py intent(in) t0
!f2py intent(hide) nt

! Input arguments
    integer nt
    real (kind=8) time(nt), vel(nt)
    real (kind=8) p, k, ecc, omega, t0

! Local variables
    integer,parameter :: sp = selected_real_kind(p=6,r=37)
    integer,parameter :: dp = selected_real_kind(p=15,r=307)

    real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
    real(dp), parameter :: twopi = 2.0_dp * pi 

    vel = rv_curve(time, P, K, ecc, omega, t0)
    

contains
  
  elemental real(dp) function rv_curve(time, P, K, ecc, omega, t0) result(vel)
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
!    write(*,*) 'M=', M, 'ecc=', ecc 
    E = ecc_anom(M, ecc)
!    write(*,*) 'E=', E
    nu = atan2(sqrt(1._dp-ecc*ecc)*sin(E),cos(E)-ecc)
  
  end function true_anom
  


  elemental real(dp) function ecc_anom(M, ecc) result(E)
    real(dp), intent(in) :: M, ecc
    real(dp), parameter :: derror = 1.0d-3
    real(dp) :: EA_new, EA_old
    
!        write(*,*) 'M=', M, 'ecc=', ecc 

    EA_old = 0.5_dp * pi 
    EA_new = newton(EA_old, ecc, M)
    do while (abs(EA_old - EA_new) >= derror)
      EA_old = EA_new
      EA_new = newton(EA_old, ecc, M)
    end do
    
!    write(*,*) 'E=',EA_new
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
        
        
  elemental real(dp) function newton(EA, e, M) result(EA_new)
  ! The following routine is used to solve kepler's equation using
  ! Newton's method. It is very fast and reliable for small values of
  ! e, but can be wildly erratic for e close to 1. See, e.g, the discussion
  ! in Jean Meeus, Astronomical Algorithms, Willmann-Bell, 1991, 181-193.
    real(dp), intent(in) :: EA, e, M
    EA_new = EA + (M + e*sin(EA) - EA)/(1 - e*cos(EA))
  end function newton
  
  
  
end subroutine



        
