SUBROUTINE GLOMBSCARGLE(T, Y, SIG, W, NCPU, P, M, NT, NW)
    !f2py threadsafe
    !f2py intent(in) T(NT)
    !f2py intent(in) Y(NT)
    !f2py intent(in) SIG(NT)
    !f2py intent(in) W(NW)
    !f2py integer intent(in), optional :: NCPU = -1
    !f2py intent(out) P(NW)
    !f2py intent(out) M
    !f2py intent(hide) NT
    !f2py intent(hide) NW
    use omp_lib
    IMPLICIT NONE

    !*  Input arguments
    INTEGER NT, NW, M
    INTEGER,optional :: NCPU
    REAL (KIND=8) T(NT), Y(NT), SIG(NT)
    REAL (KIND=8) W(NW), P(NW)

    !*  Calculate the Generalized Lomb-Scargle (GLS) periodogram as defined by
    !*  Zechmeister, M. & Kürster, M., A&A 496, 577-584, 2009 (ZK2009) 
    !*
    !*  Arguments
    !*  =========
    !* 
    !*  T   (input) REAL (KIND=8) array, dimension (NT)
    !*      Sample times
    !* 
    !*  Y   (input) REAL (KIND=8) array, dimension (NT)
    !*      Measurement values
    !* 
    !* SIG  (input) REAL (KIND=8) array, dimension (NT)
    !*      Measured uncertainties
    !*
    !*  W   (input) REAL (KIND=8) array, dimension (NW)
    !*      Angular frequencies for output periodogram
    !* 
    !* NCPU (input) INTEGER, OPTIONAL
    !*      Number of CPUs to distribute calculation
    !*
    !*  P   (output) REAL (KIND=8) array, dimension (NW)
    !*      Lomb-Scargle periodogram
    !* 
    !*  NT (input) integer
    !*      Dimension of input arrays
    !* 
    !*  NW (output) integer
    !*      Dimension of frequency and output arrays

    !*  Local variables
    INTEGER I, J
    REAL (KIND=8) sig2(NT), ww(NT), th(NT), yh(NT)
    REAL (KIND=8) upow(NW)!, a(NW), b(NW), off(NW)
    REAL (KIND=8) pi, pi2, pi4
    REAL (KIND=8) Y_, YY, C, S, YC, YS, CCh, CSh, SSh, CC, SS, CS, D
    REAL (KIND=8) x(NT), cosx(NT), sinx(NT), wcosx(NT), wsinx(NT)

    pi2 = 6.2831853071795862d0
    pi4 = 12.566370614359172d0
          
    th = T

    sig2 = sig*sig
    ww = (1.d0 / sum(1.d0/sig2)) / sig2  ! w of ZK2009, the normalized weights

    Y_ = sum(ww*Y)              ! Eq. (7)
    yh = Y - Y_                 ! Subtract weighted mean

    YY = sum(ww * yh**2)        ! Eq. (16)

    ! Unnormalized power
    upow = 0.d0
    !a = 0.d0
    !b = 0.d0
    !off = 0.d0

    ! If NCPU is present, distribute the calculation using OpenMP
    if (present(NCPU).and.(.not.(NCPU.eq.-1))) then
        call OMP_SET_NUM_THREADS(NCPU)
    else
        call OMP_SET_NUM_THREADS(1)
    endif

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(x,cosx,sinx,wcosx,wsinx,C,S,YC,YS,CCh,CSh,SSh,CC,SS,CS,D)
    do I=1,NW

        x = W(I) * th
        cosx = cos(x)
        sinx = sin(x)
        wcosx = ww*cosx         ! attach weights
        wsinx = ww*sinx         ! attach weights

        C = sum(wcosx)         ! Eq. (8)
        S = sum(wsinx)         ! Eq. (9)

        YC = sum(yh*wcosx)     ! Eq. (17)
        YS = sum(yh*wsinx)     ! Eq. (18)
        CCh = sum(wcosx*cosx)  ! Eq. (13)
        CSh = sum(wcosx*sinx)  ! Eq. (15)
        SSh = 1.d0 - CCh
        CC = CCh - C*C         ! Eq. (13)
        SS = SSh - S*S         ! Eq. (14)
        CS = CSh - C*S         ! Eq. (15)
        D = CC*SS - CS*CS      ! Eq. (6)

        !a(I) = (YC*SS-YS*CS) / D
        !b(I) = (YS*CC-YC*CS) / D
        !off(I) = -a(I)*C - b(I)*S

        upow(I) = (SS*YC*YC + CC*YS*YS - 2.d0*CS*YC*YS) / (YY*D) ! Eq. (5)

    end do
    !$OMP END PARALLEL DO

    ! An ad-hoc estimate of the number of independent frequencies 
    ! see discussion following Eq. (24)
    M = (maxval(W/pi2) - minval(W/pi2)) * (maxval(th) - minval(th))

    P = upow
        
END SUBROUTINE GLOMBSCARGLE


SUBROUTINE GLOMBSCARGLE_EXTRA_OUT(T, Y, SIG, W, NCPU, P, M, A, B, OFFSET, NT, NW)
    !f2py threadsafe
    !f2py intent(in) T(NT)
    !f2py intent(in) Y(NT)
    !f2py intent(in) SIG(NT)
    !f2py intent(in) W(NW)
    !f2py integer intent(in), optional :: NCPU = -1
    !f2py intent(out) P(NW)
    !f2py intent(out) M
    !f2py intent(out) a(NW)
    !f2py intent(out) b(NW)
    !f2py intent(out) offset(NW)
    !f2py intent(hide) NT
    !f2py intent(hide) NW
    use omp_lib
    IMPLICIT NONE

    !*  Input arguments
    INTEGER NT, NW, M
    INTEGER,optional :: NCPU
    REAL (KIND=8) T(NT), Y(NT), SIG(NT)
    REAL (KIND=8) W(NW), P(NW)
    REAL (KIND=8) a(NW), b(NW), offset(NW)

    !*  Calculate the Generalized Lomb-Scargle (GLS) periodogram as defined by
    !*  Zechmeister, M. & Kürster, M., A&A 496, 577-584, 2009 (ZK2009) 
    !*
    !*  Arguments
    !*  =========
    !* 
    !*  T   (input) REAL (KIND=8) array, dimension (NT)
    !*      Sample times
    !* 
    !*  Y   (input) REAL (KIND=8) array, dimension (NT)
    !*      Measurement values
    !* 
    !* SIG  (input) REAL (KIND=8) array, dimension (NT)
    !*      Measured uncertainties
    !*
    !*  W   (input) REAL (KIND=8) array, dimension (NW)
    !*      Angular frequencies for output periodogram
    !* 
    !* NCPU (input) INTEGER, OPTIONAL
    !*      Number of CPUs to distribute calculation
    !*
    !*  P   (output) REAL (KIND=8) array, dimension (NW)
    !*      Lomb-Scargle periodogram
    !* 
    !*  NT (input) integer
    !*      Dimension of input arrays
    !* 
    !*  NW (output) integer
    !*      Dimension of frequency and output arrays

    !*  Local variables
    INTEGER I, J
    REAL (KIND=8) sig2(NT), ww(NT), th(NT), yh(NT)
    REAL (KIND=8) upow(NW)
    REAL (KIND=8) pi, pi2, pi4
    REAL (KIND=8) Y_, YY, C, S, YC, YS, CCh, CSh, SSh, CC, SS, CS, D
    REAL (KIND=8) x(NT), cosx(NT), sinx(NT), wcosx(NT), wsinx(NT)

    pi2 = 6.2831853071795862d0
    pi4 = 12.566370614359172d0
          
    th = T

    sig2 = sig*sig
    ww = (1.d0 / sum(1.d0/sig2)) / sig2  ! w of ZK2009, the normalized weights

    Y_ = sum(ww*Y)              ! Eq. (7)
    yh = Y - Y_                 ! Subtract weighted mean

    YY = sum(ww * yh**2)        ! Eq. (16)

    ! Unnormalized power
    upow = 0.d0
    a = 0.d0
    b = 0.d0
    offset = 0.d0

    ! If NCPU is present, distribute the calculation using OpenMP
    if (present(NCPU).and.(.not.(NCPU.eq.-1))) then
        call OMP_SET_NUM_THREADS(NCPU)
    else
        call OMP_SET_NUM_THREADS(1)
    endif

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(x,cosx,sinx,wcosx,wsinx,C,S,YC,YS,CCh,CSh,SSh,CC,SS,CS,D)
    do I=1,NW

        x = W(I) * th
        cosx = cos(x)
        sinx = sin(x)
        wcosx = ww*cosx         ! attach weights
        wsinx = ww*sinx         ! attach weights

        C = sum(wcosx)         ! Eq. (8)
        S = sum(wsinx)         ! Eq. (9)

        YC = sum(yh*wcosx)     ! Eq. (17)
        YS = sum(yh*wsinx)     ! Eq. (18)
        CCh = sum(wcosx*cosx)  ! Eq. (13)
        CSh = sum(wcosx*sinx)  ! Eq. (15)
        SSh = 1.d0 - CCh
        CC = CCh - C*C         ! Eq. (13)
        SS = SSh - S*S         ! Eq. (14)
        CS = CSh - C*S         ! Eq. (15)
        D = CC*SS - CS*CS      ! Eq. (6)

        a(I) = (YC*SS-YS*CS) / D
        b(I) = (YS*CC-YC*CS) / D
        offset(I) = -a(I)*C - b(I)*S

        upow(I) = (SS*YC*YC + CC*YS*YS - 2.d0*CS*YC*YS) / (YY*D) ! Eq. (5)

    end do
    !$OMP END PARALLEL DO

    ! An ad-hoc estimate of the number of independent frequencies 
    ! see discussion following Eq. (24)
    M = (maxval(W/pi2) - minval(W/pi2)) * (maxval(th) - minval(th))

    P = upow
        
END SUBROUTINE GLOMBSCARGLE_EXTRA_OUT