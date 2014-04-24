SUBROUTINE BLOMBSCARGLE_STUB(T, X, W, P, H2BAR, CC, SS, NT, NW)
!f2py threadsafe
!f2py intent(in) T(NT)
!f2py intent(in) X(NT)
!f2py intent(in) W(NW)
!f2py intent(out) P(NW)
!f2py intent(out) H2BAR(NW)
!f2py intent(out) CC(NW)
!f2py intent(out) SS(NW)
!f2py intent(hide) NT
!f2py intent(hide) NW
    IMPLICIT NONE
!* 
!*  Input arguments
!* 
    INTEGER NT, NW
    REAL (KIND=8) T(NT), X(NT), W(NW)
    REAL (KIND=8) P(NW), H2BAR(NW), CC(NW), SS(NW)

!* 
!*  Purpose
!*  =======
!* 
!*  Arguments
!*  =========
!* 
!*  T   (input) REAL (KIND=8) array, dimension (NT)
!*      Sample times
!* 
!*  X   (input) REAL (KIND=8) array, dimension (NT)
!*      Measurement values
!* 
!*  W   (input) REAL (KIND=8) array, dimension (NT)
!*      Angular frequencies for output periodogram
!* 
!*  P   (output) REAL (KIND=8) array, dimension (NW)
!*      Lomb-Scargle periodogram
!* 
!*  NT (input) integer
!*      Dimension of input arrays
!* 
!*  NW (output) integer
!*      Dimension of output array
!* 
!*  Local variables
!* 
    INTEGER I, J
    REAL (KIND=8) pi, pi2, pi4
    REAL (KIND=8) Xmean, Nd2bar
    REAL (KIND=8) theta, twopif, fourpif, cos2pift, sin2pift, cos4pift, sin4pift
    REAL (KIND=8) RR, II, difference

    pi2 = 6.2831853071795862
    pi4 = 12.566370614359172
            
    Xmean = sum(X) / max(1,size(X))

    X = X - Xmean
    Nd2bar = sum(X*X)

    ! step through each frequency 
    do j=1,NW

        twopif  = pi2 * W(j)
        fourpif = pi4 * W(j)

        ! calculate theta
        cos4pift = 0.0
        sin4pift = 0.0
        do i=1,NT
            cos4pift = cos4pift + cos(fourpif*t(i))
            sin4pift = sin4pift + sin(fourpif*t(i))
        end do
        
        theta = 0.5*atan2(sin4pift, cos4pift)

        ! lomb-scargle and p(f|D,I) at this frequency
        RR = 0.0
        II = 0.0
        CC(j) = 0.0
        SS(j) = 0.0
        
        do i=1,NT
        
            cos2pift = cos(twopif*t(i) - theta)
            sin2pift = sin(twopif*t(i) - theta)
            RR = RR + X(i) * cos2pift
            II = II + X(i) * sin2pift
            CC(j) = CC(j) + cos2pift*cos2pift
            SS(j) = SS(j) + sin2pift*sin2pift
        end do 

        ! this is the Lomb-Scargle periodogram:
        h2bar(j) = RR*RR/CC(j) + II*II/SS(j)

        difference = (Nd2bar - h2bar(j))
        ! this is the Bayesian generalization 
        ! see Bretthorst (2000, 2001) or Gregory (2005)
        p(j) = (1.0/sqrt(CC(j)*SS(j))) * difference**((2-NT)/2)
    
    end do

        
        
END SUBROUTINE BLOMBSCARGLE_STUB

