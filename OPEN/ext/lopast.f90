!      Fortran 90 Toolkit to Accompany "A Class of Fast Methods for
!      Processing Irregularly Sampled or Otherwise Inhomogeneous
!      One-Dimensional Data" by George B. Rybicki and William H. Press
!
! The function lopast(y,x,f) returns a vector that is the low-pass
! filtered version of the data points y, at the abcissas x, with a 3dB
! cutoff frequency of f.  This is pretty much exactly as described in
! the paper.
!
! High pass filtering is obtained by changing the line
!	lopast=y-real(ans)
! to
!	lopast=real(ans)
! and (if you want f to be the 3dB point again) changing the
! constant 5.5380744 to 3.56427.  (If the filters were characterized
! by 6dB cutoffs instead of 3dB cutoffs, the constants would of course
! be the same -- this is only semantics in the definition of cutoff.)

SUBROUTINE lopast(y, x, f, y_lopast, ny)
    IMPLICIT NONE
!f2py intent(in) y(ny)
!f2py intent(in) x(ny)
!f2py intent(in) f
!f2py intent(out) y_lopast(ny)
!f2py intent(hide) ny
    REAL(kind=8), DIMENSION(ny), INTENT(IN) :: y, x
    REAL(kind=8), INTENT(IN) :: f
    REAL(kind=8), DIMENSION(ny), INTENT(OUT)  :: y_lopast
    INTEGER :: ny
    INTEGER :: n, rc
    COMPLEX(kind=8) :: fac
    COMPLEX(kind=8), DIMENSION(size(y)) :: b,rhs,ans
    COMPLEX(kind=8), DIMENSION(size(y)-1) :: r,e,w,re

    n=size(y)
    if (size(x) /= n) STOP 'bad sizes in lopast'
    !fac=abs(f)*5.5380744*cmplx(1.,1.)  ! low-pass filtering
    fac=abs(f)*3.56427*cmplx(1.,1.)  ! high-pass filtering
    w=fac*abs(x(2:n)-x(1:n-1))
    r=exp(-w)  ! note this is a complex exponential !
    e=1./(1./r - r)
    re=r*e
    b=1.
    b(1:n-1)=b(1:n-1)+re
    b(2:n)=b(2:n)+re
    w = 0.5/w
    rhs(1:n-1)=w*(y(1:n-1)-y(2:n))
    rhs(n)=0.
    rhs(2:n)=rhs(2:n)+w*(y(2:n)-y(1:n-1))
    e=-e
    rc = ctridag(e,b,e,rhs,ans)
    !y_lopast=y-real(ans)  ! low-pass filtering
    y_lopast=real(ans)  ! high-pass filtering
    return

CONTAINS

    FUNCTION ctridag(a,b,c,r,u)
        IMPLICIT NONE
        COMPLEX(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b,c,r,u
        INTEGER :: ctridag
        INTEGER :: n,j
        COMPLEX(kind=8), DIMENSION(size(b)) :: gam
        COMPLEX(kind=8) :: bet
        n=size(b)
        if (size(a) /= n-1 .or. size(c) /= n-1 .or. size(r) /= n &
            .or. size(u) /= n) STOP 'bad sizes in ctridag'
        bet=b(1)
        if(bet==(0.,0.)) STOP 'Error 1 in ctridag'
        u(1)=r(1)/bet
        do j=2,n
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j-1)*gam(j)
            if(bet==(0.,0.)) STOP 'Error 2 in ctridag'
            u(j)=(r(j)-a(j-1)*u(j-1))/bet
        enddo
        do j=n-1,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
        enddo
        ctridag = 0
    END FUNCTION ctridag

END SUBROUTINE lopast