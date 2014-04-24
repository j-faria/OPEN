module problem1

  real(kind=8), parameter :: pi = 3.1415926535897932384626433832795029d0
  real(kind=8), parameter :: twopi = 2.d0 * pi 

  integer (kind=4) :: likelihood_calls

! radial velocity data
  integer(kind=4) n ! I'm not sure if N is going to be needed outside this module...
  real(kind=8), allocatable, dimension(:) :: times, rvs, errors
! nuisance parameters and auxiliary variables
  integer, allocatable, dimension(:) :: observ
  real(kind=8), allocatable, dimension(:) :: ss, alpha, tau
  real(kind=8), allocatable, dimension(:) :: vel, dist
  real(kind=8), allocatable, save, dimension(:) :: last_vel
  real(kind=8), allocatable, dimension(:,:) :: covmat, inv_covmat
! parameter limits for prior
  real(kind=8), allocatable, dimension(:,:) :: parameter_limits
! number of planets in the model
  integer(kind=4) :: np

contains

  subroutine problem_size (chain_num, cr_num, gen_num, pair_num, par_num)
    ! sets the information related to the dimensions of the problem.
    !  Parameters:
    !
    !    Output, integer (kind = 4) CHAIN_NUM, the total number of chains. (>= 3)
    !    Output, integer (kind = 4) CR_NUM, the total number of CR values. (>= 1)
    !    Output, integer (kind = 4) GEN_NUM, the total number of generations. (>= 2)
    !    Output, integer (kind = 4) PAIR_NUM, the number of pairs of crossover chains. (>= 0)
    !    Output, integer (kind = 4) PAR_NUM, the total number of parameters. (>= 1)
    !
    implicit none

    integer (kind = 4) chain_num
    integer (kind = 4) cr_num
    integer (kind = 4) gen_num
    integer (kind = 4) pair_num
    integer (kind = 4) par_num

    integer (kind = 4) iou, ierr
    namelist /config_size/ chain_num, cr_num, gen_num, pair_num, par_num

    ! read configuration values from namelist
    call get_unit(iou)
    open(unit=iou, file="./OPEN/demc/namelist1", status='old', action='read', delim='quote', iostat=ierr)
    if (ierr /= 0) then
      stop 'Failed to open namelist file namelist1'
    else
      read(iou, nml=config_size)
      close(iou)
    end if

    if (read_data() /= 0) stop 'Error reading data file'
    likelihood_calls = 0

    return
  end subroutine problem_size

  subroutine problem_value (chain_filename, gr_filename, gr_threshold, &
    jumpstep, limits, par_num, printstep, restart_read_filename, &
    restart_write_filename )
    ! sets further config information, including numeric data.
    !  Parameters:
    !
    !    Output, character (len = *) CHAIN_FILENAME, the "base" filename
    !    to be used for the chain files.  If this is the empty string '',
    !    then the chain files will not be written.  This name should 
    !    include a string of 0's which will be replaced by the chain 
    !    indices.  For example, "chain000.txt" would work as long as the
    !    number of chains was 1000 or less.
    !
    !    Output, character (len = *) GR_FILENAME, the name of the file
    !    in which values of the Gelman-Rubin statistic will be recorded,
    !    or '' if this file is not to be written.
    !
    !    Output, real (kind = 8) GR_THRESHOLD, the convergence tolerance for the Gelman-Rubin statistic.
    !
    !    Output, integer (kind = 4) JUMPSTEP, forces a "long jump" every JUMPSTEP generations.
    !
    !    Output, real (kind = 8) LIMITS(2,PAR_NUM), lower and upper bounds for each parameter.
    !
    !    Input, integer (kind = 4) PAR_NUM, the total number of parameters. ( >= 1)
    !
    !    Output, integer (kind = 4) PRINTSTEP, the interval between generations 
    !    on which the Gelman-Rubin statistic will be computed and written to a file.
    !
    !    Output, character (len = *) RESTART_READ_FILENAME, the name of the file
    !    containing restart information, or '' if this is not a restart run.
    !
    !    Output, character (len = *) RESTART_WRITE_FILENAME, the name of the file
    !    to be written, containing restart information, or '' if a restart file 
    !    is not to be written.
    !
    implicit none

    integer (kind = 4) par_num

    character (len = *) chain_filename
    character (len = *) gr_filename
    real (kind = 8) gr_threshold
    integer (kind = 4) jumpstep
    real (kind = 8) limits(2,par_num)
    logical km_per_second
    integer (kind = 4) printstep
    character (len = *) restart_read_filename
    character (len = *) restart_write_filename

    integer (kind=4) i

    integer (kind = 4) iou, ierr
    namelist /config_value/ chain_filename, gr_filename, gr_threshold, jumpstep, &
                            limits, km_per_second, &
                            printstep, restart_read_filename, restart_write_filename

    ! read configuration values from namelist
    call get_unit(iou)
    open(unit=iou, file="./OPEN/demc/namelist1", action='read', delim='quote', iostat=ierr)
    if (ierr /= 0) then
      stop 'Failed to open namelist file namelist1'
    else
      read(iou, nml=config_value)
      close(iou)
    end if

    ! since in the namelist we only specify limits(:, 1:6) we have to do this 
    ! in order to apply these limits to two or more planets
    ! move systematic velocity to last position
    limits(1,par_num) = limits(1,6)
    limits(2,par_num) = limits(2,6)

    if (km_per_second) then
      limits(1,2) = limits(1,2) / 1000.d0
      limits(2,2) = limits(2,2) / 1000.d0
    endif

    ! repeat the other parameters np times
    np = (par_num - 1) / 5 ! number of planets
    do i=0,np-1
      limits(1, 5*i+1:5*(i+1)) = limits(1,1:5)
      limits(2, 5*i+1:5*(i+1)) = limits(2,1:5)
    end do


    allocate(parameter_limits(2, par_num))
    parameter_limits = limits

    return
  end subroutine problem_value


  function prior_density (par_num, zp)
    ! evaluates the prior density function.
    !  Parameters:
    !
    !    Input, integer (kind = 4) PAR_NUM, the total number of parameters. ( >= 1)
    !
    !    Input, real (kind = 8) ZP(PAR_NUM), the argument of the density function.
    !
    !    Output, real (kind = 8) PRIOR_DENSITY, the value of the prior density function.
    implicit none

    integer (kind = 4) par_num
    real (kind = 8) zp(par_num)
    
    real (kind = 8) :: a
    real (kind = 8) :: b
    integer (kind = 4) i
    real (kind = 8) prior_density
    real (kind = 8) r8_uniform_pdf, r8_beta_pdf

    prior_density = 1.0D+00
    
    do i = 1, par_num
      if (mod(i,5) == 3) then
        ! informative beta prior for the eccentricity, based on Kipping (2013)
        prior_density = prior_density * r8_beta_pdf(0.867d0, 3.03d0, zp(i))
      else
        a = parameter_limits(1, i)
        b = parameter_limits(2, i)
        prior_density = prior_density * r8_uniform_pdf (a, b, zp(i))
      end if
    end do

    return
  end function prior_density

  subroutine prior_sample (par_num, zp)
    ! samples from the prior distribution.
    !    Input, integer (kind = 4) PAR_NUM, the total number of parameters. ( >= 1)
    !    Output, real (kind = 8) ZP(PAR_NUM), the sample from the distribution.
    implicit none

    integer (kind = 4) par_num
    real (kind = 8) zp(par_num)

    real (kind = 8) :: a
    real (kind = 8) :: b
    integer (kind = 4) i
    real (kind = 8) r8_uniform_sample, r8_beta_sample

    do i = 1, par_num
      a = parameter_limits(1, i)
      b = parameter_limits(2, i)
      zp(i) = r8_uniform_sample (a, b) 
    end do

    ! informative beta prior for the eccentricity, based on Kipping (2013)
    ! aa = 0.867 (+ 0.044 - 0.044)
    ! bb = 3.03 (+ 0.17 - 0.16)
    do i = 3, par_num-1, 5
      zp(i) = r8_beta_sample(0.867d0, 3.03d0)
    end do

    return
  end subroutine prior_sample

  integer function read_data()
    ! reads data file and sets problem dimensions that depend on it
    ! returns status
    character(len=10) line
    integer i

    !! load data
    open(unit=15, file='input.rv', status="old")
    read(15, *) line ! header
    read(15, '(a)') line ! number of values
    line = trim(line(index(line, "#")+1:))
    read(line, *) n ! n is total number of points in data-set
    
    ! allocate storage
    allocate(times(N), rvs(N), errors(N))
    ! read data-set
    do i = 1, n
        read(15, *) times(i), rvs(i), errors(i)
    end do
    close(15)
    
    ! allocate storage for nuisance parameters and auxiliary variables
    allocate(ss(N), alpha(N), tau(N))
    allocate(observ(N))
    allocate(vel(N), dist(N))
    allocate(last_vel(N))
    allocate(covmat(N,N), inv_covmat(N,N))

    read_data = 0
  end function

  subroutine deallocate_problem()
    deallocate(parameter_limits)
    deallocate(times, rvs, errors)
    deallocate(ss, alpha, tau)
    deallocate(observ)
    deallocate(vel, dist)
    deallocate(last_vel)
    deallocate(covmat, inv_covmat)
  end subroutine deallocate_problem

  real(kind=8) function sample_likelihood(par_num, zp)
    ! computes the log-likelihood 
    !    Input, integer PAR_NUM, the total number of parameters. ( >= 1)
    !    Input, real(kind=8) ZP(PAR_NUM), a sample.
    !    Output, real(kind=8) SAMPLE_LIKELIHOOD, the log likelihood function for the sample.
    implicit none

    integer (kind=4) par_num
    real (kind=8) zp(par_num)

    real(kind=8) :: lhood, lhood_test(1,1), det
    integer (kind=4) :: i

    ! times, rvs and errors are initialized/read in READ_DATA
    ! P, K, ecc, omega, t0, Vsys
    observ = 1
    ss = 1.d0
    alpha = 1.d0
    tau = 1.d0

    sample_likelihood=-huge(1.d0)*epsilon(1.d0)

    ! get the radial velocity model with these parameters (in vel)
    !write(*,'(6f10.2)') zp
    !print *, zp(1:par_num-1:5)
    call get_rvN(times, &
                 zp(1:par_num-1:5), & ! periods for all planets
                 zp(2:par_num-1:5), & ! K for all planets
                 zp(3:par_num-1:5), & ! ecc for all planets
                 zp(4:par_num-1:5), & ! omega for all planets
                 zp(5:par_num-1:5), & ! t0 for all planets
                 zp(par_num), vel, n, 1)
    !write(*,*) 'out of get_rvN'
    dist = rvs - vel
    !call get_covmat(times, errors, observ, ss, alpha, tau, covmat, det, inv_covmat)

    !lhood_test = -0.5d0 * matmul(matmul(reshape(dist, (/1,n/)), covmat), reshape(dist, (/n,1/)))
    !lhood = lhood_test(1,1) - 0.5d0*log(twopi**n * det)
    lhood = -0.5d0 * sum(dist**2 / errors**2) - 0.5d0*log(twopi**n * product(errors))

    sample_likelihood = logSumExp(sample_likelihood, lhood)
    likelihood_calls = likelihood_calls + 1

  end function

  subroutine get_covmat(times, sigma, observ, ss, alpha, tau, covmat, det, inv_covmat)
    ! Calculates the covariance matrix of the observations and returns its 
    ! determinant and inverse.

    ! vectors with times and uncertainties
    real(kind=8), intent(in), dimension(:) :: times, sigma
    ! vector with flags for points comming from each observatory
    integer, intent(in), dimension(:) :: observ
    ! nuisance parameters for each observatory
    real(kind=8), intent(in), dimension(:) :: ss, alpha, tau
    ! on output, the covariance matrix of the observations and its inverse
    real(kind=8), intent(out), dimension(:,:) :: covmat, inv_covmat
    real(kind=8), intent(out) :: det ! determinant of covmat

    ! local variables
    integer :: i, j, nt

    covmat = 0.d0; inv_covmat = 0.d0
    nt = size(sigma)
    do i=1,nt
      do j=1,nt
        ! Kronecker delta on the times
        if (i==j) covmat(i,j) = (sigma(j)/alpha(j))**2
        ! Kronecker delta on the observatories
        if (observ(i)==observ(j)) then
          covmat(i,j) = covmat(i,j) + ss(j)**2 * exp(-abs(times(i)-times(j))/tau(j))
        endif
      end do
    end do  
    det = determinant(covmat)
    call inverse(covmat, inv_covmat)
  end subroutine get_covmat     


  subroutine inverse(matrix, inverse_matrix)
    !*******************************************************************************
    ! inversion of a matrix (uses LAPACK dgetri and dgetrf routines)
    ! doesn ot do proper error checking
    !*******************************************************************************
    real(kind=8), intent(in)    :: matrix(:, :)
    real(kind=8), intent(out)   :: inverse_matrix(:, :)
    integer, allocatable    :: ipiv(:)
    real(kind=8), allocatable   :: work(:)
    integer :: n, lwork, lda, info
    integer :: i, j

    n = size(matrix, dim=1)
    if (size(matrix, dim=2) /= n) then
        stop 'inverse: matrix must be square'
    end if
    
    lda = n
    lwork = n*n
    allocate(ipiv(n))
    allocate(work(lwork))

    ! store MATRIX in INVERSE_MATRIX to prevent it from being overwritten by LAPACK
    inverse_matrix = matrix
    ! DGETRF computes an LU factorization of a general M-by-N matrix, using 
    ! partial pivoting with row interchanges
    call dgetrf(n, n, inverse_matrix, lda, ipiv, info )
    ! DGETRI computes the inverse of a matrix using the LU factorization 
    ! computed by DGETRF
    call dgetri(n, inverse_matrix, n, ipiv, work, lwork, info)

    deallocate (ipiv)
    deallocate (work)
    
    return
  end subroutine inverse

  real(kind=8) function determinant(matrix) result(det)
    !*******************************************************************************
    ! Calculate the determinant of a real square matrix A(n,n) by Gauss method 
    ! with full pivoting.
    !
    ! Ref: "Algèbre - Algorithmes et programmes en Pascal By Jean-Louis Jardrin, 
    !      Dunod - Bordas Paris, 1988 p. 76-79" [BIBLI 10].
    !*******************************************************************************
    real(kind=8) :: matrix(:,:)
    integer :: n !size of matrix
    real(kind=8) :: eps !desired precision
    eps = epsilon(1.d0)
    n = size(matrix, dim=1)

    det = DMGT(eps, n, matrix)
    
    contains
    
      !The subroutine TSRGT applies to input real square matrix A(n,n) the upper
      !triangularization algorithm of Gauss method with full pivoting and keeps
      !trace of successive transformations done in integer vectors KP and LP.
      !-------------------------------------------------------------------------
      !  Input parameters:
      !  eps        precision (real*8)
      !  n          size of A matrix (integer)
      !  A          pointer to input real square matrix (real*8)
      !  Output parameters:
      !  it         flag=1 if A matrix ok, =0 if A matrix is singular (integer)
      !  C          pointer to table storing main diagonal elements and supra-
      !             diagonal elements of upper triangular matrix and the multi-
      !             plying coefficients used during triangularization process
      !  KP         table storing informations concerning the column exchanges
      !             during process (integer)
      !  LP         table storing informations concerning the line exchanges
      !             during process (integer)
      !-------------------------------------------------------------------------
      !The table C is first initialized to A matrix, then receives at each step 
      !k of the triangularization process, usefull elements of A matrix at step 
      !k for
      !k=1,2,...n.
      !The variables po(real*8), lo and ko(integer) store respectively pivot at 
      !step k, its line number and its column number.
      !-------------------------------------------------------------------------
      subroutine TSRGT(eps, n, A, it, C, Kp, Lp)
        real(kind=8) :: eps
        integer :: n, it
        real(kind=8), dimension(n,n) :: A, C
        integer, dimension(n) :: Kp, Lp
        integer :: i, j, k, ko, lo
        real(kind=8) :: po,t0
        C=A; it=1; k=1
        do while (it==1.and.k<n)
          po=C(k,k); lo=k; ko=k
          do i=k, n
            do j=k, n
              if (abs(C(i,j))>abs(po)) then
                po=C(i,j); lo=i; ko=j
              end if
            end do
          end do
          Lp(k)=lo; Kp(k)=ko
          if (abs(po)<eps) then
            it=0
          else
            if (lo.ne.k) then
              do j=k, n
                t0=C(k,j); C(k,j)=C(lo,j); C(lo,j)=t0
              end do
            end if
            if (ko.ne.k) then
              do i=1, n
                t0=C(i,k); C(i,k)=C(i,ko); C(i,ko)=t0
              end do
            end if 
            do i=k+1, n
              C(i,k)=C(i,k)/po
              do j=k+1, n
                C(i,j)=C(i,j)-C(i,k)*C(k,j)
              end do 
            end do
            k=k+1
          end if
        end do
        if (it==1.and.abs(C(n,n))<eps)  it=0
        return
      end subroutine TSRGT !TSRGT

      !The function DMGT returns the determinant of a real square matrix
      !A(n,n) by Gauss method with full pivoting.
      !----------------------------------------------------------------------------
      !  Input parameters:
      !  eps        precision (real*8)
      !  n          size of A matrix (integer)
      !  A          pointer to input real square matrix
      !  Output parameters:
      !  None
      !-----------------------------------------------------------------------------
      !The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
      !Output variables are it(integer), C(n,n), Kp(n) and Lp(n).
      !If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
      !at location i,j (j>=i) the corresponding element of the upper triangular matrix.
      !Tables Lp and Kp contain informations relative to exchanges of line or column
      !that occured during the process. For instance, the element number k of Lp is
      !an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
      !The number of exchanges of lines and columns is stored in integer L. the
      !determinant of A matrix is stored in d0 (real*8).
      !-----------------------------------------------------------------------------
      real(kind=8) function DMGT(eps, n, A)
        integer :: n
        real(kind=8) :: eps, A(n,n)
        real(kind=8) :: d0

        real(kind=8), pointer :: C(:,:)
        integer,pointer :: Kp(:), Lp(:)
        integer :: it, k, l, i, j

        !allocate local matrix C and vectors Kp, Lp
        allocate(C(n,n))
        allocate(Kp(n))
        allocate(Lp(n))

        call TSRGT(eps,n,A,it,C,Kp,Lp)  !call triangularization subroutine
        if (it==0) then
          d0=0.d0  !matrix singular, det=0
        else       !matrix regular, det<>0
          d0=1.d0
          do k=1, n
          d0=d0*C(k,k)
          end do
          l=0
          do k=1, n-1
            if (Lp(k).ne.k)  l=l+1
            if (Kp(k).ne.k)  l=l+1
          end do
          if (MOD(l,2).ne.0) d0=-d0  !l is odd
        end if
        DMGT=d0   !return determinant
        return
      end function DMGT
  end function determinant

  real(kind=8) function LogSumExp(x,y)
    ! LogSumExp(x,y)=log(exp(x)+exp(y))
    real (kind=8) x,y

    if (x.gt.y) then
       LogSumExp=x+log(1+exp(y-x))
    else
       LogSumExp=y+log(1+exp(x-y))
    end if
    return
  end function LogSumExp
end module problem1