module gputils

    use array_utils
    use random, only: setgmnd, genmnd
    implicit none

    real(kind=8), parameter, private :: pi = 3.1415926535897932384626433832795029d0
    real(kind=8), parameter, private :: const1 = 0.9189385332046727417803297364056176398d0 ! 0.5*ln(2*pi)

    ! extra integers to be saved between calls to mean_fun_keplerian
    integer, save :: gp_n_planets
    integer, save :: gp_n_planet_pars ! this should be 5*nplanets + 1*nobservatories

    ! arrays holding variables on which the RVs depend linearly - used in mean_fun_keplerian_plus_lin
    real(kind=8), allocatable, dimension(:), save :: linvar1, linvar2, linvar3


    ! base type for Gaussian Process kernels
    type Kernel
        real(kind=8), dimension(:), allocatable :: pars
    contains
        !procedure(evaluate_kernel), pointer :: get_matrix_pointer => NULL()
        procedure :: evaluate_kernel
        procedure :: sample_prior
        procedure :: set_kernel_pars, get_kernel_pars
        !procedure :: add_kernel_to_kernel
        !generic :: operator(+) => add_kernel_to_kernel
    end type Kernel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Implemented kernels follow !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! This kernel returns a constant over the whole covariance matrix
    type, extends(Kernel) :: ConstantKernel
    contains
        procedure :: evaluate_kernel => evaluate_kernel_ConstantKernel
    end type ConstantKernel

    ! This kernel returns a constant along the diagonal: "white noise" (homoskedastic)
    type, extends(Kernel) :: DiagonalKernel
    contains
        procedure :: evaluate_kernel => evaluate_kernel_DiagonalKernel
    end type DiagonalKernel

    ! This is a "white noise" (heteroskedastic) kernel
    type, extends(Kernel) :: WhiteKernel
    contains
        procedure :: evaluate_kernel => evaluate_kernel_WhiteKernel
    end type WhiteKernel    

    ! Exponential kernel
    type, extends(Kernel) :: ExpKernel
    contains
        procedure :: evaluate_kernel => evaluate_kernel_ExpKernel
    end type ExpKernel

    ! Exponential squared kernel
    type, extends(Kernel) :: ExpSquaredKernel
    contains
        procedure :: evaluate_kernel => evaluate_kernel_ExpSquaredKernel
    end type ExpSquaredKernel

    ! Exponential sine squared kernel
    type, extends(Kernel) :: ExpSineSquaredKernel
    contains
        procedure :: evaluate_kernel => evaluate_kernel_ExpSineSquaredKernel
    end type ExpSineSquaredKernel

    ! Sum of kernels
    ! this is a hack to substitute operator overloading (which I can't do...)
    ! it should be pretty slow 
    type, extends(Kernel) :: SumKernels
        ! nullify to prevent from needing initialization
        class(Kernel), pointer :: kernel1 => NULL()
        class(Kernel), pointer :: kernel2 => NULL()
        class(Kernel), pointer :: kernel3 => NULL()
        !class(Kernel), pointer :: kernel4 => NULL()
        !class(Kernel), pointer :: kernel5 => NULL()
    contains
        procedure :: evaluate_kernel => evaluate_kernel_sum
        procedure :: set_kernel_pars => set_sub_kernel_pars_sum
        procedure :: get_kernel_pars => get_sub_kernel_pars_sum
    end type SumKernels

    ! Product of kernels
    type, extends(Kernel) :: ProductKernels
        ! nullify to prevent from needing initialization
        class(Kernel), pointer :: kernel1 => NULL()
        class(Kernel), pointer :: kernel2 => NULL()
        class(Kernel), pointer :: kernel3 => NULL()
    contains
        procedure :: evaluate_kernel => evaluate_kernel_product
        procedure :: set_kernel_pars => set_sub_kernel_pars_prod
        procedure :: get_kernel_pars => get_sub_kernel_pars_prod
    end type ProductKernels



    ! user-defined constructor for the Kernel type
    interface Kernel
        module procedure new_Kernel
    end interface 

    ! generic interface for kernel addition subroutines (up to 3 kernels for now)
    interface evaluate_kernel_sum_kernels
        procedure evaluate_kernel_sum_2_kernels
        procedure evaluate_kernel_sum_3_kernels
    end interface

    ! generic interface for kernel multiplication subroutines (up to 3 kernels for now)
    interface evaluate_kernel_product_kernels
        procedure evaluate_kernel_product_2_kernels
        procedure evaluate_kernel_product_3_kernels
    end interface


    ! base type for a Gaussian Process    
    type GP
        ! the mean function associated with this GP
        procedure(mean_fun_template), pointer, nopass :: mean_fun => NULL()
        ! the covariance function associated with this GP
        real(kind=8), dimension(:, :), allocatable :: cov
        ! the kernel associated with this GP
        class(Kernel), pointer :: gp_kernel
        ! possible sub-kernels (we only allow one level down for now)
        class(Kernel), pointer :: sub_kernel1
        class(Kernel), pointer :: sub_kernel2
    contains 
        procedure :: is_posdef
        procedure :: get_lnlikelihood
        procedure :: predict
        procedure :: sample_conditional
    end type GP

    ! user-defined constructor for the GP type
    interface GP
        module procedure new_GP
    end interface 


contains

    function new_Kernel(parameters)
        ! this is the Kernel derived-type constructor
        real(kind=8), dimension(:), intent(in) :: parameters
        type(Kernel) new_Kernel
        integer :: size_pars

        size_pars = size(parameters)

        allocate(new_Kernel%pars(size_pars))
        new_Kernel%pars = parameters

    end function new_Kernel


    function new_GP(cov_matrix, kernel_ptr, subkernel_ptr1, subkernel_ptr2)
        ! this is the GP derived-type constructor
        real(kind=8), dimension(:, :) :: cov_matrix
        class(Kernel), pointer :: kernel_ptr
        class(Kernel), pointer, optional :: subkernel_ptr1, subkernel_ptr2
        type(GP) new_GP
        integer :: shape_cov(2), N

        shape_cov = shape(cov_matrix)
        if (shape_cov(1) /= shape_cov(2)) STOP 'covariance matrix must be square'
        N = shape_cov(1)

        allocate(new_GP%cov(N, N))
        new_GP%cov = cov_matrix
        new_GP%gp_kernel => kernel_ptr

        if(present(subkernel_ptr1)) new_GP%sub_kernel1 => subkernel_ptr1
        if(present(subkernel_ptr2)) new_GP%sub_kernel2 => subkernel_ptr2


    end function new_GP

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Kernel functions and subroutines follow !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function evaluate_kernel(self, x1, x2) result(matrix)
        ! dummy method: subclasses should implement this themselves
        class(Kernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1),size(x2)) :: matrix
        STOP 'Not implemented by parent class'
    end function evaluate_kernel


    function evaluate_kernel_ConstantKernel(self, x1, x2) result(matrix)
        ! Value of the kernel evaluated at a given pair of coordinates x1, x2
        class(ConstantKernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix

        if (size(self%pars) /= 1) STOP 'Wrong parameter dimension - ConstantKernel(1 hyper)'
        call fill_matrix_with_value(matrix, self%pars(1)**2)

    end function evaluate_kernel_ConstantKernel


    function evaluate_kernel_DiagonalKernel(self, x1, x2) result(matrix)
        ! Value of the kernel evaluated at a given pair of coordinates x1, x2
        class(DiagonalKernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix

        if (size(self%pars) /= 1) STOP 'Wrong parameter dimension - DiagonalKernel(1 hyper)'
        ! pars(1) --> amplitude
        matrix = 0.d0
        call add_value_to_diagonal(matrix, self%pars(1))

    end function evaluate_kernel_DiagonalKernel


    function evaluate_kernel_WhiteKernel(self, x1, x2) result(matrix)
        ! Value of the kernel evaluated at a given pair of coordinates x1, x2
        class(WhiteKernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix

        if (size(self%pars) /= size(x1)) STOP 'Wrong parameter dimension - WhiteKernel(diag array)'
        ! pars(:) --> uncertainties on the data points, added in quadrature to the diagonal of the covariance matrix
        matrix = 0.d0
        call add_array_to_diagonal(matrix, self%pars**2)

    end function evaluate_kernel_WhiteKernel


    function evaluate_kernel_ExpKernel(self, x1, x2) result(matrix)
        ! Value of the kernel evaluated at a given pair of coordinates x1, x2
        class(ExpKernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix

        !if (size(self%pars) /= 1) STOP 'Wrong parameter dimension'
        matrix = 0.d0
        matrix = get_radial_coordinates(x1, x2)
        matrix = exp(-abs(matrix))

    end function evaluate_kernel_ExpKernel


    function evaluate_kernel_ExpSquaredKernel(self, x1, x2) result(matrix)
        ! Value of the kernel evaluated at a given pair of coordinates x1, x2
        class(ExpSquaredKernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix

        if (size(self%pars) /= 1) STOP 'Wrong parameter dimension - ExpSquaredKernel(1 hyper)'
        ! pars(1) --> lengthscale
        matrix = 0.d0
        matrix = get_radial_coordinates(x1, x2)
        matrix = exp(-0.5d0 * (matrix / self%pars(1))**2)

    end function evaluate_kernel_ExpSquaredKernel


    function evaluate_kernel_ExpSineSquaredKernel(self, x1, x2) result(matrix)
        ! Value of the kernel evaluated at a given pair of coordinates x1, x2
        class(ExpSineSquaredKernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) :: gamma, omega

        if (size(self%pars) /= 2) STOP 'Wrong parameter dimension - ExpSineSquaredKernel(2 hyper)'
        ! pars(1) --> period
        ! pars(2) --> correlation scale
        omega = pi / self%pars(1)
        gamma = 2.d0 / self%pars(2)**2

        matrix = 0.d0
        !call add_value_to_diagonal(matrix, 1.d0)
        matrix = get_radial_coordinates(x1, x2)
        matrix = exp( -gamma * sin( omega * abs(matrix) )**2 )

    end function evaluate_kernel_ExpSineSquaredKernel


    function evaluate_kernel_sum_2_kernels(k1, k2, x1, x2, a1, a2) result(matrix)
        ! Return the sum of covariance matrices for k1 and k2 at given independent coordinates
        class(Kernel), intent(in) :: k1, k2
        real(kind=8), intent(in), optional :: a1, a2
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) aa1, aa2

        aa1 = 1.d0
        aa2 = 1.d0
        if (present(a1)) aa1 = a1
        if (present(a2)) aa2 = a2
        matrix = aa1 * k1%evaluate_kernel(x1, x2) + aa2 * k2%evaluate_kernel(x1, x2)
        !print *, 'kernel evaluated'

    end function evaluate_kernel_sum_2_kernels


    function evaluate_kernel_sum_3_kernels(k1, k2, k3, x1, x2, a1, a2, a3) result(matrix)
        ! Return the sum of covariance matrices for k1, k2 and k3 at given independent coordinates
        class(Kernel), intent(in) :: k1, k2, k3
        real(kind=8), intent(in), optional :: a1, a2, a3
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) aa1, aa2, aa3

        aa1 = 1.d0
        aa2 = 1.d0
        aa3 = 1.d0
        if (present(a1)) aa1 = a1
        if (present(a2)) aa2 = a2
        if (present(a3)) aa3 = a3

        matrix = aa1 * k1%evaluate_kernel(x1, x2) + aa2 * k2%evaluate_kernel(x1, x2) + aa3 * k3%evaluate_kernel(x1, x2)

    end function evaluate_kernel_sum_3_kernels


    function evaluate_kernel_product_2_kernels(k1, k2, x1, x2, a1, a2) result(matrix)
        ! Return the product of covariance matrices for k1 and k2 at given independent coordinates
        class(Kernel), intent(in) :: k1, k2
        real(kind=8), intent(in), optional :: a1, a2
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) aa1, aa2

        aa1 = 1.d0
        aa2 = 1.d0
        if (present(a1)) aa1 = a1
        if (present(a2)) aa2 = a2

        matrix = aa1 * k1%evaluate_kernel(x1, x2) * aa2 * k2%evaluate_kernel(x1, x2)

    end function evaluate_kernel_product_2_kernels


    function evaluate_kernel_product_3_kernels(k1, k2, k3, x1, x2, a1, a2, a3) result(matrix)
        ! Return the product of covariance matrices for k1, k2 and k3 at given independent coordinates
        class(Kernel), intent(in) :: k1, k2, k3
        real(kind=8), intent(in), optional :: a1, a2, a3
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) aa1, aa2, aa3

        aa1 = 1.d0
        aa2 = 1.d0
        aa3 = 1.d0
        if (present(a1)) aa1 = a1
        if (present(a2)) aa2 = a2
        if (present(a3)) aa3 = a3

        matrix = aa1 * k1%evaluate_kernel(x1, x2) * aa2 * k2%evaluate_kernel(x1, x2) * aa3 * k3%evaluate_kernel(x1, x2)

    end function evaluate_kernel_product_3_kernels


    function evaluate_kernel_sum(self, x1, x2) result(matrix)
        ! Return the sum of covariance matrices for self's associated kernels at given independent coordinates
        class(SumKernels), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) :: aa1, aa2, aa3

        if (.not. associated(self%kernel1)) STOP 'kernel1 is not associated'
        
        if (associated(self%kernel1) .and. associated(self%kernel2) .and. associated(self%kernel3)) then
            ! three kernels are associated
            if (size(self%pars) /= 3) STOP 'Need 3 multiplying constants for 3 associated kernels'
            aa1 = self%pars(1)
            aa2 = self%pars(2)
            aa3 = self%pars(3)
            matrix = evaluate_kernel_sum_3_kernels(self%kernel1, self%kernel2, self%kernel3, x1, x2, a1=aa1, a2=aa2, a3=aa3)

        else if (associated(self%kernel1) .and. associated(self%kernel2)) then
            ! only two kernels are associated
            if (size(self%pars) /= 2) STOP 'Need 2 multiplying constants for 2 associated kernels'
            aa1 = self%pars(1)
            aa2 = self%pars(2)
            matrix = evaluate_kernel_sum_2_kernels(self%kernel1, self%kernel2, x1, x2, a1=aa1, a2=aa2)
        end if

    end function evaluate_kernel_sum


    function evaluate_kernel_product(self, x1, x2) result(matrix)
        ! Return the product of covariance matrices for self's associated kernels at given independent coordinates
        class(ProductKernels), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: matrix
        real(kind=8) :: aa1, aa2, aa3

        if (.not. associated(self%kernel1)) STOP 'kernel1 is not associated'
        
        if (associated(self%kernel1) .and. associated(self%kernel2) .and. associated(self%kernel3)) then
            ! three kernels are associated
            if (size(self%pars) /= 3) STOP 'Need 3 multiplying constants for 3 associated kernels'
            aa1 = self%pars(1)
            aa2 = self%pars(2)
            aa3 = self%pars(3)
            matrix = evaluate_kernel_product_3_kernels(self%kernel1, self%kernel2, self%kernel3, x1, x2, a1=aa1, a2=aa2, a3=aa3)

        else if (associated(self%kernel1) .and. associated(self%kernel2)) then
            ! only two kernels are associated
            if (size(self%pars) /= 2) STOP 'Need 2 multiplying constants for 2 associated kernels'
            aa1 = self%pars(1)
            aa2 = self%pars(2)
            matrix = evaluate_kernel_product_2_kernels(self%kernel1, self%kernel2, x1, x2, a1=aa1, a2=aa2)
        end if

    end function evaluate_kernel_product


    subroutine set_kernel_pars(self, k, p)
        ! dummy method: subclasses should implement this themselves
        class(Kernel), intent(inout) :: self
        integer, intent(in) :: k
        real(kind=8), dimension(:), intent(in) :: p
    end subroutine set_kernel_pars

    subroutine set_sub_kernel_pars_sum(self, k, p)
        class(SumKernels), intent(inout) :: self
        integer, intent(in) :: k
        real(kind=8), dimension(:), intent(in) :: p

        if (.not. associated(self%kernel1)) STOP 'set_sub_kernel_pars_sum : kernel1 is not associated'
        if (k == 1) then
            !print *, self%kernel1%pars
            self%kernel1%pars = p
            !print *, self%kernel1%pars
            return
        end if

        if (.not. associated(self%kernel2)) STOP 'set_sub_kernel_pars_sum : kernel2 is not associated'
        if (k == 2) then
            !print *, self%kernel2%pars
            self%kernel2%pars = p
            !print *, self%kernel2%pars
            return
        end if
    end subroutine set_sub_kernel_pars_sum

    subroutine set_sub_kernel_pars_prod(self, k, p)
        class(ProductKernels), intent(inout) :: self
        integer, intent(in) :: k
        real(kind=8), dimension(:), intent(in) :: p

        if (.not. associated(self%kernel1)) STOP 'set_sub_kernel_pars_prod : kernel1 is not associated'
        if (k == 1) then
            !print *, self%kernel1%pars
            self%kernel1%pars = p
            !print *, self%kernel1%pars
            return
        end if

        if (.not. associated(self%kernel2)) STOP 'set_sub_kernel_pars_prod : kernel2 is not associated'
        if (k == 2) then
            !print *, self%kernel2%pars
            self%kernel2%pars = p
            !print *, self%kernel2%pars
            return
        end if
    end subroutine set_sub_kernel_pars_prod


    subroutine get_kernel_pars(self, k)
        ! dummy method: subclasses should implement this themselves
        class(Kernel), intent(inout) :: self
        integer, intent(in) :: k
    end subroutine get_kernel_pars


    subroutine get_sub_kernel_pars_sum(self, k)
        ! dummy method: subclasses should implement this themselves
        class(SumKernels), intent(inout) :: self
        integer, intent(in) :: k
        if (k==1) print *, self%kernel1%pars
        if (k==2) print *, self%kernel2%pars
    end subroutine get_sub_kernel_pars_sum

    subroutine get_sub_kernel_pars_prod(self, k)
        ! dummy method: subclasses should implement this themselves
        class(ProductKernels), intent(inout) :: self
        integer, intent(in) :: k
        if (k==1) print *, self%kernel1%pars
        if (k==2) print *, self%kernel2%pars
    end subroutine get_sub_kernel_pars_prod


    function sample_prior(self, t) result(sample)
        ! draw samples from the prior at coordinates t
        class(Kernel), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: t
        real(kind=8), dimension(size(t)) :: sample

        real(kind=8), dimension(size(t), size(t)) :: cov
        real(kind=8), dimension(size(t)) :: mean

        cov = self%evaluate_kernel(t, t)
        sample = multivariate_normal_sample(mean, cov)

    end function sample_prior


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! GP functions and subroutines follow !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function mean_fun_template(x, args) result(mean)
        ! dummy method: subclasses should implement this themselves 
        ! by assigning a suitable function to GP%mean_fun
        real(kind=8), dimension(:), intent(in) :: x, args
        real(kind=8), dimension(size(x)) :: mean
    end function mean_fun_template

    function mean_fun_constant(x, args) result(mean)
        ! Constant value mean function - returns args(1)
        real(kind=8), dimension(:), intent(in) :: x, args
        real(kind=8), dimension(size(x)) :: mean

        mean = args(1)
    end function mean_fun_constant

    function mean_fun_keplerian(x, args) result(mean)
        ! Keplerian mean function - args should contain [ (P,K,ecc,omega,t0)*nplanets, vsys ]
        real(kind=8), dimension(:), intent(in) :: x, args
        real(kind=8), dimension(size(x)) :: mean
        integer :: n
        n = size(x)
        print *, 'here'
        ! ATTENTION: the systematic velocity argument below only allows for one observatory!!!
        call get_rvN(x, &
                     args(1:gp_n_planet_pars-1:5), & ! periods for all planets
                     args(2:gp_n_planet_pars-1:5), & ! K for all planets
                     args(3:gp_n_planet_pars-1:5), & ! ecc for all planets
                     args(4:gp_n_planet_pars-1:5), & ! omega for all planets
                     args(5:gp_n_planet_pars-1:5), & ! t0 for all planets
                     args(gp_n_planet_pars), & ! systematic velocity
                     mean, n, gp_n_planets)

    end function mean_fun_keplerian

    function mean_fun_keplerian_plus_lin(x, args) result(mean)
        ! Keplerian mean function plus linear dependence on other variables
        !  args should contain [ (P,K,ecc,omega,t0)*nplanets, vsys ]
        real(kind=8), dimension(:), intent(in) :: x, args
        real(kind=8), dimension(size(x)) :: mean
        integer :: n
        n = size(x)

        ! ATTENTION: the systematic velocity argument below only allows for one observatory!!!
        call get_rvN(x, &
                     args(1:gp_n_planet_pars-1:5), & ! periods for all planets
                     args(2:gp_n_planet_pars-1:5), & ! K for all planets
                     args(3:gp_n_planet_pars-1:5), & ! ecc for all planets
                     args(4:gp_n_planet_pars-1:5), & ! omega for all planets
                     args(5:gp_n_planet_pars-1:5), & ! t0 for all planets
                     args(gp_n_planet_pars), & ! systematic velocity
                     mean, n, gp_n_planets)

        mean = mean + 1.0d0 * linvar1

    end function mean_fun_keplerian_plus_lin


    function is_posdef(self) result(is)
        class(GP), intent(in) :: self
        logical :: is
        real(kind=8), dimension(:, :), allocatable :: cov_cho_factor
        real(kind=8), dimension(:), allocatable :: yy, xsol
        integer :: N, rc

        N = minval(shape(self%cov))
        allocate(yy(N), xsol(N))
        allocate(cov_cho_factor(N,N))
        cov_cho_factor = self%cov

        call choly(1, N, cov_cho_factor, yy, xsol, rc)
        is = .true.
        if (rc == 2) is = .false.

    end function is_posdef

    function get_lnlikelihood(self, x, y, args, yerr) result(lnlike)
        ! Compute the log likelihood of a set of observations x,y under this GP
        ! args is passed directly to GP%mean_fun
        class(GP), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x, y, args
        real(kind=8), dimension(:), intent(in), optional :: yerr
        real(kind=8) :: lnlike
        real(kind=8), dimension(size(y), size(y)) :: cov_cho_factor
        real(kind=8), dimension(size(y)) :: yy, xsol
        real(kind=8), dimension(size(y), 1) :: yy_m  ! 2-dimensional array to go into DPOTRF, DPOTRS
        integer :: N, rc
        !real(kind=8) :: time1, time2

        if (size(x) /= size(y)) STOP 'Dimension mismatch in get_lnlikelihood'
        if (.not. associated(self%mean_fun)) STOP 'GP%mean_fun is not associated'

        N = size(y)
!         call cpu_time(time1)
        cov_cho_factor = self%gp_kernel%evaluate_kernel(x, x)
!         call cpu_time(time2)
!         print *, 'Took ', time2-time1

        ! optionally add data uncertainties in quadrature to covariance matrix
        if (present(yerr)) then
            if (size(yerr) /= N) STOP 'Dimension mismatch in get_lnlikelihood'
            call add_array_to_diagonal(cov_cho_factor, yerr**2)
        end if
        !cov_cho_factor_copy = cov_cho_factor

        ! solve the linear system using the cholesky method
        if (size(args) == 1 .and. args(1) == 0.d0) then
            yy = y ! skip calling the mean function if it's zero
        else
            yy = y - self%mean_fun(x, args)  ! rh
        end if
        yy_m(:, 1) = yy

!         call cpu_time(time1)
        call DPOTRF( 'L', N, cov_cho_factor, N, rc )
        call DPOTRS( 'L', N, 1, cov_cho_factor, N, yy_m, N, rc )
!         call cpu_time(time2)
!         print '("Time for cholesky = ",f6.3," seconds.")',time2-time1

        !call choly(0, N, cov_cho_factor, yy, xsol, rc)
            
        if (rc /= 0) STOP 'Error in cholesky. Matrix is not positive definite?'

        xsol = yy_m(:, 1)
        lnlike = - sum(log(get_diagonal(cov_cho_factor))) ! -0.5*log(det K)
        lnlike = lnlike - N*const1 ! -(N/2)*ln(2*pi)
        lnlike = lnlike - 0.5d0 * dot_product(yy, xsol) ! -0.5*r.T*K.I*r with K.I*r from DPOTRS

    end function get_lnlikelihood


    subroutine predict(self, x, y, args, t, mu, cov, yerr)
        ! Predictive distribution of the GP, conditional on data x, y
        ! (and optional data uncertainties yerr) and calculated on 
        ! coordinates t. args is passed directly to self%mean_fun. 
        ! Predictive mean and covariance in mu and cov, on output.
        class(GP), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x, y, args, t
        real(kind=8), dimension(:), intent(in), optional :: yerr
        real(kind=8), dimension(:), intent(inout) :: mu
        real(kind=8), dimension(:,:), intent(inout) :: cov
        real(kind=8), dimension(size(y)) :: xsol
        real(kind=8), dimension(size(t), size(x)) :: Kxs
        real(kind=8), dimension(size(x), size(t)) :: Kxs_trans
        real(kind=8), dimension(size(y), size(y)) :: cov_cho_factor
        integer :: N, ns, rc

        if (size(x) /= size(y)) STOP 'Dimension mismatch in predict'
        if (size(t) /= size(mu)) STOP 'Dimension mismatch in predict'

        !cov_cho_factor = self%cov
        cov_cho_factor = self%gp_kernel%evaluate_kernel(x, x)
        N = size(y)
        ns = size(t)

        ! optionally add data uncertainties in quadrature to covariance matrix
        if (present(yerr)) then
            if (size(yerr) /= N) STOP 'Dimension mismatch in predict'
            call add_array_to_diagonal(cov_cho_factor, yerr**2)
        end if

        ! solve the linear system using the cholesky method
        call choly(0, N, cov_cho_factor, y-self%mean_fun(x, args), xsol, rc)
        if (rc /= 0) STOP 'Error in choly. Matrix is not positive definite?'

        ! calculate predictive mean
        !print *, 'evaluate kernel with t,x'
        Kxs = self%gp_kernel%evaluate_kernel(t, x)
        mu = matmul(Kxs, xsol) + self%mean_fun(t, args)

        ! calculate predictive covariance
        Kxs_trans = transpose(Kxs)
        !print *, 'evaluate kernel with t,t'
        cov = self%gp_kernel%evaluate_kernel(t, t)
        ! LAPACK routine - DPOTRS solves a system of linear equations A*X = B
        !                  using the Cholesky factorization (for matrix B of rhs)
        call DPOTRS( 'L', N, ns, cov_cho_factor, N, Kxs_trans, N, rc )
        !if (rc /= 0) STOP 'Error in choly. Matrix is not positive definite?'
        ! on exit from DPOTRS, Kxs_trans contains the solution matrix
        cov = cov - matmul(Kxs, Kxs_trans)

    end subroutine predict


    function sample_conditional(self, x, y, args, t) result(csamples)
        ! Draw samples (at t) from the predictive distribution, conditional on x,y.
        ! args is passed directly to self%mean_fun
        class(GP), intent(in) :: self
        real(kind=8), dimension(:), intent(in) :: x, y, args, t
        real(kind=8), dimension(size(t)) :: csamples

        real(kind=8), dimension(size(t)) :: mu
        real(kind=8), dimension(size(t), size(t)) :: cov

        mu = 0.d0
        cov = 0.d0

        call self%predict(x, y, args, t, mu, cov)
        csamples = multivariate_normal_sample(mu, cov)

    end function sample_conditional

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! general functions and subroutines follow !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function get_radial_coordinates(x1, x2) result(r)
        ! Get radial coordinates matrix needed by specific kernels
        real(kind=8), dimension(:) :: x1, x2
        real(kind=8), dimension(size(x1), size(x2)) :: r
        integer :: n1, n2

        !if (size(x1) /= size(x2)) STOP 'Dimension mismatch in get_radial_coordinates'
        n1 = size(x1)
        n2 = size(x2)

        r = (spread(x1, 2, n2) - spread(x2, 1, n1))
    end function get_radial_coordinates


    function multivariate_normal_sample(mean, cov) result(s)
        ! Generate samples from a multidimensional Gaussian with covariance matrix cov
        real(kind=8), dimension(:), intent(in) :: mean
        real(kind=8), dimension(:,:), intent(in) :: cov

        real(kind=8), dimension(size(mean)) :: s
        real(kind=8), dimension(size(mean)) :: work  ! scratch array
        real(kind=8), dimension(:), allocatable :: param

        integer :: p
        p = size(mean)
        allocate(param(p*(p+3)/2 + 1))

        ! set up for genmn
        CALL setgmnd(mean, cov, p, param)
        ! actual call to generate sample (in S)
        CALL genmnd(param, s, work)
    end function multivariate_normal_sample


    !
    ! Following code originally from J-P Moreau, Paris (www.jpmoreau.fr)
    !

    subroutine choly(mode, n, a, b, x, rc)
    !======================================================================
    !*                                                                    *
    !*  This procedure solves linear systems :  a * x = b                 *
    !*  for positive definite symmetric n x n matrices a using the        *
    !*  Cholesky method.                                                  *
    !*                                                                    *
    !*  a must be symmetric and positive definite, or the method will     *
    !*  fail. i.e. for all nonzero vectors y we must have y' * a * y > 0  *
    !*                                                                    *
    !*  cholesky uses only the lower triangle of a.                       *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Application:                                                     *
    !*   ============                                                     *
    !*                                                                    *
    !*      Solve linear systems with positive definite system matrices   *
    !*      efficiently.                                                 *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Control parameter:                                               *
    !*   ==================                                               *
    !*      mode     integer                                              *
    !*       = 0     Factor matrix a and solve linear system              *
    !*       = 1     Compute factorization only                           *
    !*       = 2     Solve system only; for this call the factorization   *
    !*               must be stored in matrix a from chodec. Used to      *
    !*               solve for many right hand sides.                     *
    !*                                                                    *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      n        integer (n > 0)                                      *
    !*               Dimension of matrix a, size of vectors b and x       *
    !*      a        REAL matrix (n,n)                                    *
    !*                 mode = 0, 1: Matrix of the linear system           *
    !*                 mode = 2   : Cholesky factor                       *
    !*      b        REAL vector (n)        (for mode = 0, 2)             *
    !*               Right hand side of system of equations               *
    !*                                                                    *
    !*   Output parameters:                                               *
    !*   ==================                                               *
    !*      a        REAL matrix (n,n)      (for mode = 0, 1)             *
    !*               Cholesky decomposition of input matrix mat           *
    !*      x        REAL vector (n)        (for mode = 0, 2)             *
    !*               solution vector                                      *
    !*                                                                    *
    !*   Return value rc:                                                 *
    !*   ===============                                                  *
    !*      = 0      all ok                                               *
    !*      = 1      n < 1 or false input parameter                       *
    !*      = 2      Matrix is not positive definite                      *
    !*      = 3      wrong value for mode                                 *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Procedures in use:                                               *
    !*   =================                                                *
    !*      chodec()  : Factorization                                     *
    !*      chosol()  : Solving system                                    *
    !*                                                                    *
    !======================================================================
    integer mode, rc, n
    real(kind=8) a(0:n-1,0:n-1), b(0:n-1), x(0:n-1)

      rc=3

      if (n < 1) then
        rc=1
        return
      end if

      if (mode.eq.0) then   !Factor matrix and solve system
        call chodec(n, a, rc)
        if (rc.eq.0) then
          call chosol(n, a, b, x, rc)
        end if
      end if

      if (mode.eq.1) then    !Factor only
        call chodec (n, a, rc)
      end if

      if (mode.eq.2) then    !solve only
        call chosol (n, a, b, x, rc)
      end if

      return
    end subroutine choly


    subroutine chodec(n, a, rc)
    !======================================================================
    !*                                                                    *
    !*  chodec decomposes the symmetric positive definite matrix a.       *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      n        integer (n > 0)                                      *
    !*               Dimension of matrix a                                *
    !*      a        REAL matrix (n,n)                                    *
    !*               Matrix of left hand coefficients                     *
    !*                                                                    *
    !*   Output parameter:                                                *
    !*   ================                                                 *
    !*      a        REAL matix (n,n)                                     *
    !*               Cholesky decomposition in the lower triangle         *
    !*                                                                    *
    !*   Return value rc:                                                 *
    !*   ===============                                                  *
    !*      = 0      all ok                                               *
    !*      = 1      n < 1                                                *
    !*      = 2      Matrix not  positive definite                        *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Functions in use:   dsqrt (square root in double precision)      *
    !*   ================                                                 *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Constants in use:   EPSQUAD (small number)                       *
    !*   ================                                                 *
    !*                                                                    *
    !======================================================================
    integer rc, n
    real(kind=8)  a(0:n-1,0:n-1)
    real(kind=8) EPSQUAD, sum
    integer i, j, k

      EPSQUAD = 1.d-12

      if (n < 1) then
        rc=1                           ! n < 1  error
        return
      end if

      if (a(0,0) < EPSQUAD) then       ! matrix a not positive definite
        rc=2
        return
      end if
      a(0,0) = dsqrt(a(0,0))
      do j = 1, n-1
        a(j,0) = a(j,0) / a(0,0)
      end do

      do i = 1, n-1
        sum = a(i,i)
        do j = 0, i-1
          sum = sum - a(i,j)*a(i,j)
        end do

        if (sum < EPSQUAD) then         ! matrix a not positive definite
          rc=2
          return
        end if
        a(i,i) = dsqrt(sum)
        do j = i+1, n-1
          sum = a(j,i)
          do k = 0, i-1
            sum = sum - a(i,k) * a(j,k)
          end do
          a(j,i) = sum / a(i,i)
        end do
      end do

      rc = 0  ! all Ok
      return
    end subroutine chodec


    subroutine chosol(n, a, b, x, rc)
    !======================================================================
    !*                                                                    *
    !*  chosol finds the solution x of the linear system B' *  B * x = b  *
    !*  for a lower triangular nonsingular matrix B as supplied in chodec.*
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      n        integer  (n > 0)                                     *
    !*               Dimension of matrix a, size of vectors b and x       *
    !*      a        REAL  matrix (n,n)                                   *
    !*               lower triangular matrix as supplied by  chodec       *
    !*      b        REAL  vector (n)                                     *
    !*               Right hand side                                      *
    !*                                                                    *
    !*   Output parameter:                                                *
    !*   =================                                                *
    !*      x        REAL  vector (n)                                     *
    !*               solution vector                                      *
    !*                                                                    *
    !*   Return value rc:                                                 *
    !*   ===============                                                  *
    !*      = 0      all ok                                               *
    !*      = 1      improper lwer triangular matrix or  n < 1            *
    !*                                                                    *
    !======================================================================
    integer rc, n
    real(kind=8) a(0:n-1,0:n-1), b(0:n-1), x(0:n-1)
    real(kind=8) ZERO, sum
    integer j, k

      ZERO = 0.d0

      if (n < 1) then   ! n < 1 error
        rc=1
        return
      end if

      if (a(0,0).eq.ZERO) then               ! improper factor matrix
        rc=1
        return
      end if

      x(0) = b(0) / a(0,0)                   ! update right hand side
      do k = 1, n-1
        sum = ZERO
        do j = 0, k-1
          sum = sum + a(k,j) * x(j)
        end do
        if (a(k,k).eq.ZERO) then
          rc=1
          return
        end if
        x(k) = (b(k) - sum) / a(k,k)
      end do

      x(n-1) = x(n-1) / a(n-1,n-1)            ! back substitution
      do k=n-2, 0, -1
        sum = ZERO
        do j = k+1, n-1
          sum = sum + a(j,k) * x(j)
        end do
        x(k) = (x(k) - sum) / a(k,k)
      end do

      rc = 0  !all ok
      return

    end subroutine chosol

end module gputils