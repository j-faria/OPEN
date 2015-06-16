module params
    use gputils, only: Kernel, GP, WhiteKernel, DiagonalKernel, ExpKernel, & 
                       ExpSquaredKernel, ExpSineSquaredKernel, SumKernels, ProductKernels
    implicit none

! Debug
    logical :: doing_debug

! Data
    real(kind=8), allocatable, dimension(:) :: times, rvs, errors
    real(kind=8), allocatable, dimension(:) :: train_var

! Model Parameters
   	
    !dimensionality
    integer sdim
    !number of planets,
    !number of extra parameters besides the 5*nplanets
    !number of observatories
    integer nplanets, nextra, nobserv
    !considering a jitter parameter
    !using the gaussian process thing
    logical using_jitter, using_gp
    !train the gaussian process beforehand
    logical training
    real(kind=8), dimension(5) :: trained_parameters, trained_std
    !linear/quadratic trend in the model
    logical trend
    integer trend_degree
    !extra linear dependence in the model
    !how many variables
    logical lin_dep
    integer n_lin_dep

    !priors on the parameters are set in main.f90
    real(kind=8), allocatable, dimension(:,:) :: spriorran  !(sdim,2)

    ! nuisance parameters and auxiliary variables
    integer, allocatable, dimension(:) :: observ
    real(kind=8), allocatable, dimension(:) :: ss, alpha, tau
    real(kind=8), allocatable, dimension(:) :: vel, r, sigma
    real(kind=8), allocatable, dimension(:) :: times_oversampled, vel_oversampled
    real(kind=8), allocatable, dimension(:, :) :: last_vel_oversampled
    real(kind=8), allocatable, dimension(:,:) :: covmat, inv_covmat
    real(kind=8) :: mean_time

    integer nest_context

! Gaussian process variables
    type(WhiteKernel), target :: k2
    type(DiagonalKernel), target :: k3
    type(ExpKernel), target :: k4
    type(ExpSquaredKernel), target :: k5
    type(ExpSineSquaredKernel), target :: k6
    type(SumKernels), target :: k7
    type(ProductKernels), target :: k8
    type(GP) gp1
    class(Kernel), pointer :: kernel_to_pass, subk1, subk2


! Parameters for MultiNest
! some of these are set in a namelist for easier access

    !whether to use Nested Importance Sampling
	logical nest_IS
	
    !whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.true.)
	
    !sample with constant efficiency
    ! If ceff is set to True, then the enlargement factor of the 
    ! bounding ellipsoids are tuned so that the sampling efficiency 
    ! is as close to the target efficiency (set by efr) as possible. 
    ! This does mean however, that the evidence value may not be accurate.
	logical nest_ceff
 	parameter(nest_ceff=.false.)
	
    !max no. of live points
    integer nest_nlive
! 	parameter(nest_nlive=500) !300
      
    !total number of parameters, 
    !should be sdim in most cases but if you need to store some 
    !additional parameters with the actual parameters then
    !you need to pass them through the likelihood routine
	integer nest_nPar  ! =sdim
      
    !seed for MultiNest, negative value means take it from sys clock
	integer nest_rseed 
! 	parameter(nest_rseed=1548)
    parameter(nest_rseed=-1)
      
    !evidence tolerance factor
    double precision nest_tol 
    parameter(nest_tol=0.5)
      
    !enlargement factor reduction parameter
    ! defines the sampling efficiency. 0.8 and 0.3 are recommended for 
    ! parameter estimation & evidence evaluation respectively.
    double precision nest_efr
    parameter(nest_efr=0.3d0)
      
    !root for saving posterior files
    character*100 nest_root
	
	!after how many iterations feedback is required and the output
    !files should be updated 
	!posterior files are updated and dumper routine is called
    !after every updInt*10 iterations
	integer nest_updInt ! =10
	
	!null evidence 
    !(set it to very high negative number if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
    !max modes expected, for memory allocation
    integer nest_maxModes 
    parameter(nest_maxModes=10)
      
    !number of parameters to cluster (for mode detection)
    integer nest_nClsPar  ! usually set to 3

    !whether to resume from a previous run
    logical nest_resume
    !parameter(nest_resume=.false.)

    !whether to write output files
    logical nest_outfile
    parameter(nest_outfile=.true.)
      
    !initialize MPI routines?, relevant only if compiling with MPI
	!set to False if you want your main program to handle MPI initialization
    logical nest_initMPI
    parameter(nest_initMPI=.false.)
      
    !points with loglike < nest_logZero will be ignored by MultiNest
    double precision nest_logZero
    parameter(nest_logZero=-huge(1d0))
      
    !max no. of iterations, a non-positive value means infinity. 
    !MultiNest will terminate if either it has done max number of 
    !iterations or convergence criterion (defined through nest_tol) 
    !has been satisfied
    integer nest_maxIter
    !parameter(nest_maxIter=100)
	
	!parameters to wrap around (0 is False; non-zero is True) meaning
    !if nest_pWrap(1) = 0, then parameter 1 is NOT wrap around
	integer, allocatable, dimension(:) :: nest_pWrap  !(sdim)
	
    !feedback on the sampling progress?
    logical nest_fb 

    !feedback on current MAP parameters?
    logical nest_MAPfb
    
    !live plot with pgplot?
    logical nest_liveplot


end module params
