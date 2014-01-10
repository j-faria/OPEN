module params
    implicit none

! Debug
    logical :: doing_debug

! Data
    real(kind=8), allocatable, dimension(:) :: times, rvs, errors

! Model Parameters
   	
    !dimensionality
    integer sdim
   	parameter( sdim = 6 )

    !priors on the parameters are set in main.f90
    double precision spriorran(sdim,2)

    ! nuisance parameters and auxiliary variables
    integer, allocatable, dimension(:) :: observ
    real(kind=8), allocatable, dimension(:) :: ss, alpha, tau
    real(kind=8), allocatable, dimension(:) :: vel, dist
    real(kind=8), allocatable, dimension(:,:) :: covmat, inv_covmat

    !real(kind=8), dimension(3,119) :: context
    integer context      


! Parameters for MultiNest

    !whether to do use Nested Importance Sampling
	logical nest_IS
 	parameter(nest_IS=.true.)
	
    !whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.false.)
	
    !sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.true.)
	
    !max no. of live points
    integer nest_nlive
	parameter(nest_nlive=100) !300
      
    !total number of parameters, 
    !should be sdim in most cases but if you need to store some 
    !additional parameters with the actual parameters then
    !you need to pass them through the likelihood routine
	integer nest_nPar 
	parameter(nest_nPar=sdim)
      
    !seed for MultiNest, negative value means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-1)
      
    !evidence tolerance factor
    double precision nest_tol 
    parameter(nest_tol=0.5)
      
    !enlargement factor reduction parameter
    double precision nest_efr
    parameter(nest_efr=0.8d0)
      
    !root for saving posterior files
    character*100 nest_root
	parameter(nest_root='chains/nest-')
	
	!after how many iterations feedback is required and the output
    !files should be updated 
	!posterior files are updated and dumper routine is called
    !after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=10)
	
	!null evidence 
    !(set it to very high negative number if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
    !max modes expected, for memory allocation
    integer nest_maxModes 
    parameter(nest_maxModes=10)
      
    !number of parameters to cluster (for mode detection)
    integer nest_nClsPar
    parameter(nest_nClsPar=2)

    !whether to resume from a previous run
    logical nest_resume
    parameter(nest_resume=.false.)

    !whether to write output files
    logical nest_outfile
    parameter(nest_outfile=.true.)
      
    !initialize MPI routines?, relevant only if compiling with MPI
	!set to False if you want your main program to handle MPI initialization
    logical nest_initMPI
    parameter(nest_initMPI=.true.)
      
    !points with loglike < nest_logZero will be ignored by MultiNest
    double precision nest_logZero
    parameter(nest_logZero=-huge(1d0))
      
    !max no. of iterations, a non-positive value means infinity. 
    !MultiNest will terminate if either it has done max number of 
    !iterations or convergence criterion (defined through nest_tol) 
    !has been satisfied
    integer nest_maxIter
    parameter(nest_maxIter=10000)
	
	!parameters to wrap around (0 is False; non-zero is True) meaning
    !if nest_pWrap(1) = 0, then parameter 1 is NOT wrap around
	integer nest_pWrap(sdim)
	
    !feedback on the sampling progress?
    logical nest_fb 
    parameter(nest_fb=.true.)


end module params
