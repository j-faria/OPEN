subroutine main
!f2py threadsafe
!main program for DREAM.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!
!  Author:
!    Original FORTRAN90 version by Guannan Zhang.
!    Modifications by John Burkardt.
!
!  Reference:
!
!    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
!    Dave Higdon,
!    Accelerating Markov Chain Monte Carlo Simulation by Differential 
!    Evolution with Self-Adaptive Randomized Subspace Sampling,
!    International Journal of Nonlinear Sciences and Numerical Simulation,
!    Volume 10, Number 3, March 2009, pages 271-288.
!
!  Local parameters:
!
!    Local, character (len = 255) CHAIN_FILENAME, the "base" filename
!    to be used for the chain files.  If this is the empty string '',
!    then the chain files will not be written.  This name should 
!    include a string of 0's which will be replaced by the chain 
!    indices.  For example, "chain000.txt" would work as long as the
!    number of chains was 1000 or less.
!
!    Local, integer (kind = 4) CHAIN_NUM, the total number of chains.
!    3 <= CHAIN_NUM.
!
!    Local, integer (kind = 4) CR_NUM, the total number of CR values.
!    1 <= CR_NUM.
!
!    Local, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
!    each sample.
!
!    Local, integer (kind = 4) GEN_NUM, the total number of generations.
!    2 <= GEN_NUM.
!
!    Local, real (kind = 8) GR(PAR_NUM,GR_NUM), 
!    the Gelman-Rubin R statistic.
!
!    Local, logical GR_CONV, the Gelman-Rubin convergence flag.
!
!    Local, integer (kind = 4) GR_COUNT, counts the number of generations
!    at which the Gelman-Rubin statistic has been computed.
!
!    Local, character (len = 255) GR_FILENAME, the name of the file
!    in which values of the Gelman-Rubin statistic will be recorded,
!    or '' if no such file is to be created.
!
!    Local, integer (kind = 4) GR_NUM, the number of times the Gelman-Rubin
!    statistic may be computed.
!
!    Local, real (kind = 8) GR_THRESHOLD, the convergence tolerance for
!    the Gelman-Rubin statistic.
!
!    Local, real (kind = 8) JUMPRATE_TABLE(PAR_NUM), the jumprate table.
!
!    Local, integer (kind = 4) JUMPSTEP, forces a "long jump" every
!    JUMPSTEP generations.
!
!    Local, real (kind = 8) LIMITS(2,PAR_NUM), lower and upper bounds
!    for each parameter.
!
!    Local, integer (kind = 4) PAIR_NUM, the number of pairs of 
!    crossover chains.
!    0 <= PAIR_NUM.
!
!    Local, integer (kind = 4) PAR_NUM, the total number of parameters.
!    1 <= PAR_NUM.
!
!    Local, integer (kind = 4) PRINTSTEP, the interval between generations on 
!    which the Gelman-Rubin statistic will be computed and written to a file.
!
!    Local, character (len = 255) RESTART_READ_FILENAME, the name of the file
!    containing restart information, or '' if this is not a restart run.
!
!    Local, character (len = 255) RESTART_WRITE_FILENAME, the name of the file
!    to be written, containing restart information, or '' if a restart file 
!    is not to be written.
!
!    Local, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
!    sample data.
!
  use problem1, only: problem_size, problem_value, deallocate_problem
  implicit none

  character (len = 255) chain_filename
  integer (kind = 4) chain_num
  integer (kind = 4) cr_num
  real (kind = 8), allocatable :: fit(:,:)
  integer (kind = 4) gen_num
  real (kind = 8), allocatable :: gr(:,:)
  logical gr_conv
  integer (kind = 4) gr_count
  character (len = 255) gr_filename
  integer (kind = 4) gr_num
  real (kind = 8) gr_threshold
  real (kind = 8), allocatable :: jumprate_table(:)
  integer (kind = 4) jumpstep
  real (kind = 8), allocatable :: limits (:, :)
  integer (kind = 4) pair_num
  integer (kind = 4) par_num
  integer (kind = 4) printstep
  real (kind = 8) r8_uniform_01_sample
  character (len = 255) restart_read_filename
  character (len = 255) restart_write_filename
  real (kind = 8), allocatable :: z(:,:,:)

  ! OpenMP variables
  integer (kind=4) OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

  integer pgopen

  call timestamp ()

  write (*, '(a)') ' '
  write (*, '(a)') 'DREAM'
  write (*, '(a)') '  MCMC acceleration by Differential Evolution.'
!$OMP PARALLEL
  if (OMP_GET_THREAD_NUM() == 0) then
    write (*, '(a, i1, a)') '  Using ', OMP_GET_NUM_THREADS(), &
                            ' threads for parallelization'
  end if
!$OMP END PARALLEL
  
  !  Get the problem sizes.
  call problem_size (chain_num, cr_num, gen_num, pair_num, par_num)

  !  Decide if the problem sizes are acceptable.
  if (chain_num < 3) then
    write (*, '(a)') ' '
    write (*, '(a)') 'DREAM - Fatal error!'
    write (*, '(a)') '  CHAIN_NUM < 3.'
    stop
  end if

  if (cr_num < 1) then
    write (*, '(a)') ' '
    write (*, '(a)') 'DREAM - Fatal error!'
    write (*, '(a)') '  CR_NUM < 1.'
    stop
  end if

  if (gen_num < 2) then
    write (*, '(a)') ' '
    write (*, '(a)') 'DREAM - Fatal error!'
    write (*, '(a)') '  GEN_NUM < 2.'
    stop
  end if

  if (pair_num < 0) then
    write (*, '(a)') ' '
    write (*, '(a)') 'DREAM - Fatal error!'
    write (*, '(a)') '  PAIR_NUM < 0.'
    stop
  end if

  if (par_num < 1) then
    write (*, '(a)') ' '
    write (*, '(a)') 'DREAM - Fatal error!'
    write (*, '(a)') '  PAR_NUM < 1.'
    stop
  end if

  !IF (PGOPEN('/XSERVE') .LE. 0) STOP 'Error creating plot window'
  !CALL PGPAGE
  !CALL PGBBUF

  !  Get the problem parameter values.
  chain_filename = ''
  gr_filename = ''
  allocate (limits(1:2,1:par_num))
  restart_read_filename = ''
  restart_write_filename = ''

  call problem_value (chain_filename, gr_filename, gr_threshold, &
    jumpstep, limits, par_num, printstep, restart_read_filename, &
    restart_write_filename )

  !  Print the problem sizes and parameters.
  call input_print (chain_filename, chain_num, cr_num, gr_filename, &
    gr_threshold, jumpstep, limits, gen_num, pair_num, par_num, printstep, &
    restart_read_filename, restart_write_filename)

  !  Allocate memory.
  gr_num = gen_num / printstep

  allocate (fit(1:chain_num,1:gen_num))
  allocate (gr(1:par_num,1:gr_num))
  allocate (jumprate_table(1:par_num))
  allocate (z(1:par_num,1:chain_num,1:gen_num))

  !  Zero out memory.
  fit(1:chain_num,1:gen_num) = 0.0D+00
  gr(1:par_num,1:gr_num) = 0.0D+00
  jumprate_table(1:par_num) = 0.0D+00
  z(1:par_num,1:chain_num,1:gen_num) = 0.0D+00

  !  Set the jump rate table.
  call jumprate_table_init (jumprate_table, pair_num, par_num)

  call jumprate_table_print (jumprate_table, pair_num, par_num)

  !  Initialize the Gelman-Rubin data.
  call gr_init (gr, gr_conv, gr_count, gr_num, par_num)

  write (*, '(a)') ' '
  write (*, '(a)') 'GR_PRINT'
  write (*, '(a,l1)') '  GR_CONV  = ', gr_conv
  write (*, '(a,i6)') '  GR_COUNT = ', gr_count
  write (*, '(a,i6)') '  GR_NUM   = ', gr_num

  !  Set the first generation of the chains from restart data, or by sampling.
  if (0 < len_trim (restart_read_filename)) then
    call restart_read (chain_num, fit, gen_num, par_num, restart_read_filename, z)
  else
    call chain_init (chain_num, fit, gen_num, par_num, z)
  end if

  call chain_init_print (chain_num, fit, gen_num, par_num, restart_read_filename, z)

  !  Carry out the DREAM algorithm.
  call dream_algm (chain_num, cr_num, fit, gen_num, gr, gr_conv, &
    gr_count, gr_num, gr_threshold, jumprate_table, jumpstep, limits, &
    pair_num, par_num, printstep, z)

  !  Save Gelman-Rubin statistics to a file.
  if (0 < len_trim (gr_filename)) then
    call gr_write (gr, gr_filename, gr_num, par_num, printstep)
  end if

  !  Save parameter values for all chains at last generation.
  if (0 < len_trim (restart_write_filename)) then
    call restart_write (chain_num, fit, gen_num, par_num, restart_write_filename, z)
  end if

  !  Write each chain to a separate file.
  if (0 < len_trim (chain_filename)) then
    call chain_write (chain_filename, chain_num, fit, gen_num, par_num, z)
  end if

  !  Free memory.
  deallocate (fit)
  deallocate (gr)
  deallocate (jumprate_table)
  deallocate (limits)
  deallocate (z)
  call deallocate_problem()

  !  Terminate.
  write (*, '(a)') ' '
  write (*, '(a)') 'DREAM:'
  write (*, '(a)') '  Normal end of execution.'

  write (*, '(a)') ' '
  call timestamp ()

end subroutine main

subroutine chain_init (chain_num, fit, gen_num, par_num, z)
  ! starts Markov chains from a prior distribution.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Output, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Output, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  use problem1, only: sample_likelihood, prior_sample
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) i
  !real (kind = 8) sample_likelihood
  real (kind = 8) z(par_num,chain_num,gen_num)
  real (kind = 8), allocatable :: zp(:)

  allocate (zp(1:par_num))

  !$OMP PARALLEL
  !$OMP DO SCHEDULE(DYNAMIC,1)
  do i = 1, chain_num

    call prior_sample (par_num, zp)

    z(1:par_num,i,1) = zp(1:par_num)

    fit(i,1) = sample_likelihood (par_num, zp(1:par_num))

  end do
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate (zp)

  return
end subroutine chain_init

subroutine chain_init_print (chain_num, fit, gen_num, par_num, restart_read_filename, z)
  ! prints the initial values for Markov chains.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, character (len = *) RESTART_READ_FILENAME, the name of the file
  !    containing restart information, or '' if this is not a restart run.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  implicit none

  integer (kind = 4) gen_num
  integer (kind = 4) chain_num
  integer (kind = 4) par_num

  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) j
  character (len = *) restart_read_filename
  real (kind = 8) z(par_num,chain_num,gen_num)

  write (*, '(a)') ' '
  write (*, '(a)') 'CHAIN_INIT_PRINT'
  write (*, '(a)') '  Display initial values of Markov chains.'

  if (0 < len_trim (restart_read_filename)) then
    write (*, '(a)') '  Initialization from restart file "' &
      // trim (restart_read_filename) // '".'
  else
    write (*, '(a)') '  Initialization by sampling prior density.'
  end if

  do j = 1, chain_num
    write (*, '(a)') ' '
    write (*, '(a,i4)') '  Chain ', j
    write (*, '(a,g14.6)') '  Fitness ', fit(j,1)
    write (*, '(5g14.6)') z(1:par_num,j,1)    
  end do

  return
end subroutine chain_init_print

subroutine chain_outliers (chain_num, gen_index, gen_num, par_num, fit, z)
  ! identifies and modifies outlier chains during burn-in.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the index of the current generation.
  !    2 <= GEN_INDEX <= GEN_NUM.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input/output, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input/output, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov
  !    chain sample data.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  real (kind = 8), allocatable :: avg(:)
  real (kind = 8) avg_max
  real (kind = 8), allocatable :: avg_sorted(:)
  integer (kind = 4) best
  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) gen_index
  integer (kind = 4) i
  integer (kind = 4) ind1
  integer (kind = 4) ind3
  integer (kind = 4) j
  integer (kind = 4) klo
  integer (kind = 4) knum
  integer (kind = 4) outlier_num
  real (kind = 8) q1
  real (kind = 8) q3
  real (kind = 8) qr
  integer (kind = 4) r8_round_i4
  real (kind = 8) z(par_num,chain_num,gen_num)

  klo = gen_index / 2
  knum = gen_index + 1 - klo

  allocate (avg(1:chain_num))

  do j = 1, chain_num
    avg(j) = sum (fit(j,klo:gen_index)) / real (knum, kind = 8)
  end do 

  !  Set BEST to be the index of the chain with maximum average.
  best = 1
  avg_max = avg(1)
  do j = 2, chain_num
    if (avg_max < avg(j)) then
      best = j
      avg_max = avg(j)
    end if
  end do

  !  Determine the indices of the chains having averages 1/4 "above" 
  !  and "below" the average.
  allocate (avg_sorted(1:chain_num))

  call r8vec_copy (chain_num, avg, avg_sorted)

  call r8vec_sort_heap_a (chain_num, avg_sorted)

  ind1 = 1 + r8_round_i4 (0.25D+00 * real (chain_num, kind = 8))
  ind3 = 1 + r8_round_i4 (0.75D+00 * real (chain_num, kind = 8))

  q1 = avg_sorted(ind1)
  q3 = avg_sorted(ind3)
  qr = q3 - q1

  deallocate (avg_sorted)

  !  Identify outlier chains, and replace their later samples
  !  with values from the "best" chain.
  outlier_num = 0
  do j = 1, chain_num
    if (avg(j) < q1 - 2.0D+00 * qr) then
      outlier_num = outlier_num + 1
      z(1:par_num,j,gen_index) = z(1:par_num,best,gen_index)
      fit(j,klo:gen_index) = fit(best,klo:gen_index)
    end if
  end do

  if (0 < outlier_num) then

    write (*, '(a)') ' '
    write (*, '(a)') 'CHAIN_OUTLIERS:'
    write (*, '(a,i4,a,i4,a)') &
      '  At iteration ', gen_index, ' found ', outlier_num, ' outlier chains'
    !write (*, '(a)') '  whose indices appear below, and for which samples'
    !write (*, '(a)') &
    !  '  from the chain with the largest log likelihood function,'
    !write (*, '(a,i4,a)') '  index number ', best, '  will be substituted.'

    !do j = 1, chain_num
    !  if (avg(j) < q1 - 2.0D+00 * qr) then 
    !    write (*, '(2x,i4)') j
    !  end if
    !end do

  end if

  deallocate (avg)

  return
end subroutine chain_outliers

subroutine chain_write (chain_filename, chain_num, fit, gen_num, par_num, z)
  ! writes samples of each chain to separate files.
  !  Parameters:
  !
  !    Input, character (len = *) CHAIN_FILENAME, the "base" filename
  !    to be used for the chain files.  If this is the empty string '',
  !    then the chain files will not be written.  This name should 
  !    include a string of 0's which will be replaced by the chain 
  !    indices.  For example, "chain000.txt" would work as long as the
  !    number of chains was 1000 or less.
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  character (len = *) chain_filename
  integer (kind = 4) chain_unit
  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) flag
  integer (kind = 4) ind1
  integer (kind = 4) ind2
  integer (kind = 4) j
  integer (kind = 4) k
  real (kind = 8) z(par_num,chain_num,gen_num)

  !  Write parameter samples of all chains.
  write (*, '(a)') ' '
  write (*, '(a)') 'CHAIN_WRITE:'

  do j = 1, chain_num

    call filename_inc (chain_filename)

    call get_unit (chain_unit)

    open (unit = chain_unit, file = chain_filename, status = 'replace', &
      iostat = flag)

    if (flag /= 0) then
      write (*, '(a)') ' '
      write (*, '(a)') 'CHAIN_WRITE - Fatal error!'
      write (*, '(a)') &
        '  Could not open the file "' // trim (chain_filename) // '".'
      stop
    end if

    write (chain_unit, '(a,i2.2)') &
      'DREAM.F90:Parameters_and_log_likelihood_for_chain_#', j

    do k = 1, gen_num
      write (chain_unit, '(1x,i7,6x,es14.7,6x,1000(es14.7,2x))') &
        k, fit(j,k), z(1:par_num,j,k) 
    end do

    close (unit = chain_unit)

    write (*, '(a)') '  Created file "' // trim (chain_filename) // '".'

  end do

  return
end subroutine chain_write

subroutine cr_dis_update (chain_index, chain_num, cr_dis, cr_index, cr_num, &
  cr_ups, gen_index, gen_num, par_num, z) 
  ! updates the CR distance.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_INDEX, the index of the chain.
  !    1 <= CHAIN_INDEX <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input/output, real (kind = 8) CR_DIS(CR_NUM), the CR distances.
  !
  !    Input, integer (kind = 4) CR_INDEX, the index of the CR.
  !    1 <= CR_INDEX <= CR_NUM.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Input/output, integer (kind = 4) CR_UPS(CR_NUM), the number of updates 
  !    for each CR.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the index of the generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) cr_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  integer (kind = 4) chain_index
  real (kind = 8) cr_dis(cr_num)
  integer (kind = 4) cr_index
  integer (kind = 4) cr_ups(cr_num)
  integer (kind = 4) gen_index
  integer (kind = 4) i
  real (kind = 8) std(par_num) 
  real (kind = 8) z(par_num,chain_num,gen_num) 

  !  Compute the standard deviations.
  call std_compute (chain_num, gen_index, gen_num, par_num, z, std)

  !  Increment the update count.
  cr_ups(cr_index) = cr_ups(cr_index) + 1

  !  Update the CR distance.
  do i = 1, par_num 

    cr_dis(cr_index) = cr_dis(cr_index) &
      + ((z(i,chain_index,gen_index) - z(i,chain_index,gen_index-1)) &
      / std(i)) ** 2

  end do 

  return
end subroutine cr_dis_update

subroutine cr_index_choose (cr_index, cr_num, cr_prob)
  ! chooses a CR value. Index I is chosen with probability CR_PROB(I).
  !  Parameters:
  !
  !    Output, integer (kind = 4) CR_INDEX, the index of the CR.
  !    1 <= CR_INDEX <= CR_NUM.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Input, real (kind = 8) CR_PROB(CR_NUM), the probability of each CR.
  !
  implicit none

  integer (kind = 4) cr_num

  integer (kind = 4) cr_index
  real (kind = 8) cr_prob(cr_num)
  integer (kind = 4) i
  integer (kind = 4) n
  integer (kind = 4), allocatable :: tmp_index(:) 

  if (cr_num == 1) then

    cr_index = 1

  else
 
    allocate (tmp_index(1:cr_num))

    n = 1
    call i4vec_multinomial_sample (n, cr_prob, cr_num, tmp_index)

    do i = 1, cr_num
      if (tmp_index(i) == 1) then
        cr_index = i
        exit
      end if
    end do

    deallocate (tmp_index)

  end if

  return
end subroutine cr_index_choose

subroutine cr_init (cr, cr_dis, cr_num, cr_prob, cr_ups)
  ! initializes the crossover probability values.
  !  Parameters:
  !
  !    Output, real (kind = 8) CR(CR_NUM), the CR values.
  !
  !    Output, real (kind = 8) CR_DIS(CR_NUM), the CR distances.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Output, real (kind = 8) CR_PROB(CR_NUM), the probability of each CR.
  !
  !    Output, integer (kind = 4) CR_UPS(CR_NUM), the number of updates
  !    for each CR.
  !
  implicit none

  integer (kind = 4) cr_num

  real (kind = 8) cr(cr_num)
  real (kind = 8) cr_dis(cr_num)
  real (kind = 8) cr_prob(cr_num)
  integer (kind = 4) cr_ups(cr_num)
  integer (kind = 4) i

  do i = 1, cr_num
    cr(i) = real (i, kind = 8) / real (cr_num, kind = 8)  
  end do

  cr_dis(1:cr_num) = 1.0D+00
  cr_prob(1:cr_num) = 1.0D+00 / real (cr_num, kind = 8)
  cr_ups(1:cr_num) = 1

  return
end subroutine cr_init

subroutine cr_prob_update (cr_dis, cr_num, cr_prob, cr_ups)
  ! updates the CR probabilities.
  !  Parameters:
  !
  !    Input, real (kind = 8) CR_DIS(CR_NUM), the CR distances.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Output, real (kind = 8) CR_PROB(CR_NUM), the updated CR probabilities.
  !
  !    Input, integer (kind = 4) CR_UPS(CR_NUM), the number of updates 
  !    for each CR.
  !
  implicit none

  integer (kind = 4) cr_num

  real (kind = 8) cr_dis(cr_num)
  real (kind = 8) cr_prob(cr_num)
  real (kind = 8) cr_prob_sum
  integer (kind = 4) cr_ups(cr_num)
  integer (kind = 4) i

  do i = 1, cr_num
    cr_prob(i) = cr_dis(i) / real (cr_ups(i), kind = 8)
  end do

  cr_prob_sum = sum (cr_prob(1:cr_num))

  cr_prob(1:cr_num) = cr_prob(1:cr_num) / cr_prob_sum

  return
end subroutine cr_prob_update


subroutine diff_compute (chain_num, gen_index, gen_num, jump_dim, jump_num, &
  pair_num, par_num, r, z, diff) 
  ! computes the differential evolution.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the index of the current generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) JUMP_DIM(JUMP_NUM), the dimensions in which
  !    a jump is to be made.
  !
  !    Input, integer (kind = 4) JUMP_NUM, the number of dimensions in which
  !    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
  !
  !    Input, integer (kind = 4) PAIR_NUM, the number of pairs of 
  !    crossover chains.
  !    0 <= PAIR_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, integer (kind = 4) R(2,PAIR_NUM), pairs of chains used
  !    to compute differences.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  !    Output, real (kind = 8) DIFF(PAR_NUM), the vector of pair differences.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) jump_num
  integer (kind = 4) pair_num
  integer (kind = 4) par_num

  real (kind = 8) diff(par_num)
  integer (kind = 4) gen_index
  integer (kind = 4) i
  integer (kind = 4) j
  integer (kind = 4) jump_dim(jump_num)
  integer (kind = 4) k
  integer (kind = 4) pair
  integer (kind = 4) r(2,pair_num)
  real (kind = 8) z(par_num,chain_num,gen_num)

  !  Produce the difference of the pairs used for population evolution.
  diff(1:par_num) = 0.0D+00

  do pair = 1, pair_num
    do j = 1, jump_num
      k = jump_dim(j)
      diff(k) = diff(k) &
        + (z(k,r(1,pair),gen_index-1) - z(k,r(2,pair),gen_index-1))
    end do
  end do

  return
end subroutine diff_compute

subroutine dream_algm (chain_num, cr_num, fit, gen_num, gr, gr_conv, &
  gr_count, gr_num, gr_threshold, jumprate_table, jumpstep, limits, &
  pair_num, par_num, printstep, z)
  ! gets a candidate parameter sample.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Input, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, real (kind = 8) GR(PAR_NUM,GR_NUM), 
  !    the Gelman-Rubin R statistic.
  !
  !    Input/output, logical GR_CONV, the Gelman-Rubin convergence flag.
  !
  !    Input/output, integer (kind = 4) GR_COUNT, counts the number of 
  !    generations at which the Gelman-Rubin statistic has been computed.
  !
  !    Input, integer (kind = 4) GR_NUM, the number of times the Gelman-Rubin
  !    statistic may be computed.
  !
  !    Input, real (kind = 8) GR_THRESHOLD, the convergence tolerance for
  !    the Gelman-Rubin statistic.
  !
  !    Input, real (kind = 8) JUMPRATE_TABLE(PAR_NUM), the jumprate table.
  !
  !    Input, integer (kind = 4) JUMPSTEP, forces a "long jump" every
  !    JUMPSTEP generations.
  !
  !    Input, real (kind = 8) LIMITS(2,PAR_NUM), lower and upper bounds
  !    for each parameter.
  !
  !    Input, integer (kind = 4) PAIR_NUM, the number of pairs of 
  !    crossover chains.
  !    0 <= PAIR_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, integer (kind = 4) PRINTSTEP, the interval between generations on 
  !    which the Gelman-Rubin statistic will be computed and written to a file.
  !
  !    Output, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  !  Local parameters:
  !
  !    Local, integer (kind = 4) CHAIN_INDEX, the index of the current chain.
  !    1 <= CHAIN_INDEX <= CHAIN_NUM.
  !
  !    Local, real (kind = 8) CR(CR_NUM), the CR values.
  !
  !    Local, real (kind = 8) CR_DIS(CR_NUM), the CR distances.
  !
  !    Local, integer (kind = 4) CR_INDEX, the index of the selected CR value.
  !    1 <= CR_INDEX <= CR_NUM.
  !
  !    Local, real (kind = 8) CR_PROB(CR_NUM), the probability of each CR.
  !
  !    Local, real (kind = 8) CR_UPS(CR_NUM), the number of updates for each CR.
  !
  !    Local, integer (kind = 4) GEN_INDEX, the index of the current generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Local, real (kind = 8) ZP(PAR_NUM), a candidate sample.
  !
  !    Local, integer (kind = 4) ZP_ACCEPT, the number of candidates accepted.
  !
  !    Local, real (kind = 8) ZP_ACCEPT_RATE, the rate at which generated
  !    candidates were accepted.
  !
  !    Local, integer (kind = 4) ZP_COUNT, the number of candidates generated.
  !
  !    Local, real (kind = 8) ZP_RATIO, the Metropolis ratio for a candidate.
  !
  use problem1, only: sample_likelihood, prior_density
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) cr_num
  integer (kind = 4) gen_num
  integer (kind = 4) gr_num
  integer (kind = 4) par_num

  integer (kind = 4) chain_index
  real (kind = 8) cr(cr_num)
  real (kind = 8) cr_dis(cr_num)
  integer (kind = 4) cr_index
  real (kind = 8) cr_prob(cr_num)
  integer (kind = 4) cr_ups(cr_num)
  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) gen_index
  real (kind = 8) gr(par_num,gr_num)
  logical gr_conv
  integer (kind = 4) gr_count
  real (kind = 8) gr_threshold
  integer (kind = 4) ind1
  integer (kind = 4) ind2
  real (kind = 8) jumprate_table(par_num)
  integer (kind = 4) jumpstep
  real (kind = 8) limits(2,par_num)
  integer (kind = 4) pair_num
  real (kind = 8) pd1
  real (kind = 8) pd2
  integer (kind = 4) printstep
  !real (kind = 8) prior_density
  real (kind = 8) r
  real (kind = 8) r8_uniform_01_sample
  !real (kind = 8) sample_likelihood
  real (kind = 8) z(par_num,chain_num,gen_num)
  real (kind = 8) zp(par_num)
  integer (kind = 4) zp_accept
  real (kind = 8) zp_accept_rate
  integer (kind = 4) zp_count
  real (kind = 8) zp_fit
  real (kind = 8) zp_old(par_num)
  real (kind = 8) zp_old_fit
  real (kind = 8) zp_ratio

  zp_count = 0
  zp_accept = 0

  !  Initialize the CR values.
  call cr_init (cr, cr_dis, cr_num, cr_prob, cr_ups)

  do gen_index = 2, gen_num

    !$OMP PARALLEL
    !$OMP DO SCHEDULE(DYNAMIC,1)
    do chain_index = 1, chain_num
    
      !  Choose CR_INDEX, the index of a CR.
      call cr_index_choose (cr_index, cr_num, cr_prob)

      !  Generate a sample candidate ZP.
      call sample_candidate (chain_index, chain_num, cr, cr_index, cr_num, &
        gen_index, gen_num, jumprate_table, jumpstep, limits, pair_num, &
        par_num, z, zp)

      zp_count = zp_count + 1

      !  Compute the log likelihood function for ZP.
      zp_fit = sample_likelihood (par_num, zp)

      zp_old(1:par_num) = z(1:par_num,chain_index,gen_index-1)
      zp_old_fit = fit(chain_index,gen_index-1)

      !  Compute the Metropolis ratio for ZP versus ZP_OLD.
      pd1 = prior_density (par_num, zp)
      pd2 = prior_density (par_num, zp_old)

      zp_ratio = exp (&
        (zp_fit     + log (pd1)) - &
        (zp_old_fit + log (pd2)))

      zp_ratio = min (zp_ratio, 1.0D+00)

      !  Accept the candidate, or copy the value from the previous generation.
      r = r8_uniform_01_sample()

      if (r <= zp_ratio) then
        z(1:par_num,chain_index,gen_index) = zp(1:par_num)
        fit(chain_index,gen_index) = zp_fit
        zp_accept = zp_accept + 1
      else
        z(1:par_num,chain_index,gen_index) = zp_old(1:par_num)
        fit(chain_index,gen_index) = zp_old_fit
      end if

      !  Update the CR distance.
      if (.not. gr_conv) then
        if (1 < cr_num) then  
          call cr_dis_update (chain_index, chain_num, cr_dis, cr_index, &
            cr_num, cr_ups, gen_index, gen_num, par_num, z)
        end if
      end if

    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !!!!!!!!!!!!! ************* !!!!!!!!!!!!!!!!!!!
    !call one_chain_plot(z, gen_index, gen_num, par_num, chain_num)
    !!!!!!!!!!!!! ************* !!!!!!!!!!!!!!!!!!!

    !  Update the multinomial distribution of CR.
    if (.not. gr_conv) then
      if (1 < cr_num) then
        if (mod (gen_index, 10) == 0) then
          call cr_prob_update (cr_dis, cr_num, cr_prob, cr_ups)
        end if
      end if
    end if

    !  Every PRINTSTEP interval,
    !  compute the Gelman Rubin R statistic for this generation,
    !  and determine if convergence has occurred.
    if (mod (gen_index, printstep) == 0) then
      call gr_compute (chain_num, gen_index, gen_num, &
        gr, gr_conv, gr_count, gr_num, gr_threshold, par_num, z)

      write(*,'(i4)', advance='no') gen_index
      write(*,'(6f14.4)') z(1:par_num, 1, gen_index)
      !if (gr_conv) exit ! this was not here
    end if

    !  Check for outlier chains.
    if (.not. gr_conv) then ! only during burn-in
      if (mod (gen_index, 10) == 0) then
        call chain_outliers (chain_num, gen_index, gen_num, par_num, fit, z)
      end if
    end if

  end do        

  !  CALL PGEBUF

  !  Compute the acceptance rate.
  zp_accept_rate = real (zp_accept, kind = 8) / real (zp_count, kind = 8)

  write (*, '(a)') ' '
  write (*, '(a,g14.6)') '  The acceptance rate is ', zp_accept_rate

  return
end subroutine dream_algm

!subroutine one_chain_plot(z, gen_index, gen_num, par_num, chain_num)
!  use problem1, only: n, times, rvs, vel, last_vel
!  integer(kind=4) gen_index, gen_num
!  integer(kind=4) par_num, chain_num
!  real (kind=8) z(par_num,chain_num,gen_num)

!  call PGEX8(gen_num, 2, gen_index, z(1:par_num-1:5, 1, gen_index), z(3:par_num-1:5, 1, gen_index))

  !**********************      PLOT      ******************

    ! Define the Viewport
!    CALL PGSVP(0.1, 0.6, 0.45, 0.99)
    ! Define the Window
    !CALL PGSWIN(2449000., 2453000., -150.0, 150.0)
!    CALL PGSWIN(13000., 16000., -150.0, 150.0)
    ! Draw a box
!    CALL PGSCI(4) ! blue
!    CALL PGBOX ('BCTS', 0., 0, 'BCTSV', 0.0, 0)
    !CALL PGERAS
!    CALL PGSCI (6) ! magenta
!    CALL PGPT (n, REAL(times), REAL(rvs), 17)
!
!    CALL PGSCI (0) ! white
!    CALL PGPT (n, REAL(times), REAL(last_vel), 17)
!    CALL PGSCI (1) ! white
!    CALL PGPT (n, REAL(times), REAL(vel), 17)
!    last_vel = vel
!
!end subroutine one_chain_plot

subroutine filename_inc (filename)
  ! increments a partially numeric filename.
  !
  !  Discussion:
  !
  !    It is assumed that the digits in the name, whether scattered or
  !    connected, represent a number that is to be increased by 1 on
  !    each call.  If this number is all 9's on input, the output number
  !    is all 0's.  Non-numeric letters of the name are unaffected.
  !
  !    If the name is empty, then the routine stops.
  !
  !    If the name contains no digits, the empty string is returned.
  !
  !  Example:
  !
  !      Input            Output
  !      -----            ------
  !      'a7to11.txt'     'a7to12.txt'
  !      'a7to99.txt'     'a8to00.txt'
  !      'a9to99.txt'     'a0to00.txt'
  !      'cat.txt'        ' '
  !      ' '              STOP!
  !  Parameters:
  !
  !    Input/output, character (len = *) FILENAME.
  !    On input, a character string to be incremented.
  !    On output, the incremented string.
  !
  implicit none

  character c
  integer (kind = 4) change
  integer (kind = 4) digit
  character (len = *) filename
  integer (kind = 4) i
  integer (kind = 4) lens

  lens = len_trim (filename)

  if (lens <= 0) then
    write (*, '(a)') ' '
    write (*, '(a)') 'FILENAME_INC - Fatal error!'
    write (*, '(a)') '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = filename(i:i)

    if (lge (c, '0') .and. lle (c, '9')) then

      change = change + 1

      digit = ichar (c) - 48
      digit = digit + 1

      if (digit == 10) then
        digit = 0
      end if

      c = char (digit + 48)

      filename(i:i) = c

      if (c /= '0') then
        return
      end if

    end if

  end do

  !  No digits were found.  Return blank.
  if (change == 0) then
    filename = ' '
    return
  end if

  return
end subroutine filename_inc

subroutine get_unit (iunit)
  ! returns a free FORTRAN unit number.
  !
  !  Discussion:
  !
  !    A "free" FORTRAN unit number is a value between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5, 6 and 9, which
  !    are commonly reserved for console I/O).
  !
  !    Otherwise, IUNIT is a value between 1 and 99, representing a
  !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  !    are special, and will never return those values.
  !
  !  Parameters:
  !
  !    Output, integer (kind = 4) IUNIT, the free unit number.
  !
  implicit none

  integer (kind = 4) i
  integer (kind = 4) ios
  integer (kind = 4) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if (i /= 5 .and. i /= 6 .and. i /= 9) then

      inquire (unit = i, opened = lopen, iostat = ios)

      if (ios == 0) then
        if (.not. lopen) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end subroutine get_unit

subroutine gr_compute (chain_num, gen_index, gen_num, gr, gr_conv, gr_count, &
  gr_num, gr_threshold, par_num, z)
  ! computes the Gelman Rubin statistics R used to check convergence.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the index of the current generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Output, real (kind = 8) GR(PAR_NUM,GR_NUM), the Gelman-Rubin R statistic.
  !
  !    Output, logical GR_CONV, the Gelman-Rubin convergence flag.
  !
  !    Input/output, integer (kind = 4) GR_COUNT, counts the number of 
  !    generations at which the Gelman-Rubin statistic has been computed.
  !
  !    Input, integer (kind = 4) GR_NUM, the number of times the Gelman-Rubin
  !    statistic may be computed.
  !
  !    Input, real (kind = 8) GR_THRESHOLD, the convergence tolerance for the
  !    Gelman-Rubin statistic.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) gr_num
  integer (kind = 4) par_num

  real (kind = 8) b_var
  integer (kind = 4) chain_index
  integer (kind = 4) gen_index
  real (kind = 8) gr(par_num,gr_num)
  logical gr_conv
  integer (kind = 4) gr_count
  real (kind = 8) gr_threshold
  integer (kind = 4) ind0
  real (kind = 8) mean_all
  real (kind = 8) mean_chain(chain_num)
  integer (kind = 4) par_index
  real (kind = 8) rnd0
  real (kind = 8) s
  real (kind = 8) s_sum
  real (kind = 8) var
  real (kind = 8) w_var
  real (kind = 8) z(par_num,chain_num,gen_num)

  gr_count = gr_count + 1
  ind0 = gen_index / 2
  rnd0 = real (ind0, kind = 8)

  do par_index = 1, par_num

    do chain_index = 1, chain_num
      mean_chain(chain_index) = &
        sum (z(par_index,chain_index,ind0:gen_index)) / rnd0
    end do

    mean_all = sum (mean_chain(1:chain_num)) / real (chain_num, kind = 8)

    b_var = rnd0 &
      * sum ((mean_chain(1:chain_num) - mean_all)**2) &
      / real (chain_num - 1, kind = 8) 

    s_sum = 0.0D+00
    do chain_index = 1, chain_num
      s = sum ((z(par_index,chain_index,ind0:gen_index) &
        - mean_chain(chain_index))**2) 
      s_sum = s_sum + s
    end do
    s_sum = s_sum / (rnd0 - 1.0D+00)

    w_var = s_sum / real (chain_num, kind = 8)

    var = ((rnd0 - 1.0D+00) * w_var + b_var) / rnd0 

    gr(par_index,gr_count) = sqrt (var / w_var)

  end do

  !  Set the convergence flag.
  gr_conv = .true.

  do par_index = 1, par_num
    if (gr_threshold < gr(par_index,gr_count)) then
      gr_conv = .false.
      exit
    end if
  end do

  if (gr_conv) then
    write (*, '(a)') ' '
    write (*, '(a)') 'GR_COMPUTE:'
    write (*, '(a,i6)') '  GR convergence at iteration: ', gen_index
  end if

  return
end subroutine gr_compute

subroutine gr_init (gr, gr_conv, gr_count, gr_num, par_num)
  ! initializes Gelman-Rubin variables.
  !  Parameters:
  !
  !    Output, real (kind = 8) GR(PAR_NUM,GR_NUM), the Gelman-Rubin statistic.
  !
  !    Output, logical GR_CONV, the convergence flag.
  !
  !    Output, integer (kind = 4) GR_COUNT, counts the number of generations
  !    at which the Gelman-Rubin statistic has been computed.
  !
  !    Input, integer (kind = 4) GR_NUM, the number of times the Gelman-Rubin
  !    statistic may be computed.
  !
  !    Input, integer (kind = 4) PAR_NUM, the number of parameters.
  !    1 <= PAR_NUM.
  !
  implicit none

  integer (kind = 4) gr_num
  integer (kind = 4) par_num

  real (kind = 8) gr(par_num,gr_num)
  logical gr_conv
  integer (kind = 4) gr_count

  gr(1:par_num,1:gr_num) = 0.0D+00
  gr_conv = .false.
  gr_count = 0

  return
end subroutine gr_init

subroutine gr_write (gr, gr_filename, gr_num, par_num, printstep)
  ! writes the Gelman-Rubin R statistics into a file.
  !  Parameters:
  !
  !    Input, real (kind = 8) GR(PAR_NUM,GR_NUM), the Gelman-Rubin R statistic.
  !
  !    Input, character (len = *) GR_FILENAME, the name of the file
  !    in which values of the Gelman-Rubin statistic will be recorded,
  !    or '' if no such file is to be created.
  !
  !    Input, integer (kind = 4) GR_NUM, the number of times the Gelman-Rubin
  !    statistic may be computed.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, integer (kind = 4) PRINTSTEP, the interval between generations on 
  !    which the Gelman-Rubin statistic will be computed and written to a file.
  !
  implicit none

  integer (kind = 4) gr_num
  integer (kind = 4) par_num

  integer (kind = 4) flag
  real (kind = 8) gr(par_num,gr_num)
  character (len = *) gr_filename
  integer (kind = 4) gr_unit
  integer (kind = 4) j
  integer (kind = 4) printstep

  call get_unit (gr_unit)

  open (unit = gr_unit, file = gr_filename, status = 'replace', &
    iostat = flag)

  if (flag /= 0) then
    write (*, '(a)') ' '
    write (*, '(a)') 'GR_WRITE - Fatal error!'
    write (*, '(a)') &
      '  Could not open the file "' // trim (gr_filename) // '".'
    stop
  end if

  write (gr_unit, '(a)') &
    'DREAM.F90:Monitored_parameter_interchains_Gelman_Rubin_R_statistic'

  do j = 1, gr_num
    write (gr_unit, '(1x,i9,6x,1000(f14.4,2x))') &
      printstep * j, gr(1:par_num,j)
  end do

  close (unit = gr_unit)

  write (*, '(a)') ' '
  write (*, '(a)') 'GR_WRITE:'
  write (*, '(a)') '  Created the file "' // trim (gr_filename) // '".'

  return
end subroutine gr_write

subroutine i4mat_print (m, n, a, title)
  ! prints an I4MAT.
  !    An I4MAT is a rectangular array of I4's.
  !  Parameters:
  !
  !    Input, integer (kind = 4) M, the number of rows in A.
  !
  !    Input, integer (kind = 4) N, the number of columns in A.
  !
  !    Input, integer (kind = 4) A(M,N), the matrix to be printed.
  !
  !    Input, character (len = *) TITLE, a title.
  !
  implicit none

  integer (kind = 4) m
  integer (kind = 4) n

  integer (kind = 4) a(m,n)
  integer (kind = 4) ihi
  integer (kind = 4) ilo
  integer (kind = 4) jhi
  integer (kind = 4) jlo
  character (len = *) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some (m, n, a, ilo, jlo, ihi, jhi, title)

  return
end subroutine i4mat_print

subroutine i4mat_print_some (m, n, a, ilo, jlo, ihi, jhi, title)
  ! prints some of an I4MAT.
  !    An I4MAT is a rectangular array of I4's.
  !  Parameters:
  !
  !    Input, integer (kind = 4) M, N, the number of rows and columns.
  !
  !    Input, integer (kind = 4) A(M,N), an M by N matrix to be printed.
  !
  !    Input, integer (kind = 4) ILO, JLO, the first row and column to print.
  !
  !    Input, integer (kind = 4) IHI, JHI, the last row and column to print.
  !
  !    Input, character (len = *) TITLE, a title.
  !
  implicit none

  integer (kind = 4), parameter :: incx = 10
  integer (kind = 4) m
  integer (kind = 4) n

  integer (kind = 4) a(m,n)
  character (len = 8)  ctemp(incx)
  integer (kind = 4) i
  integer (kind = 4) i2hi
  integer (kind = 4) i2lo
  integer (kind = 4) ihi
  integer (kind = 4) ilo
  integer (kind = 4) inc
  integer (kind = 4) j
  integer (kind = 4) j2
  integer (kind = 4) j2hi
  integer (kind = 4) j2lo
  integer (kind = 4) jhi
  integer (kind = 4) jlo
  character (len = *) title

  write (*, '(a)') ' '
  write (*, '(a)') trim (title)

  if (m <= 0 .or. n <= 0) then
    write (*, '(a)') ' '
    write (*, '(a)') '  (None)'
    return
  end if

  do j2lo = max (jlo, 1), min (jhi, n), incx

    j2hi = j2lo + incx - 1
    j2hi = min (j2hi, n)
    j2hi = min (j2hi, jhi)

    inc = j2hi + 1 - j2lo

    write (*, '(a)') ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write (ctemp(j2), '(i8)') j
    end do

    write (*, '(''  Col '',10a8)') ctemp(1:inc)
    write (*, '(a)') '  Row'
    write (*, '(a)') ' '

    i2lo = max (ilo, 1)
    i2hi = min (ihi, m)

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write (ctemp(j2), '(i8)') a(i,j)

      end do

      write (*, '(i5,a,10a8)') i, ':', (ctemp(j), j = 1, inc)

    end do

  end do

  return
end subroutine i4mat_print_some

subroutine input_print (chain_filename, chain_num, cr_num, gr_filename, &
  gr_threshold, jumpstep, limits, gen_num, pair_num, par_num, printstep, &
  restart_read_filename, restart_write_filename)
  ! prints the data from the input file.
  !  Parameters:
  !
  !    Input, character (len = *) CHAIN_FILENAME, the "base" filename
  !    to be used for the chain files.  If this is the empty string '',
  !    then the chain files will not be written.  This name should 
  !    include a string of 0's which will be replaced by the chain 
  !    indices.  For example, "chain000.txt" would work as long as the
  !    number of chains was 1000 or less.
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Input, character (len = *) GR_FILENAME, the name of the file
  !    in which values of the Gelman-Rubin statistic will be recorded,
  !    or '' if no such file is to be created.
  !
  !    Input, real (kind = 8) GR_THRESHOLD, the convergence tolerance for the
  !    Gelman-Rubin statistic.
  !
  !    Input, integer (kind = 4) JUMPSTEP, forces a "long jump" every
  !    JUMPSTEP generations.
  !
  !    Input, real (kind = 8) LIMITS(2,PAR_NUM), lower and upper limits
  !    for each parameter.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAIR_NUM, the number of pairs of 
  !    crossover chains.
  !    0 <= PAIR_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, integer (kind = 4) PRINTSTEP, the interval between generations on 
  !    which the Gelman-Rubin statistic will be computed and written to a file.
  !
  !    Input, character (len = *) RESTART_READ_FILENAME, the name of the file
  !    containing restart information, or '' if this is not a restart run.
  !
  !    Input, character (len = *) RESTART_WRITE_FILENAME, the name of the file
  !    to be written, containing restart information, or '' if a restart file 
  !    is not to be written.
  !
  implicit none

  integer (kind = 4) par_num

  character (len = *) chain_filename
  integer (kind = 4) chain_num
  integer (kind = 4) cr_num
  character (len = *) gr_filename
  real (kind = 8) gr_threshold
  integer (kind = 4) j
  integer (kind = 4) jumpstep
  real (kind = 8) limits(2,par_num)
  integer (kind = 4) gen_num
  integer (kind = 4) pair_num
  integer (kind = 4) printstep
  character (len = *) restart_read_filename
  character (len = *) restart_write_filename

  write (*, '(a)') ' '
  write (*, '(a)') 'INPUT_PRINT:'
  write (*, '(a)') ' '
  write (*, '(a)') '  Number of parameters'
  write (*, '(a,i6)') '  PAR_NUM = ', par_num 
  write (*, '(a)') ' '
  write (*, '(a)') '  Lower and upper limits for each parameter:'
  write (*, '(a)') ' '
  write (*, '(a)') '   Index       Lower           Upper'
  write (*, '(a)') ' '
  do j = 1, par_num
    write (*, '(2x,i6,2x,g14.6,2x,g14.6)') j, limits(1:2,j)
  end do
  write (*, '(a)') ' '
  write (*, '(a)') '  Number of generations:'
  write (*, '(a,i6)') '  GEN_NUM = ', gen_num
  write (*, '(a)') ' '
  write (*, '(a)') '  Number of simultaneous chains:'
  write (*, '(a,i6)') '  CHAIN_NUM = ', chain_num
  write (*, '(a)') ' '
  write (*, '(a)') '  Chain filename (base):'
  if (len_trim (chain_filename) == 0) then
    write (*, '(a)') '  CHAIN_FILENAME = "(None)".'
  else
    write (*, '(a)') '  CHAIN_FILENAME = "' // trim (chain_filename) // '".'
  end if

  write (*, '(a)') ' '
  write (*, '(a)') '  Number of pairs of chains for crossover:'
  write (*, '(a,i6)') '  PAIR_NUM = ', pair_num
  write (*, '(a)') ' '
  write (*, '(a)') '  Number of crossover values:'
  write (*, '(a,i6)') '  CR_NUM = ', cr_num
  write (*, '(a)') ' '
  write (*, '(a)') '  Number of steps til a long jump:'
  write (*, '(a,i6)') '  JUMPSTEP = ', jumpstep
  write (*, '(a)') ' '
  write (*, '(a)') '  Interval between Gelman-Rubin computations:'
  write (*, '(a,i6)') '  PRINTSTEP = ', printstep
  write (*, '(a)') ' '
  write (*, '(a)') '  Gelman-Rubin output filename:'
  if (len_trim (gr_filename) == 0) then
    write (*, '(a)') '  GR_FILENAME = "(None)".'
  else
    write (*, '(a)') '  GR_FILENAME = "' // trim (gr_filename) // '".'
  end if
  write (*, '(a)') ' '
  write (*, '(a)') '  Gelman-Rubin convergence tolerance:'
  write (*, '(a,g14.6)') '  GR_THRESHOLD = ', gr_threshold
  write (*, '(a)') ' '
  write (*, '(a)') '  Restart read filename:'
  if (len_trim (restart_read_filename) == 0) then
    write (*, '(a)') '  RESTART_READ_FILENAME = "(None)".'
  else
    write (*, '(a)') '  RESTART_READ_FILENAME = "' // trim (restart_read_filename) // '".'
  end if
  write (*, '(a)') ' '
  write (*, '(a)') '  Restart write filename:'
  if (len_trim (restart_write_filename) == 0) then
    write (*, '(a)') '  RESTART_WRITE_FILENAME = "(None)".'
  else
    write (*, '(a)') '  RESTART_WRITE_FILENAME = "' // trim (restart_write_filename) // '".'
  end if

  return
end subroutine input_print

subroutine jumprate_choose (cr, cr_index, cr_num, gen_index, jump_dim, &
  jump_num, jumprate, jumprate_table, jumpstep, par_num)
  ! chooses a jump rate from the jump rate table.
  !  Parameters:
  !
  !    Input, real (kind = 8) CR(CR_NUM), the CR values.
  !
  !    Input, integer (kind = 4) CR_INDEX, the index of the CR.
  !    1 <= CR_INDEX <= CR_NUM.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the current generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Output, integer (kind = 4) JUMP_DIM(PAR_NUM), the indexes of the
  !    parameters to be updated.
  !
  !    Output, integer (kind = 4) JUMP_NUM, the number of dimensions in which
  !    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
  !
  !    Output, real (kind = 8) JUMPRATE, the jump rate.
  !
  !    Input, real (kind = 8) JUMPRATE_TABLE(PAR_NUM), the jump rate table.
  !
  !    Input, integer (kind = 4) JUMPSTEP, forces a "long jump" every
  !    JUMPSTEP generations.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  implicit none

  integer (kind = 4) cr_num
  integer (kind = 4) par_num

  real (kind = 8) cr(cr_num)
  integer (kind = 4) cr_index
  integer (kind = 4) gen_index
  integer (kind = 4) i 
  integer (kind = 4) jump_dim(par_num)
  integer (kind = 4) jump_num
  real (kind = 8) jumprate
  real (kind = 8) jumprate_table(par_num)
  integer (kind = 4) jumpstep
  real (kind = 8) r
  real (kind = 8) r8_uniform_01_sample

  !  Determine the dimensions that will be updated.
  jump_num = 0
  jump_dim(1:par_num) = 0

  do i = 1, par_num

    r = r8_uniform_01_sample ()

    if (1.0D+00 - cr(cr_index) < r) then
      jump_num = jump_num + 1
      jump_dim(jump_num) = i
    end if

  end do

  !  Calculate the general jump rate.
  if (jump_num == 0) then
    jumprate = 0.0D+00
  else
    jumprate = jumprate_table(jump_num)
  end if

  !  If parameter dimension is 1, 2, or 3, fix the jump rate to 0.6.
  if (par_num <= 3) then
    jumprate = 0.6D+00
  end if

  !  Determine if a long jump is forced.
  if (mod (gen_index - 1, jumpstep) == 0) then
    jumprate = 0.98D+00
  end if

  return
end subroutine jumprate_choose

subroutine jumprate_table_init (jumprate_table, pair_num, par_num)
  ! initializes the jump rate table.
  !  Parameters:
  !
  !    Output, real (kind = 8) JUMPRATE_TABLE(PAR_NUM), the jumprate table.
  !
  !    Input, integer (kind = 4) PAIR_NUM, the number of pairs of 
  !    crossover chains.
  !    0 <= PAIR_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  implicit none

  integer (kind = 4) par_num

  real (kind = 8) c
  integer (kind = 4) i
  real (kind = 8) jumprate_table(par_num)
  integer (kind = 4) pair_num

  c = 2.38D+00 / sqrt (real (2 * pair_num, kind = 8))

  do i = 1, par_num
    jumprate_table(i) = c / sqrt (real (i, kind = 8))
  end do

  return
end subroutine jumprate_table_init

subroutine jumprate_table_print (jumprate_table, pair_num, par_num)
  ! prints the jump rate table.
  !  Parameters:
  !
  !    Input, real (kind = 8) JUMPRATE_TABLE(PAR_NUM), the jumprate table.
  !
  !    Input, integer (kind = 4) PAIR_NUM, the number of pairs of 
  !    crossover chains.
  !    0 <= PAIR_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  implicit none

  integer (kind = 4) par_num

  integer (kind = 4) i
  real (kind = 8) jumprate_table(par_num)
  integer (kind = 4) pair_num
 
  write (*, '(a)') ' '
  write (*, '(a)') 'JUMPRATE_TABLE_PRINT'
  write (*, '(a)') ' '
  write (*, '(a)') '   I    Jumprate'
  write (*, '(a)') ' '
  do i = 1, par_num
    write (*, '(2x,i2,2x,g14.6)') i, jumprate_table(i)
  end do

  return
end subroutine jumprate_table_print

function r8_round_i4 (x)
  ! sets an R8 to the nearest integral value, returning an I4
  !
  !  Example:
  !
  !        X        R8_ROUND_I4
  !
  !      1.3         1
  !      1.4         1
  !      1.5         1 or 2
  !      1.6         2
  !      0.0         0
  !     -0.7        -1
  !     -1.1        -1
  !     -1.6        -2
  !
  !  Discussion:
  !
  !    In FORTRAN90, we rely on the fact that, for positive X, int (X)
  !    is the "floor" function, returning the largest integer less than
  !    or equal to X.
  !
  !  Parameters:
  !
  !    Input, real (kind = 8) X, the value.
  !
  !    Output, integer (kind = 4) R8_ROUND_I4, the rounded value.
  !
  implicit none

  integer (kind = 4) r8_round_i4
  integer (kind = 4) value
  real (kind = 8) x

  if (x < 0.0D+00) then
    value = - int (- x + 0.5D+00)
  else
    value =   int (+ x + 0.5D+00)
  end if

  r8_round_i4 = value

  return
end function r8_round_i4

subroutine r8vec_copy (n, a1, a2)
  ! copies an R8VEC.
  !  Parameters:
  !
  !    Input, integer (kind = 4) N, the length of the vectors.
  !
  !    Input, real (kind = 8) A1(N), the vector to be copied.
  !
  !    Output, real (kind = 8) A2(N), a copy of A1.
  !
  implicit none

  integer (kind = 4) n

  real (kind = 8) a1(n)
  real (kind = 8) a2(n)

  a2(1:n) = a1(1:n)

  return
end subroutine r8vec_copy

subroutine r8vec_heap_d (n, a)
  ! reorders an R8VEC into an descending heap.
  !
  !  Discussion:
  !
  !    A descending heap is an array A with the property that, for every index J,
  !    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
  !    2*J and 2*J+1 are legal).
  !
  !                  A(1)
  !                /      \
  !            A(2)         A(3)
  !          /     \        /  \
  !      A(4)       A(5)  A(6) A(7)
  !      /  \       /   \
  !    A(8) A(9) A(10) A(11)
  !
  !  Parameters:
  !
  !    Input, integer (kind = 4) N, the size of the input array.
  !
  !    Input/output, real (kind = 8) A(N).
  !    On input, an unsorted array.
  !    On output, the array has been reordered into a heap.
  !
  implicit none

  integer (kind = 4) n

  real (kind = 8) a(n)
  integer (kind = 4) i
  integer (kind = 4) ifree
  real (kind = 8) key
  integer (kind = 4) m

  !  Only nodes N/2 down to 1 can be "parent" nodes.
  do i = n / 2, 1, -1

    !  Copy the value out of the parent node.
    !  Position IFREE is now "open".
    key = a(i)
    ifree = i

    do

      !  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
      !  IFREE.  (One or both may not exist because they exceed N.)
      m = 2 * ifree

      !  Does the first position exist?
      if (n < m) then
        exit
      end if

      !  Does the second position exist?
      if (m + 1 <= n) then

        !  If both positions exist, take the larger of the two values,
        !  and update M if necessary.
        if (a(m) < a(m+1)) then
          m = m + 1
        end if

      end if

      !  If the large descendant is larger than KEY, move it up,
      !  and update IFREE, the location of the free position, and
      !  consider the descendants of THIS position.
      if (a(m) <= key) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do

    !  Once there is no more shifting to do, KEY moves into the free spot IFREE.
    a(ifree) = key

  end do

  return
end subroutine r8vec_heap_d

subroutine r8vec_sort_heap_a (n, a)
  ! ascending sorts an R8VEC using heap sort.
  !  Parameters:
  !
  !    Input, integer (kind = 4) N, the number of entries in the array.
  !
  !    Input/output, real (kind = 8) A(N).
  !    On input, the array to be sorted;
  !    On output, the array has been sorted.
  !
  implicit none

  integer (kind = 4) n

  real (kind = 8) a(n)
  integer (kind = 4) n1
  real (kind = 8) temp

  if (n <= 1) then
    return
  end if

  !  1: Put A into descending heap form.
  call r8vec_heap_d (n, a)

  !  2: Sort A.
  !  The largest object in the heap is in A(1).
  !  Move it to position A(N).
  temp = a(1)
  a(1) = a(n)
  a(n) = temp

  !  Consider the diminished heap of size N1.
  do n1 = n - 1, 2, -1

    !  Restore the heap structure of A(1) through A(N1).
    call r8vec_heap_d (n1, a)

    !  Take the largest object from A(1) and move it to A(N1).
    temp = a(1)
    a(1) = a(n1)
    a(n1) = temp

  end do

  return
end subroutine r8vec_sort_heap_a

subroutine restart_read (chain_num, fit, gen_num, par_num, restart_read_filename, z)
  ! reads parameter sample data from a restart file.
  ! Only a single generation (presumably the last one) was written to the file.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Output, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, character (len = *) RESTART_READ_FILENAME, the name of 
  !    the restart file.
  !
  !    Output, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) flag
  integer (kind = 4) j
  integer (kind = 4) k
  character (len = *) restart_read_filename
  integer (kind = 4) restart_unit
  real (kind = 8) z(par_num,chain_num,gen_num)

  call get_unit (restart_unit)

  open (restart_unit, file = restart_read_filename, status = 'old', &
    iostat = flag)

  if (flag /= 0) then
    write (*, '(a)') ' '
    write (*, '(a)') 'RESTART_READ - Fatal error!'
    write (*, '(a)') &
      '  Could not open the file "' // trim (restart_read_filename) // '".'
    stop
  end if

  read (restart_unit, *)
  do j = 1, chain_num
    read (restart_unit, *) k, fit(j,1), z(1:par_num,j,1)
  end do

  close (restart_unit)

  return
end subroutine restart_read

subroutine restart_write (chain_num, fit, gen_num, par_num, restart_write_filename, z)
  ! writes a restart file.
  ! Only data for the final generation is written.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, real (kind = 8) FIT(CHAIN_NUM,GEN_NUM), the likelihood of
  !    each sample.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, character (len = *) RESTART_WRITE_FILENAME, the name of the 
  !    restart file.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  real (kind = 8) fit(chain_num,gen_num)
  integer (kind = 4) flag
  integer (kind = 4) i
  integer (kind = 4) restart_unit
  character (len = *) restart_write_filename
  real (kind = 8) z(par_num,chain_num,gen_num)
   
  call get_unit (restart_unit)
      
  open (unit = restart_unit, file = restart_write_filename, &
    status = 'replace', iostat = flag)

  if (flag /= 0) then
    write (*, '(a)') ' '
    write (*, '(a)') 'RESTART_WRITE - Fatal error!'
    write (*, '(a)') &
      '  Could not open the file "' // trim (restart_write_filename) // '".'
    stop
  end if

  write (restart_unit, '(a)') 'DREAM.F90:Parameter_values_for_restart'

  do i = 1, chain_num
    write (restart_unit, '(i8,7x,es14.7,6x,1000(es14.7,2x))') &
      i, fit(i,gen_num), z(1:par_num,i,gen_num)
  end do
      
  close (restart_unit)

  write (*, '(a)') ' '
  write (*, '(a)') 'RESTART_WRITE:'
  write (*, '(a)') '  Created restart file "' &
    // trim (restart_write_filename) // '".'

  return
end subroutine restart_write

subroutine sample_candidate (chain_index, chain_num, cr, cr_index, cr_num, &
  gen_index, gen_num, jumprate_table, jumpstep, limits, pair_num, par_num, &
  z, zp)
  ! generates candidate parameter samples.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_INDEX, the chain index.
  !    1 <= CHAIN_INDEX <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, real (kind = 8) CR(CR_NUM), the CR values.
  !
  !    Input, integer (kind = 4) CR_INDEX, the index of the chosen CR value.
  !    1 <= CR_INDEX <= CR_NUM.
  !
  !    Input, integer (kind = 4) CR_NUM, the total number of CR values.
  !    1 <= CR_NUM.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the current generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, real (kind = 8) JUMPRATE_TABLE(PAR_NUM), the jumprate table.
  !
  !    Input, integer (kind = 4) JUMPSTEP, forces a "long jump" every
  !    JUMPSTEP generations.
  !
  !    Input, real (kind = 8) LIMITS(2,PAR_NUM), limits for the parameters.
  !
  !    Input, integer (kind = 4) PAIR_NUM, the number of pairs of 
  !    crossover chains.
  !    0 <= PAIR_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  !    Output, real (kind = 8) ZP(PAR_NUM), a candidate parameter sample.
  !
  !  Local parameters:
  !
  !    Input, integer (kind = 4) JUMP_DIM(JUMP_NUM), the dimensions in which
  !    a jump is to be made.
  !
  !    Local, integer JUMP_NUM, the number of dimensions in which
  !    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
  !
  !    Local, real (kind = 8) JUMPRATE, the jump rate.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) cr_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  real (kind = 8) av
  real (kind = 8) b
  integer (kind = 4) chain_index
  real (kind = 8) cr(cr_num)
  integer (kind = 4) cr_index
  real (kind = 8), allocatable :: diff(:)
  real (kind = 8), allocatable :: eps(:)
  integer (kind = 4) gen_index
  integer (kind = 4) i
  integer (kind = 4) jump_dim(par_num)
  integer (kind = 4) jump_num
  real (kind = 8) jumprate
  real (kind = 8) jumprate_table(par_num)
  integer (kind = 4) jumpstep
  real (kind = 8) limits(2,par_num)
  real (kind = 8), allocatable :: noise_e(:)
  integer (kind = 4) pair(2)
  integer (kind = 4) pair_num
  integer (kind = 4), allocatable :: r(:,:)
  real (kind = 8) r2(2)
  real (kind = 8) r8_normal_sample
  real (kind = 8) r8_uniform_01_sample
  real (kind = 8) sd
  real (kind = 8) z(par_num,chain_num,gen_num)
  real (kind = 8) zp(par_num)

  !  Used to calculate E following a uniform distribution on (-B,+B).
  !  if B is zero, the noise term is suppressed.
  !  the value 0.05 is suggested in Vrugt et al.
  b = 0.05D+00

  !  Pick pairs of other chains for crossover.
  allocate (r(1:2,1:pair_num))

  do i = 1, pair_num

    do

      r2(1) = r8_uniform_01_sample ()
      r2(2) = r8_uniform_01_sample ()

      pair(1:2) = int (r2(1:2) * real (chain_num, kind = 8)) + 1

      if (pair(1) /= pair(2) .and. &
           pair(1) /= chain_index .and. &
           pair(2) /= chain_index) then
        exit
      end if

    end do

    r(1:2,i) = pair(1:2)

  end do

  !  Determine the jump rate.
  call jumprate_choose (cr, cr_index, cr_num, gen_index, &
    jump_dim, jump_num, jumprate, jumprate_table, jumpstep, par_num)

  !  Calculate E in equation 4 of Vrugt.
  allocate (noise_e(1:par_num))

  do i = 1, par_num
    noise_e(i) = b * (2.0D+00 * r8_uniform_01_sample () - 1.0D+00)
  end do

  !  Get epsilon value from multinormal distribution                      
  allocate (eps(1:par_num))

  av = 0.0D+00
  sd = 1.0D-06  ! value suggested by Vrugt et al.
  do i = 1, par_num
    eps(i) = r8_normal_sample (av, sd)
  end do

  !  Generate the candidate sample ZP based on equation 4 of Vrugt.
  allocate (diff(1:par_num))

  call diff_compute (chain_num, gen_index, gen_num, jump_dim, jump_num, &
    pair_num, par_num, r, z, diff)

  zp(1:par_num) = z(1:par_num,chain_index,gen_index-1)

  zp(1:par_num) = zp(1:par_num) &
    + (1.0D+00 + noise_e(1:par_num)) * jumprate * diff(1:par_num) &
    + eps(1:par_num)

  !  Enforce limits on the sample ZP.
  call sample_limits (limits, par_num, zp)

  deallocate (diff)
  deallocate (eps)
  deallocate (noise_e)
  deallocate (r)

  return
end subroutine sample_candidate

subroutine sample_limits (limits, par_num, zp)
  ! enforces limits on a sample variable.
  !  Parameters:
  !
  !    Input, real (kind = 8) LIMITS(2,PAR_NUM), the parameter limits.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input/output, real (kind = 8) ZP(PAR_NUM), a variable, whose entries,
  !    if necessary, will be "folded" so that they lie within the limits.
  !
  implicit none
   
  integer (kind = 4) par_num

  integer (kind = 4) i
  real (kind = 8) limits(2,par_num)
  real (kind = 8) w
  real (kind = 8) zp(par_num)

  do i = 1, par_num

    w = limits(2,i) - limits(1,i)

    if (w == 0.0D+00) then

      zp(i) = limits(1,i)

    else if (w < 0.0D+00) then

      write (*, '(a)') ' '
      write (*, '(a)') 'SAMPLE_LIMITS - Fatal error!'
      write (*, '(a)') '  LIMITS(2,I) < LIMITS(1,I).'
      stop

    else

      do while (zp(i) < limits(1,i))
        zp(i) = zp(i) + w
      end do

      do while (limits(2,i) < zp(i))
        zp(i) = zp(i) - w
      end do

    end if

  end do

  return
end subroutine sample_limits

subroutine std_compute (chain_num, gen_index, gen_num, par_num, z, std)
  ! computes the current standard deviations, for each parameter.
  ! The computation encompasses all chains and generations up to the
  ! current ones.
  !  Parameters:
  !
  !    Input, integer (kind = 4) CHAIN_NUM, the total number of chains.
  !    3 <= CHAIN_NUM.
  !
  !    Input, integer (kind = 4) GEN_INDEX, the current generation.
  !    1 <= GEN_INDEX <= GEN_NUM.
  !
  !    Input, integer (kind = 4) GEN_NUM, the total number of generations.
  !    2 <= GEN_NUM.
  !
  !    Input, integer (kind = 4) PAR_NUM, the total number of parameters.
  !    1 <= PAR_NUM.
  !
  !    Input, real (kind = 8) Z(PAR_NUM,CHAIN_NUM,GEN_NUM), the Markov chain 
  !    sample data.
  !
  !    Output, real (kind = 8) STD(PAR_NUM), the standard deviations.
  !
  implicit none

  integer (kind = 4) chain_num
  integer (kind = 4) gen_num
  integer (kind = 4) par_num

  integer (kind = 4) gen_index
  integer (kind = 4) i
  integer (kind = 4) j
  integer (kind = 4) k
  real (kind = 8) mean
  real (kind = 8) std(par_num)
  real (kind = 8) z(par_num,chain_num,gen_num)

  do i = 1, par_num

    mean = 0.0D+00
    do k = 1, gen_index
      do j = 1, chain_num
        mean = mean + z(i,j,k)
      end do
    end do
    mean = mean / real (chain_num, kind = 8) / real (gen_index, kind = 8)

    std(i) = 0.0D+00
    do k = 1, gen_index
      do j = 1, chain_num
        std(i) = std(i) + (z(i,j,k) - mean) ** 2
      end do
    end do

    std(i) = sqrt (std(i) / real (chain_num * gen_index - 1, kind = 8))   

  end do

  return
end subroutine std_compute

