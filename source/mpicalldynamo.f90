!*****************************************************
program magnetizer
  use mpi
  use dynamo
  use input_params
  use global_input_parameters
  use IO
  use messages

  implicit none

  integer :: igal, nmygals, igal_first, igal_last, igal_finished
  integer, parameter :: master_rank = 0
  integer, allocatable, dimension(:) :: mygals, allgals
  character(len=100) :: command_argument
  integer :: i, j, iproc
  logical :: lstop, error
  logical :: lsingle_galaxy_mode = .false.
  logical :: start_galaxy = .false.
  logical :: lresuming_run = .false.
  character(len=8) :: date
  double precision :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, length, ncycles, flush_signal
  integer, parameter :: finished_tag = 0
  integer, parameter :: newjob_tag = 17
  integer, dimension(MPI_STATUS_SIZE) :: status

  call MPI_INIT(ierr)
  if (ierr/= MPI_SUCCESS) then
    call message('Error starting MPI program. Terminating.')
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  !Get the number of processors this job
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
  !Get the rank of the processor this thread is running on
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr/= MPI_SUCCESS) then
    call message('Error getting processor rank. Terminating.')
    call MPI_Abort(MPI_COMM_WORLD, rc, ierr)
  endif

  if (rank == master_rank) then !Only the master (rank 0)
    tstart = MPI_wtime()
  else
    tstart = -1
    tfinish = -1
  endif

  ! Tries to read the parameter filename from the command argument (or --help)
  call get_command_argument(1, command_argument)
  ! If --help is detected, prints help information
  if (trim(command_argument)=='--help' .or. trim(command_argument)=='-h') then
    call get_command_argument(0, command_argument)
    if (rank==master_rank) then
      print *, 'Magnetizer '
      print *,
      print *, 'Computes ISM properties and solves mean field dynamo equation'&
               //' for the output of a semi-analytic galaxy formation model.'
      print *,
      print *, 'Usage:'
      print *, trim(command_argument), ' <input_parameters_file> [galaxy number]'
      print *,
      print *, 'For more details pleas visit: '&
             //'https://github.com/luizfelippesr/magnetizer'
      print *,
    endif
    stop
  endif

  ! Welcome messages
  if (nproc==1) then
    call message('Magnetizer', rank=-1)
    call message(' ', rank=-1)
    call message('Runnning on a single processor')
  else
    call message('Magnetizer', rank=rank, set_info=info, &
                 master_only=.true., info=0)
    call message(' ', master_only=.true.)
    call message('Runnning on', val_int=nproc, msg_end='processors', &
                 master_only=.true., info=0)
  endif

  if (len_trim(command_argument) == 0) then
    ! Uses example parameter file if nothing was found
    command_argument = 'example/example_global_parameters.in'
    call read_global_parameters(trim(command_argument))
    call message('No parameter file provided. Using standard: '// &
                  trim(command_argument), master_only=.true.,&
                   set_info=info)
  else
    ! Uses specified parameter file
    call read_global_parameters(trim(command_argument))
    call message('Using global parameters file: '// trim(command_argument), &
                 master_only=.true., set_info=info)

    ! Checks whether single galaxy mode was activated
    call get_command_argument(2, command_argument)
    if (len_trim(command_argument) /= 0) then
      call message('Single galaxy mode. igal='// trim(command_argument), &
                 master_only=.true., set_info=info)
      lsingle_galaxy_mode = .true.
    endif
  endif

  ! Checks whether this is a new run or if one is resuming an older run
  ! Note: resuming runs currently only works for separate output files
  if (p_IO_separate_output) then
    if (rank == master_rank) then
      inquire (file=trim(output_file_name), exist=lresuming_run)
    endif
    ! Broadcast the the result
    call MPI_Bcast(lresuming_run, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  endif

  ! Initializes IO (this also reads ngals from the hdf5 input file)
  call IO_start(MPI_COMM_WORLD, MPI_INFO_NULL, lresuming_run)

  ! Allocates a arrays to for logging
  nmygals = 0
  allocate(mygals(ngals))
  mygals = -17 ! Initializes to bad value
  allocate(allgals(ngals*nproc))
  flush_signal = ngals+42 ! Arbitrary larger-than-ngals value

  if (lsingle_galaxy_mode .and. (rank == master_rank)) then
    igal = str2int(command_argument)
    ncycles = 0
    call message('Starting',gal_id=igal, rank=rank)
    ! Check whether the galaxy was processed previously
    start_galaxy = IO_start_galaxy(igal)
    if (start_galaxy) then
      ! If it is a new galaxy, runs it!
      call dynamo_run(igal, p_no_magnetic_fields_test_run, rank, error)
      nmygals = nmygals + 1
      mygals(nmygals) = igal
    endif
  else
    ! Before distributing work between workers, master will try to run a single
    ! galaxy. This helps catching any obvious problem and allows one to find
    ! out where a possible previous run has stopped.
    if (rank==master_rank) then

      do igal=1, ngals
        start_galaxy = IO_start_galaxy(igal)
        if (start_galaxy) then
          call dynamo_run(igal, p_no_magnetic_fields_test_run, rank, error)
          if (error) cycle
          nmygals = nmygals + 1
          mygals(nmygals) = igal
          exit
        endif
      enddo
    endif
    ! Calculates the number of cycles (under the assumption of
    ncycles = max(ngals/p_ncheckpoint, 1)
    call message('Total number of cycles',val_int=ncycles, master_only=.true.)
  endif


  ! ----------
  !   Master
  ! ----------
  if (rank==master_rank) then
    do j=0, ncycles
      if (lsingle_galaxy_mode) exit
      ! Finds boundaries of present cycle
      igal_first = 1+j*p_ncheckpoint
      igal_last = (j+1)*p_ncheckpoint
      call message('Cycle',val_int=j+1, master_only=.true.)
      ! If previous work had been done, accounts for it.
      if (igal>igal_last) cycle
      igal_first = max(igal_first, igal)

      ! Submit initial jobs
      do iproc=1, nproc-1
          igal = iproc-1+igal_first
          if (igal > igal_last) igal=ngals+1
          call MPI_Send(igal, 1, MPI_INTEGER, iproc, newjob_tag, MPI_COMM_WORLD, ierr)
      enddo

      ! Loops sending and receiving
      do igal=nproc+igal_first-1, igal_last
        ! Sends finished work from worker
        if (igal>ngals) exit
        if ( p_master_works_too .and. &
            modulo(igal, (nproc + int(nproc/p_master_skip))) == 0) then
          ! Call dynamo code for each galaxy in mygals
          call message('Starting',gal_id=igal, rank=rank)
          ! Check whether the galaxy was processed previously
          start_galaxy = IO_start_galaxy(igal)
          if (start_galaxy) then
            ! If it is a new galaxy, runs it!
            call dynamo_run(igal, p_no_magnetic_fields_test_run, rank, error)
            nmygals = nmygals + 1
            mygals(nmygals) = igal
          endif
        else
          ! Receives finished work from worker
          call MPI_Recv(igal_finished, 1, MPI_INTEGER, MPI_ANY_SOURCE, finished_tag, &
                        MPI_COMM_WORLD, status, ierr)
          ! Finds out which worker has just finished this work
          iproc = STATUS(MPI_SOURCE)
          ! Sends new job to that worker
          call MPI_Send(igal, 1, MPI_INTEGER, iproc, newjob_tag, MPI_COMM_WORLD, &
                        ierr)
        endif
      enddo

      ! Checks whether a STOP file exists
      inquire (file='STOP', exist=lstop)
      ! If yes, exits gently
      if (lstop) then
        call message('STOP file found. Will finish this cycle and then exit.', &
                     info=0, master_only=.true.)
      endif

      ! Goes through all workers, requesting them to either stop or flush
      ! the data to the disk
      i = 0
      do iproc=1, nproc-1
          ! Receives whatever message it is sending
          call MPI_Recv(igal_finished, 1, MPI_INTEGER, iproc, finished_tag, &
                        MPI_COMM_WORLD, status, ierr)
          if (lstop .or. (igal_last>ngals)) then
            ! If it is time to stop, tells to workers to finish
            call MPI_Send(-1, 1, MPI_INTEGER, iproc, newjob_tag, &
                          MPI_COMM_WORLD, ierr)
          else
            ! Otherwise, sends a signal requesting the worker to flush
            call MPI_Send(flush_signal, 1, MPI_INTEGER, iproc, newjob_tag, &
                          MPI_COMM_WORLD, ierr)
          endif
      enddo

      ! Tells the master to stop
      if (lstop) exit
      ! Tells the master to flush
      call IO_flush()
    enddo
  ! ----------
  !   Worker
  ! ----------
  else
      do
        ! Receives new task or a flag
        call MPI_Recv(igal, 1, MPI_INTEGER, 0, newjob_tag, MPI_COMM_WORLD, status, ierr)
        ! If received an exit flag, exits
        if (igal<0) then
          call message('Exiting', rank=rank)
          exit
        else if (igal<=ngals) then
          ! Call dynamo code for each galaxy in mygals
          call message('Starting',gal_id=igal, rank=rank)
          ! Check whether the galaxy was processed previously
          start_galaxy = IO_start_galaxy(igal)
          if (start_galaxy) then
            ! If it is a new galaxy, runs it!
            call dynamo_run(igal, p_no_magnetic_fields_test_run, rank, error)
            nmygals = nmygals + 1
            mygals(nmygals) = igal
          endif

          ! Sends the result (to mark it done)
          call MPI_Send(igal, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
        else if (igal==flush_signal) then
          call IO_flush()
          call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)

        else
          ! Invalid galaxy. Probably, more processors than galaxies!) Nothing is done.
          call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
        endif
     enddo
  endif


  call IO_flush()
  call message('All computing done', info=0)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (nmygals>0) then
    if (nmygals<10) then
      print *, trim(str(rank)),':  Finished working on galaxies:', mygals(:nmygals),'.'
    else
      call message('Finished working on', val_int=nmygals, msg_end='galaxies.')
    endif
  else
    call message('Finished without working on any galaxy.')
  endif

  call MPI_Gather(mygals, ngals, MPI_INTEGER, allgals, ngals, &
                  MPI_INTEGER, master_rank, MPI_COMM_WORLD, ierr)

  ! Gets the date
  if (rank == master_rank) then
    call date_and_time(date=date)
  end if
  ! Broadcast the date
  call MPI_Bcast(date, 8, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Prints a small report
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! Finalizes IO
  call IO_end(date)
  call message('', master_only=.true.)
  call message('IO finished', master_only=.true.)

  if (rank == master_rank) then !Only the master (rank 0)
    tfinish= MPI_WTime()
    ! Removes stop file if necessary
    if (lstop) then
      i = 17
      open(unit=i, FILE='STOP')
      close(i, status='delete')
      call message('Removed STOP file.', master_only=.true.)
    endif
  endif

  if (info>1 .and. rank==master_rank) then
    i = 1
    do j=1,ngals*nproc
      igal = allgals(j)
      if (igal>0) then
        mygals(i) = igal
        i = i+1
      endif
    enddo
    ngals = i-1
    if (ngals<200 .and. ngals>1) then
      print *,
      print *, '  Galaxies in this run:'
      print *, mygals(:ngals)
      print *,
    else
      print *,
      print *, '  Number of galaxies in this run:', ngals
      print *,
    endif
  endif

  call message('Total wall time in seconds =',tfinish-tstart, &
               master_only=.true., info=0)
  if (.not.lsingle_galaxy_mode) then
    call message('Wall time per galaxy =', (tfinish-tstart)/ngals, &
                 master_only=.true., info=0)
    call message('Average CPU per galaxy =', (tfinish-tstart)*nproc/ngals, &
                 master_only=.true., info=0)
  endif
  !Tell the MPI library to release all resources it is using
  call MPI_Finalize(ierr)
  call message('MPI finished', master_only=.true.)

end program magnetizer
