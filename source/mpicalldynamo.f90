!*****************************************************
program magnetizer
  use mpi
  use dynamo
  use input_params
  use global_input_parameters
  use IO
  use messages

  implicit none

  integer :: igal, jgal, nmygals
  integer, parameter :: master_rank = 0
  integer, allocatable, dimension(:) :: mygals
  character(len=100) :: command_argument
  integer :: i
  logical :: lstop
  logical :: lsingle_galaxy_mode = .false.

  character(len=8) :: date
  double precision :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, len
  character(len=MPI_MAX_PROCESSOR_NAME) hostname

  call MPI_INIT(ierr)
  if (ierr/= MPI_SUCCESS) then
    call message('Error starting MPI program. Terminating.')
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr) !Get the rank of the processor this thread is running on
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr) !Get the number of processors this job
  call MPI_Get_Processor_Name(hostname, len, ierr) !Get the name of this processor (usually the hostname)
  if (ierr/= MPI_SUCCESS) then
    call message('Error getting processor hostname. Terminating.')
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
  if (trim(command_argument)=='--help' .or. &
                  trim(command_argument)=='-h') then
    call get_command_argument(0, command_argument)
    print *, 'Magnetizer '
    print *,
    print *, 'Usage:'
    print *, trim(command_argument), ' <input_parameters_file> [galaxy number]'
    print *,
    print *, 'One day this will be a nice help message'
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
    call message(' ', rank=-1)
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

  ! Initializes IO
  call IO_start(MPI_COMM_WORLD, MPI_INFO_NULL)

  if (.not.lsingle_galaxy_mode) then
    ! Distributes galaxies between processors
    allocate(mygals(ngals))
    nmygals = 0
    do i=0,ngals
      igal = rank + i*nproc+1
      if (igal>ngals) exit
      mygals(i+1) = igal
      nmygals = nmygals + 1
    end do
  else
    allocate(mygals(1))
    if (rank == master_rank) then
      mygals(1) = str2int(command_argument)+1
      nmygals = 1
    else
      nmygals = 0
    endif
  endif

  if (ngals<1000) then
    call message('List of galaxies to run', rank=rank)
    ! The following doesn't work because an interface for printing vectors
    ! wasn't yet writen.
    ! call message('    mygals =',mygals(:nmygals), rank=rank)
    ! Will do it by hand. Later, if this proves useful, the messages module can
    ! be updated.
    print *, str(rank),':  ','mygals = ',mygals(:nmygals)
  endif

  if (nmygals > 0) then
    do jgal=1,nmygals
      ! Call dynamo code for each galaxy in mygals
      call message('Starting',gal_id=mygals(jgal), rank=rank)
      call dynamo_run(mygals(jgal), p_no_magnetic_fields_test_run, rank)
      ! Checks whether a STOP file exists
      inquire (file='STOP', exist=lstop)
      ! If yes, exits gently
      if (lstop) then
        call message('Found STOP file. Exiting..', info=0)
        exit
      endif
    enddo
  endif

  call message('All computing done', info=0)

  ! Gets the date
  if (rank == master_rank) then
    call date_and_time(date=date)
  end if
  ! Broadcast the date
  call MPI_Bcast(date, 8, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Finalizes IO
  call IO_end(date)
  call message('IO finished')

  !Tell the MPI library to release all resources it is using
  call MPI_Finalize(ierr)
  call message('MPI finished')

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

  call message('Total wall time in seconds =',tfinish-tstart, &
               master_only=.true., info=0)
  call message('Average (total) time per galaxy =', &
               (tfinish-tstart)*nproc/ngals, &
               master_only=.true., info=0)
end program magnetizer
