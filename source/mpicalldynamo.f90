!*****************************************************
program mpicalldynamo
  use mpi
  use dynamo
  use input_params
  use global_input_parameters
  use IO
  use messages
!
  implicit none 
!
  logical, parameter :: master_participate= .false. 
  integer :: igal, jgal, nmygals
  integer, parameter :: master=0
  integer, allocatable, dimension(:) :: mygals
  character(len=100) :: command_argument
  integer :: i
  logical :: lstop

  character(len=8) :: date
  double precision :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, len
  character(len=MPI_MAX_PROCESSOR_NAME) hostname

  call MPI_INIT(ierr)
  if (ierr/= MPI_SUCCESS) then
    call message('Error starting MPI program. Terminating.')
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) !Get the rank of the processor this thread is running on
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr) !Get the number of processors this job
  call MPI_GET_PROCESSOR_NAME(hostname, len, ierr) !Get the name of this processor (usually the hostname)
  if (ierr/= MPI_SUCCESS) then
    call message('Error getting processor hostname. Terminating.')
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

  if (rank== 0) then !Only the master (rank 0)
    tstart= MPI_wtime()
  endif

  if (nproc==1) then
    call message('Starting galform magnetizer', rank=-1)
    call message('Runnning on a single processor')
  else
    call message('Starting galform magnetizer', rank=rank, set_info=info, &
                 master_only=.true., info=0)
    call message('Runnning on', val_int=nproc, msg_end='processors', &
                 master_only=.true.)
  endif

  call get_command_argument(1, command_argument)
  if (len_trim(command_argument) == 0) then
    command_argument = 'example/example_global_parameters.in'
    call read_global_parameters(trim(command_argument))
    call message('No parameter file provided. Using standard: '// &
                  trim(command_argument), master_only=.true.,&
                   set_info=info)
  else
    call read_global_parameters(trim(command_argument))
    call message('Using global parameters file: '// trim(command_argument), &
                 master_only=.true., set_info=info)
  endif
  
  ! Initializes IO
  call IO_start(MPI_COMM_WORLD, MPI_INFO_NULL)

  allocate(mygals(ngals))
  nmygals = 0
  do i=0,ngals
    igal = rank + i*nproc+1
    if (igal>ngals) exit
    mygals(i+1) = igal
    nmygals = nmygals + 1
  end do

  call message('List of galaxies to run', rank=rank)
  ! The following doesn't work because an interface for printing vectors
  ! wasn't yet writen.
  ! call message('    mygals =',mygals(:nmygals), rank=rank)
  ! Will do it by hand. Later, if this proves useful, the messages module can
  ! be updated.
  print *, str(rank),':  ','mygals = ',mygals(:nmygals)

  if (nmygals > 0) then
    do jgal=1,nmygals
      ! Call dynamo code for each galaxy in mygals
      call message('Starting',gal_id=mygals(jgal))
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
  if (rank==0) then
    call date_and_time(date=date)
  end if
  ! Broadcast the date
  call MPI_BCAST(date, 8, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Finalizes IO
  call IO_end(date)
  call message('IO finished')

  !Tell the MPI library to release all resources it is using
  call MPI_FINALIZE(ierr)
  call message('MPI finished')

  if (rank == 0) then !Only the master (rank 0)
    tfinish= MPI_wtime()
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
end program mpicalldynamo
!*****************************************************
