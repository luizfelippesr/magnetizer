

program test_IO_read
  use grid
  use IO
  use mpi
  use messages
  use global_input_parameters
  implicit none
  character(len=100) :: command_argument
  character(len=100) :: date
  logical :: success, incomplete
  integer, parameter :: master_rank = 0
  integer :: rank, nproc, ierr, rc, ncycles, flush_signal
  integer :: igal, info_mpi
  integer,dimension(8) :: time_vals

  ! Initializes MPI
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
      print *, trim(command_argument), ' <input_parameters_file> [galaxy number] [-f]'
      print *,
      print *, 'For more details please visit: '&
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
  endif


  if (rank==master_rank) then
    call date_and_time(values=time_vals)
    ! Displays date in the format yyyy-mm-dd  hh-mm-ss timezone
    date = str(time_vals(1))//'-'//str(time_vals(2))//'-'//str(time_vals(3)) &
                          //'  '//str(time_vals(5))//':'//str(time_vals(6))&
                          //' (UTC'//str(time_vals(4))//')'
  endif
  ! It is probably a good idea to sync, to avoid different nodes having
  ! slightly different times, for any reason (the attribute creation is
  ! collective, anyway).
  call MPI_Bcast(date, len_trim(date), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Avoids MPI errors
  call MPI_Info_create(info_mpi, IERR)
  call MPI_Info_set(info_mpi, "romio_ds_write", "disable", ierr)
  call MPI_Info_set(info_mpi, "romio_ds_read", "disable", ierr)

  ! Initializes IO (this also reads ngals from the hdf5 input file)
  call IO_start(MPI_COMM_WORLD, info_mpi, .true., date)

  igal = 2
  print *, 'Starting Galaxy', igal
  incomplete = IO_start_galaxy(igal)

  IO_read_dataset_scalar

end program test_IO_read