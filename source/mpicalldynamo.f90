!*****************************************************
program mpicalldynamo
  use mpi
  use dynamo
  use input_params
  use global_input_parameters
  use IO
!
  implicit none 
!
  logical, parameter :: master_participate= .false. 
  integer :: igal, jgal, nmygals, flag
  integer, parameter :: master=0
  integer, allocatable, dimension(:) :: mygals
  character(len=32) :: command_argument
  integer :: i

!*****************************************************
  double precision :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, len
  character*(MPI_MAX_PROCESSOR_NAME) name
!
  call MPI_INIT(ierr)
  if (ierr/= MPI_SUCCESS) then
    print*,'Error starting MPI program. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) !Get the rank of the processor this thread is running on
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr) !Get the number of processors this job
  call MPI_GET_PROCESSOR_NAME(name, len, ierr) !Get the name of this processor (usually the hostname)
  if (ierr/= MPI_SUCCESS) then
    print*,'Error getting processor name. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

  if (rank== 0) then !Only the master (rank 0)
    tstart= MPI_wtime()
  endif

  call get_command_argument(1, command_argument)
  if (len_trim(command_argument) == 0) then
    write(*,*) 'No parameter file provided. Using standard global parameters.'
  else
    call read_global_parameters(trim(command_argument))
  endif

  allocate(mygals(ngals))

  nmygals = 0
  do i=0,ngals
    igal = rank + i*nproc+1
    if (igal>ngals) exit
    mygals(i+1) = igal
    nmygals = nmygals + 1
  end do
  
  print*,'rank=',rank,'    mygals=',mygals(:nmygals) !processor id, list of galaxies for that processor
  ! Initializes IO  
  call IO_start(MPI_COMM_WORLD, MPI_INFO_NULL)

  ! Call dynamo code for each galaxy in mygals
  if (nmygals > 0) then
    do jgal=1,nmygals
      call dynamo_run(info, mygals(jgal), flag, &
                      p_no_magnetic_fields_test_run)
      print*,'flag',flag,'obtained by processor',rank,'for galaxy',mygals(jgal)
    enddo
  endif
  
  print*,'rank=',rank,'    All done'
  ! Finalizes IO
  call IO_end(info)
  print*,'rank=',rank,'    IO finished'

  !Tell the MPI library to release all resources it is using
  call MPI_FINALIZE(ierr)
  
  print*,'rank=',rank,'    MPI finished'
  
  if (rank == 0) then !Only the master (rank 0)
    tfinish= MPI_wtime()
    print*,'total wall time in seconds=',tfinish-tstart
  endif
end program mpicalldynamo
!*****************************************************
