!*****************************************************
module teste
  implicit none

contains

  subroutine set_random_seed(rank)
    ! Initializes the random seed based on igal and p_random_seed
    ! (actually, any two integers!)
    ! The optional argument use_array uses a more complex seed
    ! (possibly making the random sequence more random) but
    ! a the cost of making the RNG platform dependent.
    integer, intent(in) :: rank
    integer, allocatable, dimension(:) :: seed
    integer :: i, n
    double precision :: p

    ! Discovers the size of the seed array and allocates it
    call random_seed(size = n)
    allocate(seed(n))
    ! Does *arbitrary* calculation to combine p_random_seed parameter and igal
    seed = 17*rank
    ! Changes the seed!
    call random_seed(put=seed)
  end subroutine set_random_seed

  subroutine compute(gals, rank)
    integer, dimension(:), intent(in) :: gals
    integer, intent(in) :: rank
    integer :: i
    real :: a

    do i=1, size(gals)
      call random_number(a)
      print *, rank, ': Working on', gals(i)
      call sleep(int(a*2))
    enddo
  end subroutine
end module teste
program magnetizer
  use mpi
  use teste
  implicit none

  integer, dimension(1) :: ind
  integer :: igal, jgal, nmygals, ngals_start, ngals
  integer, parameter :: master_rank = 0
  integer, allocatable, dimension(:) :: mygals, allgals
  character(len=100) :: command_argument
  integer :: i, j, finished, gal, ncycles, ncheckpoint
  double precision :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, length
  character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
  integer, dimension(MPI_STATUS_SIZE) :: status
  logical :: master_works_too = .true.
  integer :: master_skip = 3
  integer :: i_first_gal, i_last_gal, iproc
  integer, parameter :: finished_tag = 0
  integer, parameter :: newjob_tag = 17


  call MPI_INIT(ierr)
  if (ierr/= MPI_SUCCESS) then
    print *, 'Error starting MPI program. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr) !Get the rank of the processor this thread is running on
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr) !Get the number of processors this job
  call MPI_Get_Processor_Name(hostname, length, ierr) !Get the name of this processor (usually the hostname)
  if (ierr/= MPI_SUCCESS) then
    print *, 'Error getting processor hostname. Terminating.'
    call MPI_Abort(MPI_COMM_WORLD, rc, ierr)
  endif
  call set_random_seed(rank)

  if (rank == master_rank) then !Only the master (rank 0)
    tstart = MPI_wtime()
  else
    tstart = -1
    tfinish = -1
  endif

  ngals = 200

  ! Distributes galaxies between processors
  allocate(mygals(ngals))
  allocate(allgals(ngals*nproc))
  mygals = -17
  allgals= -10000

  ncheckpoint = 18
  ncycles = ngals/ncheckpoint
  nmygals = 0

  ! ----------
  !   Master
  ! ----------
  if (rank==master_rank) then
    do j=0, ncycles
      i_first_gal = 1+j*ncheckpoint
      i_last_gal = (j+1)*ncheckpoint
      print *, '--- cycle',j+1,'---'
      ! Submit initial jobs
      do iproc=1, nproc-1
          igal = iproc-1+i_first_gal
!           print *, igal, iproc
          if (igal > i_last_gal) igal=ngals+1
          print *, '  igal',igal, 'sending to', iproc
          call MPI_Send(igal, 1, MPI_INTEGER, iproc, newjob_tag, MPI_COMM_WORLD, ierr)
      enddo
      print *, '---'
      ! Loops sending and receiving
      do igal=nproc+i_first_gal-1, i_last_gal
        if (&
            (master_works_too .and. modulo(igal,(nproc+nproc/master_skip))==0)) then
          nmygals = nmygals + 1
          call compute([igal], rank)
          mygals(nmygals ) = igal
          allgals(igal) = igal
          print *, '  igal',igal, 'computing'
        else

          ! Receives finished work from worker
          call MPI_Recv(finished, 1, MPI_INTEGER, MPI_ANY_SOURCE, finished_tag, &
                        MPI_COMM_WORLD, status, ierr)
          iproc = STATUS(MPI_SOURCE)

          ! Sends new job!
          call MPI_Send(igal, 1, MPI_INTEGER, iproc, newjob_tag, MPI_COMM_WORLD, &
                        ierr)
          print *, '  igal',igal, 'sending to', iproc

        endif
      enddo

      do iproc=1, nproc-1
          call MPI_Recv(finished, 1, MPI_INTEGER, iproc, finished_tag, &
                        MPI_COMM_WORLD, status, ierr)
          if (finished>0 .and. finished<=ngals) then
            i = i+1
!             print *,'ttt',i, finished
            allgals(i) = finished
          endif

          if (i_last_gal<ngals) then
            ! Sends an invalid galaxy to keep the workers going
            call MPI_Send(ngals+1, 1, MPI_INTEGER, iproc, newjob_tag, &
                          MPI_COMM_WORLD, ierr)
            print *, '  signal',ngals+1, 'sending to', iproc
          else
            ! Tells to workers to finish
            call MPI_Send(-1, 1, MPI_INTEGER, iproc, newjob_tag, &
                          MPI_COMM_WORLD, ierr)
            print *, '  signal',-1, 'sending to', iproc
          endif
      enddo
    enddo
  ! ----------
  !   Worker
  ! ----------
  else
      do
        ! Receives new task or a flag
        call MPI_Recv(gal, 1, MPI_INTEGER, 0, newjob_tag, MPI_COMM_WORLD, status, ierr)
        ! If received an exit flag, exits
        if (gal<0) then
          print *, rank, 'Exiting'
          exit
        else if (gal<=ngals) then
          ! Does the computation for this galaxy
          call compute([gal], rank)
          nmygals = nmygals + 1
          mygals(nmygals) = gal
          ! Sends the result (to mark it done)
          call MPI_Send(gal, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
        else
          ! Invalid galaxy. (More processors than galaxies!) Nothing is done.
          call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
        endif
     enddo
  endif
  print *, rank,'out'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (nmygals>0) then
    print *, 'rank',rank,'Worked on', nmygals, 'gals:',mygals(:nmygals)
  else
    print *, 'rank',rank,'Worked on no gals'
  endif

  call MPI_Gather(mygals, ngals, MPI_INTEGER, allgals, ngals, &
                  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


  call MPI_Finalize(ierr)
  if (rank==master_rank) then
    print *, 'FIM'
    do i=1,ngals*2
      ind = maxloc(allgals)
      if (allgals(ind(1))>0 .and. allgals(ind(1))<=ngals) print *, allgals(ind(1))
      allgals(ind(1)) = -1000
    enddo
  endif
end program magnetizer
