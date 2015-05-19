!*****************************************************
program mpicalldynamo
  use mpi
  use dynamo
  use input_params
!
  implicit none 
!
  integer, parameter :: ngals=8, info=0
  logical, parameter :: master_participate= .false. 
  integer :: igal, jgal, nmygals, nslaves, order, flag
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, parameter :: master=0
  integer, allocatable, dimension(:) :: mygals
!*****************************************************
  real :: tstart,tfinish
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
!  print "('calldynamo.f: Number of tasks=',I3,' My rank=',I3,' My name=',A,'')",nproc, rank, trim(name)
!*****************************************************
!
! Assign a chunk of galaxies to each processor
!
if (rank== 0) then !Only the master (rank 0)
  tstart= MPI_wtime()
endif
!
if (master_participate) then
  nslaves=nproc
  order=rank
else
  nslaves=nproc-1
  order=rank-1
endif
!
nmygals=0
do igal=1,ngals
  if (mod(igal,nslaves) == order) then
    nmygals=nmygals+1
  endif
enddo
allocate(mygals(nmygals))
jgal=1
do igal=1,ngals
  if (mod(igal,nslaves) == order) then
    mygals(jgal)=igal
    jgal=jgal+1
  endif
enddo

print*,'rank=',rank,'    mygals=',mygals !processor id, list of galaxies for that processor

!*****************************************************
!
! Call dynamo code for each galaxy in mygals
if (nmygals>0) then
  do igal=1,nmygals
    call dynamo_run(info, igal, flag)
    print*,'flag',flag,'obtained by processor',rank,'for galaxy',mygals(igal)
  enddo
endif
!
!*****************************************************
!
! intent(in)
!  info:   how much info to print during run
!          options are 0=none, 1=min, 2=standard, 3=max 
!  gal_id: id number of galaxy
!          ranges from 1 to ~300,000
!
! intent(out)
!  flag:   status of run
!          results are -1=run successful, 1=, 2=, 3=
!
! Note: physical variables (ts_t, ts_Br, ts_Bp, ts_Bzmod, ts_alp_m) are written to files, not passed
!
! Format:  call dynamo_run(info, gal_id, flag)
!  Standard parameters: (2, igal, flag), where igal is looped over within python
!
!
  call MPI_FINALIZE(ierr) !Tell the MPI library to release all resources it is using
!
if (rank== 0) then !Only the master (rank 0)
  tfinish= MPI_wtime()
  print*,'total wall time in seconds=',tfinish-tstart
endif
!
!
end program mpicalldynamo
!*****************************************************
