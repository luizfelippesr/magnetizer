program test_IO_read
  use grid
  use IO
  use mpi
  use math_constants
  use units
  use messages
  use global_input_parameters
  implicit none
  character(len=100) :: command_argument
  character(len=100) :: date
  logical :: success, incomplete
  integer, parameter :: master_rank = 0
  integer :: rank, nproc, ierr, rc, ncycles, flush_signal
  integer :: igal, info_mpi, i, j, iz
  integer,dimension(8) :: time_vals
  double precision :: alpha, b
  double precision, allocatable, dimension(:,:) :: Br, Bp, Bzmod, Rcyl,xc2
  double precision, allocatable, dimension(:,:) :: Bpara, Bperp
  double precision, allocatable, dimension(:) :: Bpara_valid, Bperp_valid
  double precision, allocatable, dimension(:) :: Bx, By, Bz, Bmag, xc, zc
  double precision, allocatable, dimension(:) :: angle_B_LoS, tmp
  logical, allocatable, dimension(:,:) :: valid
  double precision, allocatable, dimension(:,:) :: Rpath
  integer, allocatable, dimension(:) :: js
  double precision, allocatable, dimension(:) :: test

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

  ! For testing
  igal = 2
  b = 2 ! kpc, Impact parameter
  alpha = pi/2d0  ! kpc, angle relative to the normal to the midplane


  print *, 'Starting Galaxy', igal
  incomplete = IO_start_galaxy(igal)

  allocate(Rcyl(number_of_redshifts,p_nx_ref))
  call IO_read_dataset_vector('r', igal, Rcyl, group='Output')

  allocate(Br(number_of_redshifts,p_nx_ref))
  call IO_read_dataset_vector('Br', igal, Br, group='Output')

  allocate(Bp(number_of_redshifts,p_nx_ref))
  call IO_read_dataset_vector('Bp', igal, Bp, group='Output')

  !lfsr Later, we will need to account for the sign of Bz...
  allocate(Bzmod(number_of_redshifts,p_nx_ref))
  call IO_read_dataset_vector('Bzmod', igal, Bzmod, group='Output')

  allocate(xc(number_of_redshifts))
  allocate(xc2(number_of_redshifts, 2*p_nx_ref))
  allocate(zc(number_of_redshifts))
  allocate(js(2*p_nx_ref))
  allocate(valid(number_of_redshifts,2*p_nx_ref))
  valid = .false.
  allocate(Bpara(number_of_redshifts,2*p_nx_ref))
  allocate(Bperp(number_of_redshifts,2*p_nx_ref))
  allocate(Bx(number_of_redshifts))
  allocate(By(number_of_redshifts))
  allocate(Bmag(number_of_redshifts))
  allocate(angle_B_LoS(number_of_redshifts))
!   allocate(Bz(number_of_redshifts))

  ! Construct coordinates and auxiliary indices js
  do i=1,2*p_nx_ref
    if (i<p_nx_ref+1) then
      ! Index for behind the x-z plane
      j=p_nx_ref-i+1
      js(i) = j

      ! Coordinates
      xc = -1
    else
      ! Index ahead of the x-z plane
      j = i-p_nx_ref
      js(i) = j

      ! Coordinates
      xc = 1
    end if

    ! The available values for x are x_i = sqrt(R_i^2-b^2)
    tmp = Rcyl(:,j)**2-b**2
    where (tmp>0)
      xc = xc * sqrt(tmp)
      ! Stores a mask to be used with the fortran 2003 'pack' function
      valid(:,i) = .true.
    endwhere

    ! Bx = Br * x/R - Bp * y/R
    Bx = Br(:,j) * xc/Rcyl(:,j) - Bp(:,j)* b/Rcyl(:,j)
    ! By = Br * y/R + Bp * x/R
    By = Br(:,j) * b/Rcyl(:,j) - Bp(:,j)* xc/Rcyl(:,j)
    ! Bz = Bzmod * ?
    Bz = Bzmod(:,j)

    ! Simple vector calculations
    ! |B|
    Bmag = sqrt(Bx**2 + By**2 + Bz**2)
    ! B_\parallel = dot(B,n)
    Bpara(:,i) = Bx*cos(alpha) + Bz*sin(alpha)
    ! angle = arccos(Bpara/|B|)
    angle_B_LoS = acos(Bpara(:,i)/Bmag)
    ! B_\perp = |B|*sin(angle) -- magnitude of the perpendicular component
    Bperp(:,i) = Bmag*sin(angle_B_LoS)
    xc2(:,i) = xc
  enddo

  ! Now work is done for each redshift (as the valid section of each array may
  ! may be different)
  do iz=1,number_of_redshifts
    ! Filters away invalid part of an array
    Bpara_valid = pack(Bpara(iz,:),valid(iz,:))
    print *, shape(Bpara_valid)
  enddo

  print *,'----'
  test = pack(xc2(2,:),valid(2,:))
  print *, shape(test)
  print *, test
  do i=1,size(test)
    print *, test(i)
  enddo


end program test_IO_read
