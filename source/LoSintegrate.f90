program test_IO_read
  use grid
  use IO
  use mpi
  use math_constants
  use units
  use messages
  use global_input_parameters
  use LoSintegrate_aux
  implicit none
  character(len=100) :: command_argument
  character(len=100) :: date
  logical :: success, incomplete
  integer, parameter :: master_rank = 0
  integer, parameter :: nz = 100
  integer :: rank, nproc, ierr, rc, ncycles, flush_signal
  integer :: igal, info_mpi, i, j, iz, it
  integer,dimension(8) :: time_vals
  double precision :: alpha, b
  double precision, allocatable, dimension(:,:) :: Br, Bp, Bzmod, Rcyl, h
  double precision, allocatable, dimension(:,:) :: Bpara_all, Bperp_all, xc, zc
  double precision, allocatable, dimension(:) :: Bpara, Bperp
  double precision, allocatable, dimension(:) :: Bx, By, Bz, Bmag
  double precision, allocatable, dimension(:) :: angle_B_LoS, tmp
  logical, allocatable, dimension(:,:) :: valid
  double precision, allocatable, dimension(:,:) :: Rpath
  integer, allocatable, dimension(:) :: js
  double precision, allocatable, dimension(:) :: test
  logical, parameter :: l_B_scale_with_z = .false.
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
  alpha = 1.51843644923506673192  ! kpc, angle relative to the normal to the midplane


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

  allocate(h(number_of_redshifts,p_nx_ref))
  call IO_read_dataset_vector('h', igal, h, group='Output')
  h = h * 1d-3 ! Converts from pc to kpc


  allocate(xc(number_of_redshifts, 2*p_nx_ref))
  allocate(zc(number_of_redshifts, 2*p_nx_ref))
  allocate(js(2*p_nx_ref))

  allocate(valid(number_of_redshifts,2*p_nx_ref))
  valid = .false.
  allocate(Bpara_all(number_of_redshifts,2*p_nx_ref))
  allocate(Bperp_all(number_of_redshifts,2*p_nx_ref))
  allocate(Bmag(number_of_redshifts))
  allocate(angle_B_LoS(number_of_redshifts))
  allocate(tmp(number_of_redshifts))

  ! i -> index for quantities in a cartesian box
  ! j -> index for axi-symmetric quantities
  do i=1,2*p_nx_ref
    ! Constructs coordinates and auxiliary indices js
    if (i<p_nx_ref+1) then
      ! Index for behind the x-z plane
      j=p_nx_ref-i+1
      js(i) = j
      ! sign for x coordinate
      xc(:,i) = -1
    else
      ! Index ahead of the x-z plane
      j = i-p_nx_ref
      js(i) = j
      ! sign for x coordinate
      xc(:,i) = 1
    end if

    ! The available values for x are x_i = sqrt(R_i^2-b^2)
    tmp = Rcyl(:,j)**2-b**2 ! NB allocated on-the-fly: Fortran2003 feature
    where (tmp>0)
      xc(:,i) = xc(:,i) * sqrt(tmp)
      ! Stores a mask to be used with the Fortran 2003 'pack' intrinsic function
      valid(:,i) = .true.
    endwhere

    ! Sets z-coord
    zc(:,i) = xc(:,i) / tan(alpha)

    ! Bx = Br * x/R - Bp * y/R  (NB allocated on-the-fly)
    Bx = Br(:,j) * xc(:,i)/Rcyl(:,j) - Bp(:,j)* b/Rcyl(:,j)
    ! By = Br * y/R + Bp * x/R (NB allocated on-the-fly)
    By = Br(:,j) * b/Rcyl(:,j) - Bp(:,j)* xc(:,i)/Rcyl(:,j)
    ! Bz = Bzmod * ? (NB allocated on-the-fly)
    Bz = Bzmod(:,j)

    ! Simple vector calculations
    ! |B|
    Bmag = sqrt(Bx**2 + By**2 + Bz**2)
    ! B_\parallel = dot(B,n), where n is the LoS direction
    Bpara_all(:,i) = Bx*cos(alpha) + Bz*sin(alpha)
    ! angle = arccos(Bpara/|B|)
    angle_B_LoS = acos(Bpara_all(:,i)/Bmag)
    ! B_\perp = |B|*sin(angle) -- magnitude of the perpendicular component
    Bperp_all(:,i) = Bmag*sin(angle_B_LoS)

    tmp = abs(zc(:,i))/h(:,j)
    if (l_B_scale_with_z) then
      ! Scale with z (coordinate)
      Bperp_all(:,i) = Bperp_all(:,i) * exp(-tmp)
      Bpara_all(:,i) = Bpara_all(:,i) * exp(-tmp)
    else
      ! Constant for |z|<h, zero otherwise
      where (tmp>1)
        Bperp_all(:,i) = 0d0
        Bpara_all(:,i) = 0d0
      endwhere
    endif
    j = 2
!     print *, xc(j,i), zc(j,i), Bperp_all(j,i), Bpara_all(j,i), tmp(j), valid(j,i)
  enddo
! stop
  ! Now work is done for each redshift (as the valid section of each array may
  ! may be different)
  do it=1,number_of_redshifts


    ! Filters away invalid part of an array
    Bpara = pack(Bpara_all(it,:),valid(it,:))
    ! Includes z dependence

    test = pack(xc(it,:),valid(it,:))
    print *, shape(test)
    print *,'----'
    print *, '   xc        '
    do i=1,size(test)
      print *, test(i), Bpara(i)
    enddo
      print *,'----'
  enddo



end program test_IO_read
