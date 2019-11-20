program Observables_single
  use grid
  use IO
  use mpi
  use math_constants
  use units
  use messages
  use global_input_parameters
  use LoSintegrate_aux
  implicit none
  character(len=300) :: command_argument, image_dir
  character(len=50) :: run_type
  character(len=100) :: date
  logical :: incomplete
  integer, parameter :: master_rank = 0
  integer :: rank, nproc, ierr, rc
  integer :: info_mpi
  integer,dimension(8) :: time_vals
  double precision :: error
  type(Galaxy_Properties) :: props
  type(LoS_data) :: data
  double precision, allocatable, dimension(:,:) :: buffer
  double precision :: res, impact_y, impact_z
  integer :: it, i

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

  ! Tries to read the parameter filename from the command argument
  call get_command_argument(2, command_argument)
  if (len_trim(command_argument) == 0) then
    ! Exits if nothing was found
    call get_command_argument(0, command_argument)

    call message('Usage: '//trim(command_argument)//' <type> <parameters_file> <gal_id> <iz>'// &
                  ' <theta> <ignore_small_scale_field> [y] [z] [image_dir]',&
                  master_only=.true.)
    call message('', master_only=.true.)
    call message('type ', master_only=.true.)
    call message('      "RM": compute backlit RM along LoS; ', master_only=.true.)
    call message('      "I": total integrated synchrotron emission;',&
                    master_only=.true.)
    call message('      "Image": outputs images of Stokes I,U,Q and RM to image_dir',&
                    master_only=.true.)
    call message('      "RM_study": computes RM for various LoS with the same theta.',&
                    master_only=.true.)
    call message('parameters_file - Magnetizer parameters file used for the run', &
                 master_only=.true.)
    call message('gal_id - galaxy index', master_only=.true.)
    call message('iz - redshift index', master_only=.true.)
    call message('theta - angle of the line of sight (radians)', master_only=.true.)
    call message('wavelength - in meters', master_only=.true.)
    call message('ignore_small_scale_field - If 1, only the mean field will be used', &
                 master_only=.true.)
    call message('y - if "RM", position at which the LoS intercepts the x-z ' &
                //'plane (in units of rmax), if "Image", maximum y-intercept (kpc)',&
                 master_only=.true.)

    call message('z - if "RM", position at which the LoS intercepts the x-y ' &
                //'plane (in units of rmax), if "Image", maximum z-intercept (kpc)',&
                 master_only=.true.)
    stop
  else
    ! Uses specified parameter file
    call read_global_parameters(trim(command_argument))
    call message('Using global parameters file: '// trim(command_argument), &
                 master_only=.true., set_info=info)
  endif

  call get_command_argument(1, command_argument)
  run_type = trim(command_argument)

  ! Welcome messages
  if (nproc==1) then
!     call message('Magnetizer', rank=-1)
    call message(' ', rank=-1)
!     call message('Runnning on a single processor')
  else
    stop
!     call message('Magnetizer', rank=rank, set_info=info, &
!                  master_only=.true., info=0)
!     call message(' ', master_only=.true.)
!     call message('Runnning on', val_int=nproc, msg_end='processors', &
!                  master_only=.true., info=0)
  endif

!     call message('', master_only=.true., set_info=info)


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

  print *, '   Mode = ', run_type
  ! Reads the other command arguments
  ! <gal_id>
  call get_command_argument(3, command_argument)
  props%igal = str2int(command_argument)
  print *, '   igal = ', props%igal
  ! <iz> - Index of the selected redshift
  call get_command_argument(4, command_argument)
  it = str2int(command_argument)
  print *, '   iz = ', it
  ! <theta> - Relative orientation of the galaxy/LoS
  ! The line-of-sight (LOS) is assumed to be parallel to y-z plane
  ! Angle relative to the normal to the midplane
  call get_command_argument(5, command_argument)
  data%theta = str2dbl(command_argument)
  print *, '   theta = ', data%theta

  call get_command_argument(6, command_argument)
  data%wavelength = str2dbl(command_argument)
  print *, '   wavelenth = ', data%wavelength

  ! Sets other parameters (currently, hard-coded)
  data%alpha = 3d0 ! Spectral index of the cr energy distribution
  ! Only one RM to be computed in this mode (Observables_single)
  props%n_RMs = 1
  data%B_scale_with_z = .false.
  call get_command_argument(7, command_argument)
  if (str2int(command_argument)==1) then
    data%ignore_small_scale_field = .true.
  else
    data%ignore_small_scale_field = .false.
  endif
  print *, '   ignore_small_scale_field = ', data%ignore_small_scale_field

  ! Initializes Magnetizer IO and reads relevant galaxy properties
  incomplete = IO_start_galaxy(props%igal)
  if (incomplete) then
    call error_message('','Galaxy not complete', abort=.true.)
  endif
  ! Prepares a derived data type carrying the galaxy properties
  call alloc_Galaxy_Properties(number_of_redshifts,p_nx_ref, props)
  ! The reading below requires the use of a buffer variable, possibly
  ! due to something in the HDF5 library
  allocate(buffer(number_of_redshifts,p_nx_ref))
  call IO_read_dataset_vector('r', props%igal, buffer, group='Output')
  props%Rcyl = buffer
  call IO_read_dataset_vector('Br', props%igal, buffer, group='Output')
  props%Br = buffer
  call IO_read_dataset_vector('Bp', props%igal, buffer, group='Output')
  props%Bp = buffer
  call IO_read_dataset_vector('Bzmod', props%igal, buffer, group='Output')
  props%Bz = buffer
  call IO_read_dataset_vector('h', props%igal, buffer, group='Output')
  props%h = buffer/1d3 ! Converts from pc to kpc
  call IO_read_dataset_vector('n', props%igal, buffer, group='Output')
  props%n = buffer


  ! Prepares image if y and zarguments are present
  call get_command_argument(8, command_argument)
  if (len_trim(command_argument) /= 0) then
    ! [ymax]
    impact_y = str2dbl(command_argument)
  else
    impact_y = -1000d0
  endif

  call get_command_argument(9, command_argument)
  if (len_trim(command_argument) /= 0) then
    ! [zmax]
    impact_z = str2dbl(command_argument)
  else
    impact_z = -1000d0
  endif

   ! [image_dir]
   call get_command_argument(10, command_argument)
   if (len_trim(command_argument) /= 0) then
     i = len_trim(command_argument)
     if (command_argument(i:i)/='/') then
        command_argument = trim(command_argument)//'/'
     endif
     image_dir = command_argument
   endif

  select case (trim(run_type))
    case ('Image')
      call message('Preparing images for Q, U, I and RM', gal_id=props%igal, &
                   msg_end='saving to dir: '//image_dir)
      ! [image_dir]
      call print_image(props, data, image_dir, impact_y, impact_z, &
                        nprint=70, isnap=it)
    case ('I')
      call message('Computed integrated synchrotron intensity for this choie of theta', gal_id=props%igal)
      res = IntegrateImage('I', props, data, it, 1700,'VEGAS', error)
      print *, 'Total synchrotron intensity:'
      print *, res

    case ('PI')
      call message('Computed integrated synchrotron intensity for this choie of theta', gal_id=props%igal)
      res = IntegrateImage('PI', props, data, it, 1700,'VEGAS', error)
      print *, 'Total synchrotron polarised intensity:'
      print *, res

    case ('RM')
      call message('Computing RM for this choice of y, z and theta', gal_id=props%igal)
      call LoSintegrate(props, impact_y, impact_z, data, it, RM_out=.true., &
                        I_out=.false., Q_out=.false., U_out=.false.)
      print *, 'RM:'
      print *, data%RM

    case ('RM_study')
      call message('Computing RM for this of theta', gal_id=props%igal)
      do i=1,200000
        call random_number(impact_y)
        impact_y = (impact_y*2d0-1d0)!/2.
        call random_number(impact_z)
        impact_z = (impact_z*2d0-1d0)!/2.
        res = sqrt(impact_z**2 + impact_y**2)
        impact_z = impact_z/sin(data%theta)
        call LoSintegrate(props, impact_y, impact_z, data, it, RM_out=.true., &
                          I_out=.false., Q_out=.false., U_out=.false.)
        print *, data%RM(it), res, impact_y, impact_z
      enddo
    case ('RM_study_r')
      call message('Computing RM for this choice of y, z and theta', gal_id=props%igal)
      do i=1,100000
        call random_number(res)
        call random_number(impact_y)
        impact_y = (impact_y*2d0-1d0)*res
        call random_number(impact_z)
        impact_z = (impact_z*2d0-1d0)
        impact_z = sign(sqrt(res**2-impact_y**2), impact_z)
        impact_z = impact_z/sin(data%theta)
        call LoSintegrate(props, impact_y, impact_z, data, it, RM_out=.true., &
                          I_out=.false., Q_out=.false., U_out=.false.)
        print *, data%RM(it), res, impact_y, impact_z
      enddo
    case default
      stop '?'
  end select
end program Observables_single
