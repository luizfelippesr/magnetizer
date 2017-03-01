! Contains two modules which deal with input parameters
module input_params
  ! Reads and sets input parameters
  use global_input_parameters
  use math_constants

  implicit none

  character (len=200) :: header_time_indep, header_time_dep
  integer :: iread=0

  !
  logical :: last_output = .false.

  ! Time-stepping parameters
  integer, protected :: n1 = -1 !Total number of snapshots (set to -1 to flag it is uninitialized)
  integer, protected :: init_it = -1 ! Index of the initial snapshot (set to -1 to flag it is uninitialized)
  integer, protected :: max_it = -1 ! Index of the final snapshot (set to -1 to flag it is uninitialized)
  double precision, protected :: tsnap !Time between successive snapshots
  integer, protected :: nsteps=20  !Number of timesteps in between snapshots
  double precision, protected :: dt ! Timestep
  double precision :: t=0,first=0.d0 !Timestep variables (lfsr: this looks a bit dangerous..)

  double precision, private :: time_between_inputs=0
  
  ! Galaxy parameters
  double precision, protected :: t_Gyr
  double precision, protected :: r_max_kpc_history
  double precision, protected :: R_kappa
  double precision, protected :: l_sol_kpc, r_l_kpc
  double precision, protected :: v_sol_kms, r_v_kpc
  double precision, protected :: Uz_sol_kms, r_Uz_kpc
  double precision, protected :: r_disk, v_disk
  double precision, protected :: r_bulge, v_bulge
  double precision, protected :: r_halo, v_halo, nfw_cs1
  double precision, protected :: Mstars_disk, Mgas_disk, SFR
  double precision, protected :: Mhalo, Mstars_bulge
  ! other
  double precision :: lambda

  ! Maximum number of columns in the galaxy input files
  integer, private, parameter :: number_of_columns=11
  ! All galaxy data
  double precision, allocatable, dimension(:,:),private :: galaxy_data
  double precision, allocatable, dimension(:), private :: output_times
  integer, private ::  current_gal_id

  contains

  function set_timestep(h,v,etat,reduce_timestep) result(success)
    ! Sets the timestep
    ! if reduce_timestep=False, nsteps=nsteps_0
    ! otherwise nsteps*=2
    use units
    use grid
    use messages, only: message
    logical, intent(in), optional :: reduce_timestep
    double precision, dimension(:), optional :: h,v,etat
    logical :: reduce_ts, success
    double precision :: minh
    integer :: i

    success = .true.

    if (present(reduce_timestep)) then
      reduce_ts = reduce_timestep
    else
      reduce_ts = .false.
    endif

    tsnap = time_between_inputs/t0_Gyr

    if (p_variable_timesteps) then
      if (present(h) .and. present(v) .and. present(etat)) then
        minh = 1d10
        do i=1, size(h)
          if (h(i)>0 .and. h(i) < minh)  minh=h(i)
        enddo

        if (abs(maxval(v)) > 1d-20 .and. abs(maxval(etat)) > 1d-20) then
          dt = minval([ p_courant_v * dx/lambda/maxval(v), &
                        p_courant_eta * minh**2/maxval(etat) ])
          nsteps = int(tsnap/dt)
        else
          call message('set_timestep: max(v) or max(etat) is negligible', &
                       info=1, gal_id=current_gal_id)
          nsteps = p_nsteps_max
          success = .false.
        endif
      else
        call message('set_timestep: missing arguments. Falling back to fixed number of timesteps scheme.', &
                     info=3, gal_id=current_gal_id)
        nsteps = nsteps_0
      endif
    else
      if (.not.reduce_ts) then
        ! Initializes the number of steps to the global input value
        nsteps = nsteps_0
      else
        nsteps = 2*nsteps
      endif
      dt = tsnap/nsteps  !Timestep in units of t0=h0^2/etat0
    endif

    call message('set_timestep: dt = ',dt*t0_Gyr, msg_end='Gyr', info=3, &
                  gal_id=current_gal_id)
    call message('set_timestep: nsteps = ',val_int=nsteps, info=2, &
                  gal_id=current_gal_id)

    if (nsteps > p_nsteps_max) then
      call message('set_timestep: Maximum number of timesteps ', &
                   val_int=p_nsteps_max, msg_end=' reached.', &
                   info=2, gal_id=current_gal_id)
      nsteps = p_nsteps_max
      success = .false.
    endif
  end function set_timestep

  subroutine read_input_parameters(gal_id)
    ! Reads the input parameters file to RAM
    use iso_fortran_env
    use IO
    integer, intent(in) :: gal_id
    integer :: i

    if (.not.allocated(galaxy_data)) then
      allocate(galaxy_data(number_of_redshifts,number_of_columns))
      allocate(output_times(number_of_redshifts))
      output_times = 0
      ! Reads the output times
      call IO_read_dataset_scalar('t', gal_id, output_times, &
                                  nrows=number_of_redshifts)
    endif

    ! Saves current galaxy identifier
    current_gal_id = gal_id
    ! Resets galaxy data array
    galaxy_data(:,:) = -1

    ! Reads time dependent parameters for this galaxy
    call IO_read_dataset_scalar('r_disk', gal_id, galaxy_data(:,1))
    call IO_read_dataset_scalar('v_disk', gal_id, galaxy_data(:,2))
    call IO_read_dataset_scalar('r_bulge', gal_id, galaxy_data(:,3))
    call IO_read_dataset_scalar('v_bulge', gal_id, galaxy_data(:,4))
    call IO_read_dataset_scalar('r_halo', gal_id, galaxy_data(:,5))
    call IO_read_dataset_scalar('v_halo', gal_id, galaxy_data(:,6))
    call IO_read_dataset_scalar('nfw_cs1', gal_id, galaxy_data(:,7))
    call IO_read_dataset_scalar('Mgas_disk', gal_id, galaxy_data(:,8))
    call IO_read_dataset_scalar('Mstars_disk', gal_id, galaxy_data(:,9))
    call IO_read_dataset_scalar('SFR', gal_id, galaxy_data(:,10))
    call IO_read_dataset_scalar('Mhalo', gal_id, galaxy_data(:,11))
    call IO_read_dataset_scalar('Mstars_bulge', gal_id, galaxy_data(:,12))

    ! Determines the maximum radius for this galaxy over the whole history
    r_max_kpc_history = maxval(galaxy_data(:,1)) * p_rmax_over_rdisk

    ! Sets n1, maximum number of snapshots
    n1 = number_of_redshifts

    ! Sets the initial and final valid snapshots (uses the disk size as a marker)
    init_it = -1
    max_it = number_of_redshifts
    do i=1,n1
      if ((galaxy_data(i,1)>=p_rdisk_min) .and. (init_it<0)) then
        init_it = i
      endif
      if ((galaxy_data(i,1)<0) .and. (init_it>0)) then
        max_it = i-1
        exit
      endif
    enddo
    iread = init_it-1

  endsubroutine read_input_parameters


  subroutine reset_input_params()
    ! Resets the reading of the input parameters file
    iread = 0
    last_output = .false.
  end subroutine reset_input_params


  subroutine set_input_params(gal_id, error)
    ! Reads dimensional input parameters that must be specified and may vary
    ! from galaxy to galaxy and from snapshot to snapshot
    use messages, only: error_message
    integer, intent(in) :: gal_id
    logical, intent(out) :: error
    double precision :: next_time_input
    double precision :: current_time_input

    error = .false.

    ! Reads the whole file on first access
    if (gal_id /= current_gal_id ) then
      call read_input_parameters(gal_id)
    endif

    iread=iread+1

    ! Traps the case where there is no valid galaxy data or
    ! there is a problem reading the parameters
    if (iread<0 .or. max_it-init_it==0) then
      error = .true.
      call error_message('set_input_params', &
                         'Error while reading data for galaxy.', &
                         gal_id=gal_id, info=1, code='p')
      return
    endif

    current_time_input = output_times(iread)

    if ( iread < max_it ) then
      next_time_input = output_times(iread+1)
      time_between_inputs = next_time_input-current_time_input
      t = 0 ! At each snapshot, reset the time variable
    else
      ! Repeats time_between_inputs of the previous iteration
      ! (by not overwriting it) and tags as last snapshot
      last_output = .true.
    endif

    t_Gyr   = current_time_input
    r_disk  = galaxy_data(iread,1)
    v_disk  = galaxy_data(iread,2)
    r_bulge = galaxy_data(iread,3)
    v_bulge = galaxy_data(iread,4)
    r_halo  = galaxy_data(iread,5)
    v_halo  = galaxy_data(iread,6)
    nfw_cs1 = galaxy_data(iread,7)
    Mgas_disk = galaxy_data(iread,8)
    Mstars_disk = galaxy_data(iread,9)
    SFR = galaxy_data(iread,10)
    Mhalo = galaxy_data(iread,11)
    Mstars_bulge = galaxy_data(iread,12)
    ! Temporarily setting v_sol_kms to the turbulent speed
    v_sol_kms = p_ISM_sound_speed_km_s * p_ISM_kappa
    l_sol_kpc = p_ISM_turbulent_length

!     DIMENSIONLESS PARAMETERS THAT MUST BE SPECIFIED BUT WILL NOT NORMALLY VARY FROM GALAXY TO GALAXY
!     DIFFUSIVE MAGNETIC HELICITY FLUX (DEFAULTS)
    R_kappa=         1.0d0 !Ratio kappa_t/eta_t of turbulent diffusivities of alpha_m and B

  endsubroutine set_input_params
end module input_params
!*****************************************************
module calc_params
  !Contains parameters that are calculated from the input parameters
  use units
  use math_constants
  use input_params
  use grid
  use input_constants
!
  implicit none
!
  double precision :: etat_sol,etat_sol_kmskpc,etat_sol_cm2s,td_sol,td_sol_kpcskm,td_sol_Gyr,td_sol_s,om0_kmskpc, &
                      Ur_sol_kms,r_sol,n_sol,r_n,h_sol,r_h,l_sol,r_l,v_sol,r_v,Uz_sol,r_Uz, &
                      Ur_sol,om0,r_om,r1
  contains
    subroutine set_calc_params
!     DIMENSIONLESS PARAMETERS THAT CAN BE CALCULATED OR THAT ARE NOT NORMALLY VARIED:
!     NUMERICAL
      lambda=h0_kpc/r_max_kpc  !Typical aspect ratio of disk


  endsubroutine set_calc_params
end module calc_params


