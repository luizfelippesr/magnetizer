!# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
!#
!# This file is part of Magnetizer.
!#
!# Magnetizer is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# Magnetizer is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
!#
module Observables_aux
  use LoSintegrate_aux
  use messages
  implicit none
  private
  public Compute_I_PI_RM, set_runtype, write_observables
  type(LoS_data),public :: gbl_data
  logical :: lI = .false., lPI = .false., lU=.false., lQ=.false.
  logical :: lLoS=.false., lFRB = .false., ldust = .false.
  integer, parameter :: nRMs = 10
  double precision, parameter :: INVALID = -1d50
  ! Used to decide whether a galaxy should be neglected for being to small
  double precision, public, parameter :: minimum_disc_radius=0.15 ! kpc

  double precision, allocatable, dimension(:,:) :: ts_U, ts_Q, ts_I, ts_PI
  double precision, allocatable, dimension(:,:) :: ts_U_err, ts_Q_err, ts_I_err, ts_PI_err
  double precision, allocatable, dimension(:,:) :: ts_SM, ts_DM, ts_counts
  double precision, allocatable, dimension(:,:) :: ts_RM, ts_column
  double precision, allocatable, dimension(:,:) :: ts_theta, ts_z, ts_y
  double precision :: impact_y, impact_z
  character(len=5) :: prefix
  character(len=10) :: tag
contains
  subroutine set_runtype(run_type)
    character(len=50) :: run_type
    select case (trim(run_type))
      case ('intervening')
        lLoS = .true. ! Line of sight observables: RM, DM, SM
      case ('FRB_host')
        lFRB = .true.
        lLoS = .true. ! Line of sight observables: RM, DM, SM
      case ('synchrotron')
        lPI = .true.
        lQ = .true.
        lU = .true.
        lI = .true.
      case ('dust')
        ldust = .true.
        lQ = .true.
        lU = .true.
        lI = .true.
    end select

    prefix = ''
    if (lFRB) prefix = 'FRB_'

  end subroutine set_runtype

  function Compute_I_PI_RM(gal_id, resume, error) result(runtime)
    use messages
    use IO
    use math_constants
    use global_input_parameters
    use random
    integer, intent(in) :: gal_id
    logical :: random_theta
    logical, intent(in) :: resume
    logical, intent(out) :: error
    double precision :: runtime, err
    double precision :: cpu_time_start,cpu_time_finish
    type(Galaxy_Properties) :: props
    double precision, allocatable, dimension(:,:) :: buffer, ts_buffer
    double precision, allocatable, dimension(:) :: bufferz
    integer :: i, iz, iRM, LoS_counter

    error = .false.

    call cpu_time(cpu_time_start)


    ! Sets spectral index of the CR energy distribution
    gbl_data%alpha = p_obs_CR_alpha
    gbl_data%dust_alpha = p_obs_dust_alpha
    gbl_data%B_scale_with_z = p_obs_scale_with_z
    gbl_data%ignore_small_scale_field = p_obs_ignore_small_scale

    ! If an specific theta is provided (for testing) use it
    if (gbl_data%theta > 100) then
      ! greater than 100 signals random choice of theta
      random_theta = .true.
    else
      random_theta = .false.
    endif

    ! Character tag to be added for integrated observables
    if (.not.ldust) then
      tag = str(gbl_data%wavelength*100, 2)//'cm'
    else
      tag = 'dust'
    endif


    ! Reads everything
    props%n_RMs = nRMs
    props%igal = gal_id

    ! Reads galaxy properties from hdf5 file
    call alloc_Galaxy_Properties(number_of_redshifts,p_nx_ref, props)
    allocate(buffer(number_of_redshifts,p_nx_ref))
    allocate(bufferz(number_of_redshifts))
    call IO_read_dataset_vector('r', gal_id, buffer, group='Output')
    props%Rcyl = buffer
    call IO_read_dataset_vector('Br', gal_id, buffer, group='Output')
    props%Br = buffer
    call IO_read_dataset_vector('Bp', gal_id, buffer, group='Output')
    props%Bp = buffer
    call IO_read_dataset_vector('Bzmod', gal_id, buffer, group='Output')
    props%Bz = buffer
    call IO_read_dataset_vector('h', gal_id, buffer, group='Output')
    props%h = buffer/1d3 ! Converts from pc to kpc
    call IO_read_dataset_vector('n', gal_id, buffer, group='Output')
    props%n = buffer
    call IO_read_dataset_scalar('z', gal_id, bufferz, group='Input')
    props%z = bufferz
    call IO_read_dataset_scalar('Mgas_disk', gal_id, bufferz, group='Input')
    props%Mgas_disk = bufferz
    call IO_read_dataset_scalar('Mstars_disk', gal_id, bufferz, group='Input')
    props%Mstars_disk = bufferz
    call IO_read_dataset_scalar('r_disk', gal_id, bufferz, group='Input')
    props%r_disk = bufferz


    if (lLoS) then
      call prepare_ts(ts_RM, number_of_redshifts, props%n_RMs, &
                      trim(prefix)//'RM', gal_id)
      call prepare_ts(ts_SM, number_of_redshifts, props%n_RMs, &
                      trim(prefix)//'SM', gal_id)
      call prepare_ts(ts_DM, number_of_redshifts, props%n_RMs, &
                      trim(prefix)//'DM', gal_id)

      call prepare_ts(ts_column, number_of_redshifts, props%n_RMs, &
                      trim(prefix)//'column_density', gal_id)

      call prepare_ts(ts_y, number_of_redshifts, props%n_RMs, &
                      trim(prefix)//'LoS_y', gal_id)
      call prepare_ts(ts_z, number_of_redshifts, props%n_RMs, &
                      trim(prefix)//'LoS_z', gal_id)
      call prepare_ts(ts_counts, number_of_redshifts, 1, trim(prefix)//'LoS_counts', gal_id)
    endif

    if (lU) call prepare_ts(ts_U, number_of_redshifts, 1, 'U_'//trim(tag), gal_id)
    if (lU) call prepare_ts(ts_U_err, number_of_redshifts, 1, 'U_'//trim(tag)//'_err', gal_id)

    if (lQ) call prepare_ts(ts_Q, number_of_redshifts, 1, 'Q_'//trim(tag), gal_id)
    if (lQ) call prepare_ts(ts_Q_err, number_of_redshifts, 1, 'Q_'//trim(tag)//'_err', gal_id)

    if (lI) call prepare_ts(ts_I, number_of_redshifts, 1, 'I_'//trim(tag), gal_id)
    if (lI) call prepare_ts(ts_I_err, number_of_redshifts, 1, 'I_'//trim(tag)//'_err', gal_id)

    if (lPI) call prepare_ts(ts_PI, number_of_redshifts, 1, 'PI_'//trim(tag), gal_id)
    if (lPI) call prepare_ts(ts_PI_err, number_of_redshifts, 1, 'PI_'//trim(tag)//'_err', gal_id)

    ! Any option will have an inclination angle
    call prepare_ts(ts_theta, number_of_redshifts, props%n_RMs, &
                    trim(prefix)//'theta', gal_id)

    call message('Computing observables', gal_id=gal_id, info=3)
    do i=1, size(p_obs_redshift_indices)
      if (.not.p_obs_use_all_redshifts .or. lLoS) then
        ! If required, only a select set of redshifts will be run
        iz = p_obs_redshift_indices(i)
        ! negative values in p_obs_redshift_indices mean "skip-me"
      endif

      ! Catches invalid cases!
      if (iz<1 .or. iz>number_of_redshifts) cycle ! invalid redshifts
      if (props%h(iz,2)<=0d0) cycle ! Invalid scaleheights
      if (props%r_disk(iz)<minimum_disc_radius) cycle  ! Invalid discs
      if (props%Rcyl(iz,1)<0d0) cycle ! Invalid coordinates

      ! Chooses the random seed. This done combining galaxy id and redshift to
      ! ensure reproducibility (sometimes, one may run only a select set of
      ! redshifts and still should get the same result
      call set_random_seed(gal_id, p_random_seed+iz)

      if (random_theta) then
        ! Sets theta, (automatically) only if it was not previously set
        call set_random_theta(gbl_data%theta, ts_theta(iz,1))
      endif

      ! ------ Integrated observables -----------------------------------------
      if (lI) then
        if (ts_I(iz,1)<-1d30) then
          call message('calculating I', gal_id=gal_id, val=props%z(iz), info=2)
          ts_I(iz,1) = IntegrateImage('I', props, gbl_data,iz,dust=ldust, &
                                      error=ts_I_err(iz,1))
        endif
      endif

      if (lPI) then
        if (ts_PI(iz,1)<-1d30) then
          call message('calculating PI', gal_id=gal_id, val_int=iz, info=2)
          ts_PI(iz,1) = IntegrateImage('PI', props, gbl_data,iz,dust=ldust, &
                                       error=ts_PI_err(iz,1))
        endif
      endif

      if (lQ) then
        if (ts_Q(iz,1)<-1d30) then
          call message('calculating Q (integrated)', gal_id=gal_id, val_int=iz, info=2)
          ts_Q(iz,1) = IntegrateImage('Q', props, gbl_data,iz,dust=ldust, &
                                      error=ts_Q_err(iz,1))
        endif
      endif

      if (lU) then
        if (ts_U(iz,1)<-1d30) then
          call message('calculating U (integrated)', gal_id=gal_id, val_int=iz, info=2)
          ts_U(iz,1) = IntegrateImage('U', props, gbl_data,iz,dust=ldust, &
                                      error=ts_U_err(iz,1))
        endif
      endif

      ! ------ LoS observables -----------------------------------------
      if (lLoS) then
        if (ts_counts(iz,1)>1) cycle
        LoS_counter = 0
        call message('calculating RM/DM/SM', gal_id=gal_id,  val_int=iz, info=2)
        do iRM=1, props%n_RMs
          ! Tries different angles and lines of sights until something
          ! useful is intercepted
          do
            error=.false.
            LoS_counter = LoS_counter + 1
            if (iRM/=1) then
              if (random_theta) call set_random_theta(gbl_data%theta, ts_theta(iz,iRM))
            endif
            ! Picks up a random line of sight betwen 0 and 3/2 maximum radius
            ! NB this is done on the plane of the sky!
            call set_random_val(impact_y, ts_y(iz,iRM), 10d0)
            call set_random_val(impact_z, ts_z(iz,iRM), 10d0)

            ! Converts from the plane of the sky into actual z
            impact_z = impact_z/sin(gbl_data%theta)
            call message('  Theta', gal_id=gal_id,  val=gbl_data%theta, info=4)
            call message('  impact_y', gal_id=gal_id,  val=impact_y, info=4)
            call message('  impact_z', gal_id=gal_id,  val=impact_z*1000, info=4)
            ! Integrates
            call LoSintegrate(props, impact_y, impact_z, gbl_data, iz,     &
                              I_out=.false., Q_out=.false., U_out=.false., &
                              RM_out=.true., DM_out=.true., SM_out=.true., &
                              iRM=iRM, FRB_mode=lFRB, error=error)
            if (.not.error) then
              ! All fine, exit loop
              exit
            else
              ! Bogus LoS, reset
              ts_theta(iz,iRM) = INVALID
              ts_y(iz,iRM) = INVALID; ts_z(iz,iRM) = INVALID
            endif
          enddo
        enddo
        ts_RM(iz,:) = gbl_data%RM
        ts_SM(iz,:) = gbl_data%SM
        ts_DM(iz,:) = gbl_data%DM
        ts_column(iz,:) = gbl_data%column_density
        ts_counts(iz,1) = LoS_counter
      endif
    enddo
    call cpu_time(cpu_time_finish)
    runtime = cpu_time_finish - cpu_time_start
    call message('Finished observables after ', runtime,  gal_id=gal_id, &
                 msg_end='s  CPU time', info=1)
  end function Compute_I_PI_RM

  subroutine write_observables(gal_id, runtime)
    use data_transfer

    integer, intent(in) :: gal_id
    double precision, intent(in) :: runtime


    call send_or_write_dataset(trim(prefix)//'theta', gal_id, ts_theta, units='radians', &
                                 description='Inclination (for observables calculation)')

    if (lI) then
      print *, shape(ts_I(:,1)), shape(ts_I(:,2))
      call send_or_write_dataset('I_'//trim(tag), &
                                 gal_id, ts_I(:,1), units='arbitrary', &
                                 description='Integrated synchrotron emission')
      call send_or_write_dataset('I_'//trim(tag)//'_err', &
                                 gal_id, ts_I_err(:,1), units='arbitrary', &
                                 description='Error in integrated synchrotron emission')
    endif

    if (lPI) then
      call send_or_write_dataset('PI_'//trim(tag), &
                                 gal_id, ts_PI(:,1), units='arbitrary', &
                                 description='Integrated polarised synchrotron emission')
      call send_or_write_dataset('PI_'//trim(tag)//'_err', &
                                 gal_id, ts_PI_err(:,1), units='arbitrary', &
                                 description='Error in integrated poloarised synchrotron emission')
    endif

    if (lQ) then
      call send_or_write_dataset('Q_'//trim(tag), &
                                 gal_id, ts_Q(:,1), units='arbitrary', &
                                 description='Integrated synchrotron Stokes Q')
      call send_or_write_dataset('Q_'//trim(tag)//'_err', &
                                 gal_id, ts_Q_err(:,1), units='arbitrary', &
                                 description='Error in integrated synchrotron Stokes Q')
    endif

    if (lU) then
      call send_or_write_dataset('U_'//trim(tag), &
                                 gal_id, ts_U(:,1), units='arbitrary', &
                                 description='Integrated synchrotron Stokes U')
      call send_or_write_dataset('U_'//trim(tag)//'_err', &
                                 gal_id, ts_U_err(:,1), units='arbitrary', &
                                 description='Error in integrated synchrotron Stokes U')
    endif

    if (lLoS) then
      call send_or_write_dataset(trim(prefix)//'RM', gal_id, ts_RM, units='rad/m^2', &
                                 description='Rotation measure along a random LoS')
      call send_or_write_dataset(trim(prefix)//'DM', gal_id, ts_DM, units='pc cm^-3', &
                                   description='Dispersion measure along a random LoS')
      call send_or_write_dataset(trim(prefix)//'SM', gal_id, ts_SM, units='kpc m^-{20/3}', &
                                 description='Scattering measure along a random LoS')
      call send_or_write_dataset(trim(prefix)//'column_density', gal_id, ts_column,  &
                                 units='cm^-2', description='Column density of warm neutral gas')
      call send_or_write_dataset(trim(prefix)//'LoS_y', gal_id, ts_y, &
                                 description='Impact parameter used in the RM calculation in units of rmax')
      call send_or_write_dataset(trim(prefix)//'LoS_z', gal_id, ts_z, &
                                 description='Impact parameter used in the RM calculation in units of rmax')
      call send_or_write_dataset(trim(prefix)//'LoS_counts', gal_id, ts_counts(:,1), &
                                 description='Number of sighlines used')
    endif

    call send_end_message

  end subroutine write_observables

  subroutine prepare_ts(ts, nz, n, dataset_name, gal_id)
    use IO
    double precision, dimension(:,:), allocatable, intent(inout) :: ts
    integer, intent(in) :: nz, n, gal_id
    character(len=*), intent(in) :: dataset_name

    if (.not.allocated(ts)) allocate(ts(nz,n))

    if (IO_dataset_exists(dataset_name, 'Output')) then
      if (n>1) then
        call IO_read_dataset_vector(dataset_name, gal_id, ts, 'Output')
      else
        call IO_read_dataset_scalar(dataset_name, gal_id, ts(:,1), group='Output')
      endif
    else
      ts = INVALID
    endif
  end subroutine prepare_ts

  subroutine set_random_theta(glb_val, ts_val)
    use math_constants
    double precision, intent(inout) :: glb_val, ts_val
    ! Unless a fixed angle is signaled, selects a random inclination
    ! from a uniform distribution between -90 and 90 degrees
    call set_random_val(glb_val, ts_val, pi * 0.5d0 ) ! From -90 to 90
  end subroutine set_random_theta

  subroutine set_random_val(glb_val, ts_val, max_val)
    double precision, intent(inout) :: glb_val, ts_val
    double precision, intent(in) :: max_val
    ! Unless a fixed angle is signaled, selects a random inclination
    ! from a uniform distribution between -90 and 90 degrees
    call random_number(glb_val)

    glb_val = glb_val*2-1d0 ! -1 to 1
    glb_val = glb_val * max_val ! From -max_val to max_val

    ! If ts_val contains a valid value (from a previous run),
    ! substitutes glb_val by it, otherwise, stores glb_val in ts_val
    if (ts_val<-1d20) then
      ts_val = glb_val
      print *, 'Using the glb_val'
    else
      print *, 'Using the ts_val', ts_val
      glb_val = ts_val
    endif
  end subroutine set_random_val

end module Observables_aux

program Observables
  use mpi
  use input_params
  use global_input_parameters
  use IO
  use messages
  use jobs
  use Observables_aux
  implicit none
  character(len=300) :: command_argument
  integer, allocatable, dimension(:) :: galaxies_list
  logical :: random_theta


  ! Prints welcome message and prepares to distribute jobs
  ! By default, only previously run galaxies will be selected
  call jobs_prepare(completed=.true., incomplete=.false.)
  ! Reads the command arguments
  ! Sets the type of run (all, RM, PI, I or PI/I)
  call get_command_argument(2, command_argument)
  call set_runtype(command_argument)
  ! Sets the wavelength in m
  call get_command_argument(3, command_argument)
  gbl_data%wavelength = str2dbl(command_argument)
  ! Tries to read a list of galaxy numbers from argument 4 onwards
  galaxies_list = jobs_reads_galaxy_list(4)
  ! Computes I and PI for the galaxies in the sample
  ! gbl_data%theta = 45 * 3.14156295358d0/180d0
  ! random_theta = .false.
  random_theta = .true.
  call jobs_distribute(Compute_I_PI_RM, write_observables, random_theta, &
                       galaxies_list)
end program Observables
