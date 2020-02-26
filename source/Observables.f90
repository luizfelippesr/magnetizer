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
  logical :: lRM=.false., lI=.false., lPI=.false., lFRB=.false., lDM=.false., lSM=.false.
  integer, parameter :: nRMs = 10
  double precision, parameter :: INVALID = -1d50
  ! Used to decide whether a galaxy should be neglected for being to small
  double precision, public, parameter :: minimum_disc_radius=0.15 ! kpc

  double precision, allocatable, dimension(:) :: ts_I, ts_PI, ts_counts
  double precision, allocatable, dimension(:,:) :: ts_SM, ts_DM
  double precision, allocatable, dimension(:,:) :: ts_RM, ts_column
  double precision, allocatable, dimension(:,:) :: ts_theta, ts_z, ts_y
  double precision :: impact_y, impact_z

contains
  subroutine set_runtype(run_type)
    character(len=50) :: run_type
    select case (trim(run_type))
      case ('all')
        lRM = .true.
        lPI = .true.
        lI = .true.
      case ('RM')
        lRM = .true.
      case ('DM')
        lDM = .true.
      case ('FRBint')
        lRM = .true.
        lDM = .true.
        lSM = .true.
      case ('FRBhost')
        lFRB = .true.
        lRM = .true.
        lDM = .true.
        lSM = .true.
      case ('I/PI', 'I+PI','PI/I','PI+I')
        lPI = .true.
        lI = .true.
      case ('I')
        lPI = .true.
      case ('PI')
        lPI = .true.
    end select
  end subroutine

  function Compute_I_PI_RM(gal_id, random_theta, error) result(runtime)
    use messages
    use IO
    use math_constants
    use global_input_parameters
    use random
    integer, intent(in) :: gal_id
    logical, intent(in) :: random_theta
    logical, intent(out) :: error
    double precision :: runtime
    double precision :: cpu_time_start,cpu_time_finish
    type(Galaxy_Properties) :: props
    double precision, allocatable, dimension(:,:) :: buffer
    double precision, allocatable, dimension(:) :: bufferz
    integer :: iz, iRM, LoS_counter

    call cpu_time(cpu_time_start)
    error = .false.
    ! Reads galaxy properties from hdf5 file
    call alloc_Galaxy_Properties(number_of_redshifts,p_nx_ref, props)
    allocate(buffer(number_of_redshifts,p_nx_ref))
    allocate(bufferz(number_of_redshifts))

    ! Sets spectral index of the CR energy distribution
    gbl_data%alpha = p_obs_CR_alpha
    gbl_data%dust_alpha = p_obs_dust_alpha
    gbl_data%B_scale_with_z = p_obs_scale_with_z
    gbl_data%ignore_small_scale_field = p_obs_ignore_small_scale

    props%igal = gal_id
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

    props%n_RMs = nRMs

    if (.not.allocated(ts_RM)) then
      allocate(ts_I(number_of_redshifts))
      allocate(ts_PI(number_of_redshifts))
      allocate(ts_RM(number_of_redshifts, props%n_RMs))
      allocate(ts_SM(number_of_redshifts, props%n_RMs))
      allocate(ts_DM(number_of_redshifts, props%n_RMs))
      allocate(ts_column(number_of_redshifts, props%n_RMs))
      allocate(ts_theta(number_of_redshifts, props%n_RMs))
      allocate(ts_y(number_of_redshifts, props%n_RMs))
      allocate(ts_z(number_of_redshifts, props%n_RMs))
      allocate(ts_counts(number_of_redshifts))
    endif

    call message('Computing observables', gal_id=gal_id, info=1)
    do iz=1, number_of_redshifts
      ! Catches invalid redshifts (markers for invalid runs)
      if (props%h(iz,2)<=0d0 .or. props%r_disk(iz)<minimum_disc_radius) then
        ! Marks them as invalid (will be converted into NaN by the python API)
        ts_I(iz) = INVALID
        ts_PI(iz) = INVALID
        ts_theta(iz,:) = INVALID
        ts_y(iz,:) = INVALID
        ts_z(iz,:) = INVALID
        ts_RM(iz,:) = INVALID
        ts_DM(iz,:) = INVALID
        ts_SM(iz,:) = INVALID
        ts_column(iz,:) = INVALID
        ts_counts(iz) = INVALID
        cycle
      endif

      ! First, randomizes inclination and impact parameter
      ! Unless a fixed angle is signaled, selects a random inclination
      ! from a uniform distribution between -90 and 90 degrees
      call set_random_seed(gal_id, p_random_seed)
      if (random_theta) then
        call random_number(gbl_data%theta)
        gbl_data%theta = gbl_data%theta*2-1d0 ! -1 to 1
        gbl_data%theta = gbl_data%theta * pi * 0.5d0 ! From -90 to 90
      endif
      ! The first value of theta is the one used for integrated I and PI.
      ts_theta(iz,1) = gbl_data%theta

      if (props%Rcyl(iz,1)>0d0) then
        if (lI) then
          call message('calculating I', gal_id=gal_id, val_int=iz, info=2)
          ts_I(iz) = IntegrateImage('I', props, gbl_data,iz)
        endif
        if (lPI) then
          call message('calculating PI', gal_id=gal_id, val_int=iz, info=2)
          ts_PI(iz) = IntegrateImage('PI', props, gbl_data,iz)
        endif
        if (lRM .or. lDM .or. lSM .or. lDM) then
          LoS_counter = 0
          call message('calculating RM/DM/SM', gal_id=gal_id,  val_int=iz, info=2)
          do iRM=1, props%n_RMs
            do
              LoS_counter = LoS_counter + 1
              ! Tries different angles and lines of sights until
              if (iRM/=1) then
                if (random_theta) then
                  call random_number(gbl_data%theta)
                  gbl_data%theta = gbl_data%theta*2-1d0 ! -1 to 1
                  gbl_data%theta = gbl_data%theta * pi * 0.5d0 ! From -90 to 90
                endif
                ts_theta(iz,iRM) = gbl_data%theta
              endif
              error=.false.
              ! Picks up a random line of sight betwen 0 and 3/2 maximum radius
              ! NB this is done on the plane of the sky!
              call random_number(impact_y)
              call random_number(impact_z)
              ! from -1 to 1
              impact_y = impact_y*2-1d0
              impact_z = impact_z*2-1d0
              ! from -3/2 to 3/2
              impact_y = impact_y*3/2d0
              impact_z = impact_z*3/2d0
              ts_y(iz,iRM) = impact_y
              ts_z(iz,iRM) = impact_z
              ! Converts from the plane of the sky into actual z
              impact_z = impact_z/sin(gbl_data%theta)
              call message('  Theta', gal_id=gal_id,  val=gbl_data%theta, info=4)
              call message('  impact_y', gal_id=gal_id,  val=impact_y, info=4)
              call message('  impact_z', gal_id=gal_id,  val=impact_z, info=4)
              ! Integrates
              call LoSintegrate(props, impact_y, impact_z, gbl_data, iz,     &
                                I_out=.false., Q_out=.false., U_out=.false., &
                                iRM=iRM, FRB_mode=lFRB,       &
                                RM_out=lRM, DM_out=lDM, SM_out=lSM, &
                                error=error)
              if (.not.error) exit
            enddo
          enddo
          if (lRM) ts_RM(iz,:) = gbl_data%RM
          if (lSM) ts_SM(iz,:) = gbl_data%SM
          if (lDM) ts_DM(iz,:) = gbl_data%DM
          ts_column(iz,:) = gbl_data%column_density
          ts_counts(iz) = LoS_counter
        endif
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

      if (.not.lFRB) then
      if (lI) &
        call send_or_write_dataset('I_'//str(gbl_data%wavelength*100, 2)//'cm', &
                                   gal_id, ts_I, units='arbitrary', &
                                   description='Integrated synchrotron emission')
      if (lPI) &
        call send_or_write_dataset('PI_'//str(gbl_data%wavelength*100, 2)//'cm', &
                                   gal_id, ts_PI, units='arbitrary', &
                                   description='Integrated polarised synchrotron emission')
      if (lRM) then
        call send_or_write_dataset('RM', gal_id, ts_RM, units='rad/m^2', &
                                   description='Rotation measure along a random LoS')
        call send_or_write_dataset('column_density', gal_id, ts_column, units='cm^-2', &
                                   description='Column density of warm neutral gas')
        call send_or_write_dataset('RM_LoS_y', gal_id, ts_y, &
            description='Impact parameter used in the RM calculation in units of rmax')
        call send_or_write_dataset('RM_LoS_z', gal_id, ts_z, &
            description='Impact parameter used in the RM calculation in units of rmax')
        call send_or_write_dataset('RM_counts', gal_id, ts_counts, &
                                   description='Number of sighlines used')
      endif
      if (lDM) then
        call send_or_write_dataset('DM', gal_id, ts_DM, units='pc cm^-3', &
                                   description='Dispersion measure along a random LoS')
      endif
      if (lSM) then
        call send_or_write_dataset('SM', gal_id, ts_SM, units='kpc m^-{20/3}', &
                                   description='Scattering measure along a random LoS')
      endif
      call send_or_write_dataset('theta', gal_id, ts_theta, units='radians', &
                                 description='Inclination (for observables calculation)')
    else
      ! In FRB mode, write things to a different path
      call send_or_write_dataset('FRB_RM', gal_id, ts_RM, units='rad/m^2', &
                                 description='Rotation measure along a random LoS')
      call send_or_write_dataset('FRB_RM_LoS_y', gal_id, ts_y, &
          description='Impact parameter used in the RM calculation in units of rmax')
      call send_or_write_dataset('FRB_RM_LoS_z', gal_id, ts_z, &
          description='Impact parameter used in the RM calculation in units of rmax')
      call send_or_write_dataset('FRB_theta', gal_id, ts_theta, units='radians', &
                                description='Inclination (for FRB calculation)')
    endif

    call send_end_message
  end subroutine write_observables
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
  !gbl_data%theta = 30.d0 * 3.14156295358d0/180d0
  !random_theta = .false.
  random_theta = .true.
  call jobs_distribute(Compute_I_PI_RM, write_observables, random_theta, &
                       galaxies_list)
end program Observables
