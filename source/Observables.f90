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
  public Compute_I_PI_RM, set_runtype
  type(LoS_data),public :: gbl_data
  logical :: lRM=.false., lI=.false., lPI=.false.
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
      case ('I/PI', 'I+PI','PI/I','PI+I')
        lPI = .true.
        lI = .true.
      case ('I')
        lPI = .true.
      case ('PI')
        lPI = .true.
    end select

  end subroutine

  subroutine Compute_I_PI_RM(gal_id, random_theta, error)
    use messages
    use math_constants
    use global_input_parameters
    use IO
    use random
    integer, intent(in) :: gal_id
    logical, intent(in) :: random_theta
    logical, intent(out) :: error
    type(Galaxy_Properties) :: props
    double precision, allocatable, dimension(:,:) :: buffer
    integer :: iz
    double precision, dimension(number_of_redshifts) :: ts_I, ts_PI, ts_z, ts_y, ts_theta
    double precision :: impact_y, impact_z
    error = .false.

    ! Reads galaxy properties from hdf5 file
    call alloc_Galaxy_Properties(number_of_redshifts,p_nx_ref, props)
    allocate(buffer(number_of_redshifts,p_nx_ref))
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

    call message('Computing observables', gal_id=gal_id, info=1)
    do iz=1, number_of_redshifts
      ! First, randomizes inclination and impact parameter
      ! Unless a fixed angle is signaled, selects a random inclination
      ! from a uniform distribution between 0 and 90 degrees
      if (random_theta) then
        call set_random_seed(gal_id, p_random_seed)
        call random_number(gbl_data%theta)
        gbl_data%theta = gbl_data%theta * pi * 0.5d0
      endif
      ts_theta(iz) = gbl_data%theta

      if (lRM) then
        ! Picks up a random line of sight betwen 0 and maximum radius
        ! NB this is done on the plane of the sky!
        call random_number(impact_y)
        call random_number(impact_z)
        ts_y(iz) = impact_y
        ts_z(iz) = impact_z
        ! Converts from the plane of the sky into actual z
        impact_z = impact_z/sin(gbl_data%theta)
      endif

      if (props%Rcyl(iz,1)>0d0) then
        if (lI) then
          call message('calculating I', gal_id=gal_id, val_int=iz, info=2)
          ts_I(iz) = IntegrateImage('I', props, gbl_data,iz)
        endif
        if (lPI) then
          call message('calculating PI', gal_id=gal_id, val_int=iz, info=2)
          ts_PI(iz) = IntegrateImage('PI', props, gbl_data,iz)
        endif
        if (lRM) then
          call message('calculating RM', gal_id=gal_id,  val_int=iz, info=2)
          call LoSintegrate(props, impact_y, impact_z, gbl_data, iz, &
                            RM_out=.true., I_out=.false., Q_out=.false., U_out=.false.)

        endif
      endif
    enddo

    if (lI) &
      call IO_write_dataset('I_'//str(gbl_data%wavelength*100, 2)//'cm', &
                          gal_id, ts_I, units='arbitrary', &
                          description='Integrated synchrotron emission')
    if (lPI) &
      call IO_write_dataset('PI_'//str(gbl_data%wavelength*100, 2)//'cm', &
                            gal_id, ts_PI, units='arbitrary', &
                            description='Integrated polarised synchrotron emission')

    if (lRM) &
      call IO_write_dataset('RM', gal_id, gbl_data%RM, units='rad/m^2', &
                            description='Rotation measure along a random LoS')
      call IO_write_dataset('RM_LoS_y', gal_id, ts_y, &
          description='Impact parameter used in the RM calculation in units of rmax')
      call IO_write_dataset('RM_LoS_z', gal_id, ts_z, &
          description='Impact parameter used in the RM calculation in units of rmax')

    call IO_write_dataset('theta', gal_id, ts_theta, units='radians', &
                          description='Inclination (for observables calculation)')
  end subroutine Compute_I_PI_RM
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

  ! Prints welcome message and prepares to distribute jobs
  ! By default, only previously run galaxies will be selected
  call jobs_prepare(completed=.true., incomplete=.false.)

  ! Reads the command arguments
  ! Sets the wavelength in m
  call get_command_argument(2, command_argument)
  gbl_data%wavelength = str2dbl(command_argument)
  ! Sets spectral index of the CR energy distribution
  call get_command_argument(3, command_argument)
  gbl_data%alpha = str2dbl(command_argument)
  ! Sets the type of run (all, RM, PI, I or PI/I)
  call get_command_argument(4, command_argument)
  call set_runtype(command_argument)
  ! Hard-coded parameters (varied for testing only)
  gbl_data%B_scale_with_z = .false.
  gbl_data%ignore_small_scale_field = .false.
  ! Tries to read a list of galaxy numbers from argument 4 onwards
  galaxies_list = jobs_reads_galaxy_list(5)
  ! Computes I and PI for the galaxies in the sample
  call jobs_distribute(Compute_I_PI_RM, .true., galaxies_list)
end program Observables
