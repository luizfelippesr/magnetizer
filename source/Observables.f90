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
  private
  public Compute_I_and_PI
  type(LoS_data),public :: gbl_data

contains

  subroutine Compute_I_and_PI(gal_id, random_theta, error)
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
    double precision, dimension(number_of_redshifts) :: ts_I, ts_PI

    ! Unless a fixed angle is signaled, selects a random inclination
    ! from a uniform distribution between 0 and 90 degrees
    if (random_theta) then
      call set_random_seed(gal_id, p_random_seed)
      call random_number(gbl_data%theta)
      gbl_data%theta = gbl_data%theta * pi * 0.5d0
    endif

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

    ts_I = -99999d0; ts_PI = -99999d0
    call message('Computing observables', gal_id=gal_id, info=2)
    do iz=1, number_of_redshifts
      if (props%Rcyl(iz,1)>0d0) then
        call message('calculating I', gal_id=gal_id, val_int=iz, info=4)
        ts_I(iz) = IntegrateImage('I', props, gbl_data,iz)
        call message('calculating PI', gal_id=gal_id, val_int=iz, info=4)
        ts_PI(iz) = IntegrateImage('PI', props, gbl_data,iz)
      endif
    enddo

    call IO_write_dataset('I', gal_id, ts_I, units='arbitrary', &
                          description='Integrated synchrotron emission')
    call IO_write_dataset('PI', gal_id, ts_PI, units='arbitrary', &
                        description='Integrated polarised synchrotron emission')
    call IO_write_dataset('theta', gal_id, [gbl_data%theta], units='radians', &
                        description='Inclination')
  end subroutine Compute_I_and_PI
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
  ! Hard-coded parameters (varied for testing only)
  gbl_data%B_scale_with_z = .false.
  gbl_data%ignore_small_scale_field = .false.
  ! Tries to read a list of galaxy numbers from argument 4 onwards
  galaxies_list = jobs_reads_galaxy_list(4)
  ! Computes I and PI for the galaxies in the sample
  call jobs_distribute(Compute_I_and_PI, .true., galaxies_list)
end program Observables
