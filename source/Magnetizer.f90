!# Magnetizer
!# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
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
!*****************************************************
program magnetizer
  use mpi
  use dynamo
  use input_params
  use global_input_parameters
  use IO
  use messages
  use jobs
  implicit none
  integer, allocatable, dimension(:) :: galaxies_list

  ! Prints welcome message and prepares to distribute jobs
  ! By default, previously run galaxies will be skipped
  call jobs_prepare(completed=.false., incomplete=.true.)
  ! Tries to read a list of galaxy numbers from argument 2 onwards
  galaxies_list = jobs_reads_galaxy_list(2)
  ! Runs the Magnetizer
  call jobs_distribute(dynamo_run, p_no_magnetic_fields_test_run, galaxies_list)

end program magnetizer
