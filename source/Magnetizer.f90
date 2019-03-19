!# Magnetizer
!# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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

  ! Tries to read the parameter filename from the command argument (or --help)
!   call get_command_argument(1, command_argument)
!   ! If --help is detected, prints help information
!   if (trim(command_argument)=='--help' .or. trim(command_argument)=='-h') then
!     call get_command_argument(0, command_argument)
!     if (rank==master_rank) then
!       print *, 'Magnetizer '
!       print *,
!       print *, 'Computes ISM properties and solves mean field dynamo equation'&
!                //' for the output of a semi-analytic galaxy formation model.'
!       print *,
!       print *, 'Usage:'
!       print *, trim(command_argument), ' <input_parameters_file> [galaxy number] [-f]'
!       print *,
!       print *, 'For more details please visit: '&
!              //'https://github.com/luizfelippesr/magnetizer'
!       print *,
!     endif
!     stop
!   endif

  ! Skips previously run galaxies
  call jobs_prepare(completed=.false., incomplete=.true.)
  galaxies_list = jobs_reads_galaxy_list(2)
  print *, galaxies_list
!   stop
  call jobs_distribute(dynamo_run, p_no_magnetic_fields_test_run, galaxies_list)

end program magnetizer
