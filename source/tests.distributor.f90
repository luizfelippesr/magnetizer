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
module test_dist
  implicit none

  contains

  subroutine print_and_exit(gal_id, switch, error)
    use messages
    integer, intent(in) :: gal_id
    logical, intent(in) :: switch
    logical, intent(out) :: error
    call message('Fake work', gal_id=gal_id, info=1)
    error = .false.
  end subroutine print_and_exit

end module test_dist
program test_distributor
  use mpi
  use input_params
  use global_input_parameters
  use IO
  use messages
  use test_dist
  use jobs
  implicit none

  integer :: igal, nmygals, igal_first, igal_last, igal_finished
  integer :: info_mpi
  integer, parameter :: master_rank = 0
  integer, allocatable, dimension(:) :: galaxies_list
  character(len=100) :: command_argument
  integer :: i, j, iproc
  logical :: lstop, error
  logical :: lforce = .false.
  logical :: lsingle_galaxy_mode = .false.
  logical :: start_galaxy = .false.
  logical :: lresuming_run = .false.
  character(len=100) :: date
  double precision :: tstart,tfinish
  integer :: rank, nproc, ierr, rc, ncycles, flush_signal
  integer, parameter :: finished_tag = 0
  integer, parameter :: newjob_tag = 17
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer,dimension(8) :: time_vals


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

  call jobs_prepare(.false.,.true.)
  galaxies_list = jobs_reads_galaxy_list(2)
  call jobs_distribute(print_and_exit, .true., galaxies_list)

end program test_distributor
