!# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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
!This program runs the code in fortran without mpi
!*****************************************************
program calldynamo
  !
  ! intent(in)
  !  info:   how much info to print during run
  !          options are 0=min, 1=standard, 2=max
  !  gal_id: id number of galaxy
  !          ranges from 1 to ~300,000
  !
  ! intent(out)
  !  flag:   status of run
  !          results are -1=run successful, 1=, 2=, 3=
  !
  ! Note: calculated physical variables (ts_t, ts_Br, ts_Bp, ts_Bzmod, ts_alp_m) are written to files, not passed
  !
  ! Format:  call dynamo_run(info, gal_id, flag)
  ! Standard parameters: (0, 1, flag), where igal is looped over from 1 to ngals
  !
  use dynamo
  use global_input_parameters

  implicit none 

  integer :: igal, flag
  character(len=32) :: command_argument

  call get_command_argument(1, command_argument)
  if (len_trim(command_argument) == 0) then
    write(*,*) 'No parameter file provided. Using standard global parameters.'
  else
    call read_global_parameters(trim(command_argument))
  endif

  ! Initializes IO  
  call IO_start(trim(path_to_input_directories) // '/output/' // &
                trim(model_name), trim(output_file_name))

  
  do igal=1,ngals
    call dynamo_run(info, igal, flag)
    if (info>1) then

      print*,'flag=',flag
      print*,''
    endif
  enddo

  ! Finalizes IO
  call IO_end()

end program calldynamo
!*****************************************************
