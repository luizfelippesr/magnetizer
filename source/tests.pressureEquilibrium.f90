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
! Contains a test program for the ISM pressure equilibrium calculations
program testPressureEquilibrium
  use global_input_parameters
  use  iso_fortran_env
  use PressureEquilibrium
  use messages, only:str2dbl
  implicit none
!   double precision, dimension(:,:), allocatable :: all_roots
  double precision, dimension(:), allocatable :: Rm, Sigma_d, Sigma_star
  double precision, dimension(:), allocatable :: r, h, rho
  double precision :: rdisk, Mgas, Mstars, dx
  character(len=100) :: command_argument
  integer :: i

  call get_command_argument(1, command_argument)
  ! If --help is detected, prints help information
  if (trim(command_argument)=='--help' .or. &
                  trim(command_argument)=='-h' .or. &
                  len_trim(command_argument)==0 ) then
    call get_command_argument(0, command_argument)
    print *, 'Usage:  ',trim(command_argument), ' <rdisk> <Mgas> <Mstars> [input_parameters_file]'
    print *,
    print *, 'rdisk  -> half-mass radius of the disk (kpc)'
    print *, 'Mgas   -> gas mass of the disk (Msun)'
    print *, 'Mstars -> stellar mass of the disk (Msun)'
    stop
  endif
  rdisk = str2dbl(command_argument)
  call get_command_argument(2, command_argument)
  Mgas = str2dbl(command_argument)
  call get_command_argument(3, command_argument)
  Mstars = str2dbl(command_argument)

  ! Tries to read the parameter filename from the command argument
  call get_command_argument(4, command_argument)
  if (len_trim(command_argument) == 0) then
    ! Uses example parameter file if nothing was found
    print *, '# No parameter file provided. Using standard parameters.'
  else
    ! Uses specified parameter file
    call read_global_parameters(trim(command_argument))
    print *, '# Using global parameters file: '// trim(command_argument)
  endif

  ! Prints input information
  print *, '# Using rdisk =', rdisk
  print *, '# Using Mgas =', Mgas
  print *, '# Using Mstars =', Mstars
  print *, '#'
  print *, '#'

  ! Allocates the arrays
  allocate(r(p_nx_ref))
  allocate(h(p_nx_ref))
  allocate(rho(p_nx_ref))
  allocate(Sigma_star(p_nx_ref))
  allocate(Sigma_d(p_nx_ref))
  allocate(Rm(p_nx_ref))
!   allocate(all_roots(p_nx_ref,3))

  ! Prepares r
  dx = p_rmax_over_rdisk*rdisk/(p_nx_ref-1)
  do i=1, p_nx_ref
    r(i) = (i-1)*dx
  enddo

  ! Solves for it
  call solve_hydrostatic_equilibrium_numerical(rdisk, Mgas, Mstars, r, 0d0*r, &
                                                Sigma_d*0d0, Sigma_d*0d0, &
                                                rho, h, &
                                                Sigma_star, Sigma_d, Rm)
  ! Prints the solutions
  write (*,*) '#','       r ','|      rho ','|        h ', '|Sigma_star','| Sigma_gas', '|        Rm |'
  do i=1,p_nx_ref


    write(*,'(ES11.4, ES11.4, ES11.4, ES11.4, ES11.4, ES11.4)') r(i), rho(i), h(i), Sigma_star(i), Sigma_d(i), Rm(i)
  enddo
end program testPressureEquilibrium
