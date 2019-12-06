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
module floor_field
  ! Contains subroutines and functions related to the calculation of the
  ! floor of the large scale galactic magnetic field
  use global_input_parameters
  use random, only: random_sign
  implicit none
  private

  public initialize_floor_sign, update_floor_sign
  public compute_floor_target_field, compute_floor_source_coefficient

  double precision, public, protected :: Bfloor_sign
  double precision :: t_last_sign_choice
  double precision, allocatable, dimension(:) :: sign_array
  logical :: refresh_sign_array
contains
  function compute_floor_target_field(r,l,h,Beq,Delta_r) result (B_floor)
    use grid, only: lambda, check_allocate
    ! Computes the target floor field
    double precision, intent(in), dimension(:) :: r, l, h, Beq
    double precision, intent(in) :: Delta_r
    double precision, dimension(size(r)) :: B_floor
    double precision, dimension(size(r)) :: Ncells, brms
    double precision :: next_layer
    integer :: i

    !Number of turbulent cells in the annular volume
    Ncells = abs(3.d0*r*Delta_r*h/l**3/lambda**2)
    !Small-scale magnetic field strength
    brms = fmag*Beq
    !Floor magnetic field
    B_floor = exp(-Delta_r/2./abs(r))*brms/Ncells**(1d0/2d0)*l/Delta_r*lambda/3
    B_floor = B_floor * Bfloor_sign * C_floor

    if (p_space_varying_floor) then
      ! If called for the first time and/or the grid was extended, refreshes
      if (.not.allocated(sign_array)) then
        refresh_sign_array = .true.
      else if (size(sign_array) /= size(B_floor)) then
        refresh_sign_array = .true.
      endif

      ! If marked to refresh (usually done by update_floor_sign)
      if (refresh_sign_array) then
        ! First, check whether the sign array has the correct shape
        call check_allocate(sign_array, size(B_floor))
        ! Initializes with Bfloor_sign
        sign_array = Bfloor_sign
        ! Go through the layers of width Delta_r, switching the signs on each
        ! layer
        next_layer = Delta_r
        do i=1, size(r)
          if (r(i) > next_layer) then
            next_layer = next_layer + Delta_r ! increments layer
            sign_array(1:i-1) = sign_array(1:i-1) * random_sign() ! saves
          endif
        enddo
        ! No need to refresh in the next call
        refresh_sign_array = .false.
      endif

      B_floor = B_floor * sign_array
    endif
  end function compute_floor_target_field

  function compute_floor_source_coefficient(h,etat,Uz,Dyn_gen) result (A_floor)
    use math_constants
    use input_constants
    ! Produces the target floor field
    double precision, intent(in), dimension(:) ::  h, etat, Uz, Dyn_gen
    double precision, dimension(size(h)) :: A_floor
    double precision, dimension(size(h)) :: R_U, Dyn_crit

    !Dimensionless outflow parameter
    R_U = Uz*h/etat
    !Estimate of critical dynamo number
    Dyn_crit = -(pi/2)**5*(1.d0 +4*C_U*R_U/pi**2)**2
    !Coefficient of Bfloor for source term to get Bp=Bfloor
    A_floor = (pi/2)**2*etat/h**2*(1.d0+4*C_U*R_U/pi**2)*(1.d0-Dyn_gen/Dyn_crit)
  end function compute_floor_source_coefficient

  subroutine update_floor_sign(time, tau)
    ! Alternates the sign of the floor magnetic field
    !
    ! Input: time -> present time
    !        taus -> an array containig the values of tau for the galaxy
    !                in the same units as 'time'
    double precision, intent(in) :: time
    double precision, dimension(:), intent(in) :: tau

    if (time-t_last_sign_choice > p_floor_kappa*minval(tau)) then
      Bfloor_sign = random_sign()
      t_last_sign_choice = time
      ! Sets the sign array (for spatially varying floor) to refresh mode
      refresh_sign_array = .true.
    endif
  end subroutine

  subroutine initialize_floor_sign()
    ! Initializes the magnetic field floor sign
    t_last_sign_choice = 0d0
    ! Initializes the floor sign randomly
    Bfloor_sign = random_sign()
    ! Sets the sign array (for spatially varying floor) to refresh mode
    refresh_sign_array = .true.
  end subroutine initialize_floor_sign

end module floor_field


