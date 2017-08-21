module floor_field
  ! Contains subroutines and functions related to the calculation of the
  ! floor of the large scale galactic magnetic field
  use global_input_parameters

  implicit none
  private

  public initialize_floor_sign, update_floor_sign
  public compute_floor_target_field, compute_floor_source_coefficient

  double precision, public, protected :: Bfloor_sign
  double precision :: t_last_sign_choice
contains
  function compute_floor_target_field(r,l,h,Beq,Delta_r) result (B_floor)
    use grid, only: lambda
    ! Computes the target floor field
    double precision, intent(in), dimension(:) :: r, l, h, Beq
    double precision, intent(in) :: Delta_r
    double precision, dimension(size(r)) :: B_floor
    double precision, dimension(size(r)) :: Ncells, brms
    integer :: i

    !Number of turbulent cells in the annular volume
    Ncells = abs(3.d0*r*Delta_r*h/l**3/lambda**2)
    !Small-scale magnetic field strength
    brms = fmag*Beq
    !Floor magnetic field
    B_floor = exp(-Delta_r/2./abs(r))*brms/Ncells**(1d0/2d0)*l/Delta_r*lambda/3
    B_floor = B_floor * Bfloor_sign * C_floor

  end function compute_floor_target_field

  function compute_floor_source_coefficient(h,etat,Uz,Dyn_gen) result (A_floor)
    use math_constants
    use input_constants
    ! Produces the target floor field
    double precision, intent(in), dimension(:) ::  h, etat, Uz, Dyn_gen
    double precision, dimension(size(h)) :: A_floor
    double precision, dimension(size(h)) :: R_U, Dyn_crit
    integer :: i

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
      call random_number(Bfloor_sign)
      Bfloor_sign = sign(1d0,Bfloor_sign-0.5d0)
      t_last_sign_choice = time
    endif
  end subroutine

  subroutine initialize_floor_sign()
    ! Initializes the magnetic field floor sign
      t_last_sign_choice = 0d0
      ! Initializes the floor sign randomly
      call random_number(Bfloor_sign)
      Bfloor_sign = sign(1d0,Bfloor_sign-0.5d0)
  end subroutine initialize_floor_sign

end module floor_field


