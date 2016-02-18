! Contains subroutines which compute the outflows
module outflow
  implicit none
  private
  public outflow_speed
contains
  function outflow_speed(r, rho, h, vt, r_disk, SFR, Mgas, Mstars) result(v)
    ! Computes the pressure in midplane assuming that stars and gas
    ! follow an exponential disks and have the same vertical distribution
    ! Input: r -> radii array, in kpc
    !        rho -> densities array, in g/cm^3
    !        h -> heights array, in kpc
    !        vt -> turbulent velocities array, in km/s
    !        r_disk -> half mass radius of the disk, in kpc
    !        SFR -> total star formation rate of the galaxy, in (solar masses)/yr
    !        Mgas -> gas mass of the disk, in solar masses
    !        Mstars -> stellar mass of the disk, in solar masses
    ! Output: array containing the outflow speed at different radii

    use input_constants
    use global_input_parameters
    implicit none
    double precision, intent(in) :: r_disk, SFR, Mgas, Mstars
    double precision, dimension(:), intent(in) :: r, rho, h, vt
    double precision, dimension(size(r)) :: v, n
    double precision :: constant, rs
    logical :: no_average

    select case (trim(p_outflow_type))
      case('no_outflow')
        v = r*0.0
        return
      case('vturb')
        v = vt
        return
      case('superbubble_simple')
        no_average=.true.
      case('superbubble')
        no_average=.false.
      case default
        stop 'outflow: Unknown outflow_type'
    end select

    rs = r_disk/constDiskScaleToHalfMassRatio
    n = rho/Hmass

    constant = 51.3 * (p_outflow_Lsn/1d38)**(1./3.) * p_outflow_fOB/0.7d0 &
        * (p_outflow_etaSN/9.4e-3) * (p_tOB/3d0) * (40d0/p_N_SN1OB)

    v = constant
    v = v * (rs/3.0)**(-2) * SFR
    v = v * n**(-1./3.)
    v = v * (h/0.2)**(4./3.)
    v = v * exp(-r/rs)
    if (no_average) return

    v = (p_outflow_hot_gas_density/rho) * v

  end function outflow_speed
end module outflow
