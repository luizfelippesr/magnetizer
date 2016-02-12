! Contains subroutines which compute the outflows
module outflow
  implicit none
  private

contains
  function outflow(r, rho, h, vt, rdisk, SFR, Mgas, Mstars) result(v)
    ! Computes the pressure in midplane assuming that stars and gas
    ! follow an exponential disks and have the same vertical distribution
    ! Input: r -> radii array, in kpc
    !        rho -> densities array, in g/cm^3
    !        h -> heights array, in kpc
    !        vt -> turbulent velocities array, in km/s
    !        rdisk -> half mass radius of the disk, in kpc
    !        SFR -> total star formation rate of the galaxy, in (solar masses)/yr
    !        Mgas -> gas mass of the disk, in solar masses
    !        Mstars -> stellar mass of the disk, in solar masses
    ! Output: array containing the outflow speed at different radii

    use input_constants
    use global_input_parameters
    implicit none
    double precision, intent(in) :: r_disk, SFR, Mgas, Mstars
    double precision, dimension(:), intent(in) :: r, rho, h
    double precision, dimension(size(r)) :: v, n
    double precision :: constant, rs

    select case (trim(p_outflow_type))
    case('no_outflow'); v = r*0.0

    case('vturb'); v = vt

    case('superbubble_simple');
      rs = rdisk/constDiskScaleToHalfMassRatio
      n = rho/Hmass

      constant = 51.3 * (p_outflow_Lsn/1e38)**(1./3.) * p_outflow_fOB/0.7 &
        * (p_outflow_etaSN/9.4e-3) * (p_tOB/3.0) * (40/p_N_SN1OB)

      v = constant * n**(-1./3.) * n**(-1./3.) * (rs/3.0)**(-2) * SFR

      v = v * exp(-r/rs)

    case('superbubble');
      v = (outflow_hot_gas_density/rho)         &
          * outflow(r, rho, h, vt, rdisk, SFR, Mgas, Mstars, &
                    'superbubble_simple')

    case default; stop 'outflow: Unknown outflow_type'
    end select

end module outflow
