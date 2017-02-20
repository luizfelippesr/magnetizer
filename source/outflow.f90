! Contains subroutines which compute the outflows
module outflow
  implicit none
  private
  public outflow_speed
contains
  function outflow_speed(r, rho, h, vt, rdisk, vdisk, SFR, Rm, outflow_type)
    ! Computes the outflow speed.
    ! There are currently 5 prescriptions available:
    !   1) no_outflow -> vout=0
    !   2) vturb -> assumes is always vout=vt
    !   3) wind -> compute vout from Galform's mass outflow rate
    !   4) superbubble -> uses the superbubble formalism
    !   5) superbubble_simple -> the same, without the \rho_{hot}/\rho step
    !
    ! Input: r -> radii array, in kpc
    !        rho -> densities array, in g/cm^3
    !        h -> heights array, in kpc
    !        vt -> turbulent velocities array, in km/s
    !        rdisk -> half mass radius of the disk, in kpc
    !        vdisk -> velocity at the half mass radius of the disk, in km/s
    !        SFR -> total star formation rate of the galaxy, in (solar masses)/yr
    !        Rm -> molecular/diffuse gas ratio
    !        Mgas -> gas mass of the disk, in solar masses
    !        Mstars -> stellar mass of the disk, in solar masses
    ! Output: array containing the outflow speed at different radii

    use input_constants
    use global_input_parameters
    implicit none
    double precision, intent(in) :: rdisk, vdisk, SFR
    double precision, dimension(:), intent(in) :: r, rho, h, vt, Rm
    double precision, dimension(size(r)) :: outflow_speed, n, fm
    double precision :: constant, rs, beta
    character(len=*), optional, intent(in) :: outflow_type
    character(len=35) :: outflow_type_actual
    logical :: no_average = .true.

    if (present(outflow_type)) then
      outflow_type_actual = outflow_type
    else
      outflow_type_actual = p_outflow_type
    endif

    select case (trim(outflow_type_actual))
      case('no_outflow')
        outflow_speed = r*0.0
        return
      case('vturb')
        outflow_speed = vt
        return
      case('wind')
        ! Computes the molecular fraction
        fm = Rm/(1d0+Rm)
        ! Computes the mass loading
        beta = (vdisk/p_outflow_vhot)**p_outflow_alphahot
        ! This case assumes the outflow rates used in galform are correct
        ! This, together with a Blitz&Rosolowsky star formation law leads
        ! to:
        outflow_speed = 0.5 * p_outflow_nu0 * beta * h * fm
        return
        ! Note: at the moment, nu0 and alpha_hot and Vhot are input
        ! parameters, but later these will be read from the input hdf5 files
        ! (which in turn contains the parameters used in that particular
        ! galform run).
      case('superbubble_simple')
        no_average=.true.
      case('superbubble')
        no_average=.false.
      case default
        print *, outflow_type_actual
        stop 'outflow: Unknown outflow_type'
    end select

    ! The following is for the superbubble case

    ! First, prepares some constants
    rs = rdisk*constDiskScaleToHalfMassRatio
    n = rho/Hmass

    constant = 51.3 * (p_outflow_Lsn/1d38)**(1./3.) * p_outflow_fOB/0.7d0 &
        * (p_outflow_etaSN/9.4e-3) * (p_tOB/3d0) * (40d0/p_N_SN1OB)

    outflow_speed = constant * (rs/3.0)**(-2) * SFR
    outflow_speed = outflow_speed * n**(-1./3.)
    outflow_speed = outflow_speed * (h/0.2)**(4./3.)
    outflow_speed = outflow_speed * exp(-r/rs)
    if (no_average) return

    outflow_speed = (p_outflow_hot_gas_density/rho) * outflow_speed

  end function outflow_speed
end module outflow
