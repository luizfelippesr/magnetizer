! Contains subroutines which compute the pressure and density at the midplane
! and disk scaleheight assuming (magneto)hydrostatic equilibrium and accounting
! for the various pressure contributions.
module pressureEquilibrium
  use math_constants
  implicit none
  private

  procedure(press_dbl), pointer, public :: midplane_pressure
  procedure(dens_dbl), pointer, public :: midplane_density
  procedure(h_dbl), pointer, public :: scaleheight

  abstract interface
    function press_dbl(r, rdisk, Mgas, Mstars, sigma_scaling)
        import
        double precision, dimension(:), intent(in) :: r
        double precision, intent(in) ::rdisk, Mgas, Mstars
        double precision, dimension(size(r)) :: press_dbl
        logical, optional, intent(in) :: sigma_scaling
    end function press_dbl
  end interface

  abstract interface
    function dens_dbl(r, P, B, cs, gamma, csi, vturb)
        import
        double precision, dimension(:), intent(in) :: r
        double precision, dimension(:), intent(in) :: P, B
        double precision, dimension(:), intent(in), optional :: vturb
        double precision, intent(in) :: cs, gamma, csi
        double precision, dimension(size(r)) :: dens_dbl
    end function dens_dbl
  end interface

  abstract interface
    function h_dbl(r, rdisk, Mgas, rho)
        import
        double precision, dimension(:), intent(in) :: r, rho
        double precision, intent(in) :: rdisk, Mgas
        double precision, dimension(size(r)) :: h_dbl
    end function h_dbl
  end interface

  public set_density_procedures

contains

  subroutine set_density_procedures(densityProcedure, pressureProcedure)
    ! Sets which function will be used to compute the midplane pressure
    ! and which function is used to compute the midplane density
    character(len=*), intent(in) :: densityProcedure, pressureProcedure

    select case (trim(pressureProcedure))
    case('simple'); midplane_pressure => midplane_pressure_simple
    case default; stop 'set_density_procedures: Unknown pressureProcedure'
    end select

    select case (trim(densityProcedure))
    case('simple'); midplane_density => midplane_density_simple
    case default; stop 'set_density_procedures: Unknown densityProcedure'
    end select

    ! Other options may be included later
    scaleheight => scaleheight_simple

  end subroutine set_density_procedures

  function midplane_pressure_simple(r, rdisk, Mgas, Mstars,sigma_scaling) result(P)
    ! Computes the pressure in midplane assuming that stars and gas
    ! follow an exponential disks and have the same vertical distribution
    ! Input: r -> radii array
    !        rdisk -> half mass radius of the disk, in kpc
    !        Mg -> gas mass of the disk, in solar masses
    !        Mstars -> stellar mass of the disk, in solar masses
    ! Output: array containing the pressure in Gaussian units
    use input_constants
    use fgsl
    double precision, dimension(:), intent(in) :: r
    double precision, intent(in) :: rdisk, Mgas, Mstars
    double precision, dimension(size(r)) :: P, sigma_star_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, parameter :: convertPressureSItoGaussian=10
    double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
    double precision, parameter :: G_SI = FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
    double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
    double precision, parameter :: stellarHeightToRadiusScale = 1d0/7.3
    double precision :: sigma_gas_SI = 10d3 ! m/s
    double precision :: constant
    logical, optional, intent(in) :: sigma_scaling
    logical :: sigma_scaling_actual

    if (present(sigma_scaling)) then
      sigma_scaling_actual = sigma_scaling
    else
      sigma_scaling_actual = .false.
    endif

    if (sigma_scaling) then
      sigma_star_SI = sqrt( G_SI/2d0 &
                        * stellarHeightToRadiusScale &
                        * rs_to_r50**2/rdisk/kpc_SI  &
                        * Mstars * Msun_SI &
                        * exp(- abs(r) * rs_to_r50/rdisk ))

      where (sigma_star_SI < sigma_gas_SI)
        sigma_star_SI = sigma_gas_SI
      endwhere
    else
      sigma_star_SI = sigma_gas_SI
    endif

    P = G_SI/8d0/pi &
              * Mgas*(Mgas+ (sigma_gas_SI/sigma_star_SI)*Mstars) * Msun_SI**2 &
              * (rs_to_r50/rdisk/kpc_SI)**4      &
              * exp(- 2d0 * abs(r) * rs_to_r50/rdisk ) * convertPressureSItoGaussian

  end function midplane_pressure_simple

  function midplane_density_simple(r, P, B, cs, gamma, csi, vturb) result(rho)
    ! Computes the density in the midplane under the assumption of
    ! c_s = v_{turb} = constant
    ! Input: r -> radii array (in kpc)
    !        P -> total pressure in the midplane, erg/cm^3
    !        B -> array containing the large scale magnetic field B(r), in microGauss
    !        cs -> sound speed, in km/s
    !        gamma -> adiabatic index
    !        csi -> ratio between turbulent pressure and turbulent magnetic
    !               field pressure
    ! Output: array containing the densities in Gaussian units
    double precision, dimension(:), intent(in) :: r
    double precision, dimension(:), intent(in) :: P, B
    double precision, dimension(:), intent(in), optional :: vturb
    double precision, intent(in) :: cs, gamma, csi
    double precision, dimension(size(r)) :: rho
    double precision, parameter :: KMS_TO_CMS=1d5
    integer :: i

    if (present(vturb)) &
      print *, 'midplane_density_simple: warning, vturb will be ignored.'

    ! Computes the density, taking into account the large scale field pressur e
    rho = P/(cs*KMS_TO_CMS)**2/(csi+1d0/3d0+1d0/gamma)
!     rho = (P - (B/1d6)**2/4d0/pi)/(cs*KMS_TO_CMS)**2/(csi+1d0/3d0+1d0/gamma)

    if ( any(rho <= 0.0)) &
      print *, 'midplane_density_simple: warning, non-positive densities.'

  end function midplane_density_simple

  function scaleheight_simple(r, rdisk, Mgas, rho) result(height)
    ! Computes the pressure in midplane assuming that stars and gas
    ! follow an exponential disks and have the same vertical distribution
    ! Input: r -> radii array
    !        rdisk -> half mass radius of the disk, in kpc
    !        Mg -> gas mass of the disk, in solar masses
    !        Mstars -> stellar mass of the disk, in solar masses
    ! Output: array containing the pressure in Gaussian units
    use input_constants
    use fgsl
    double precision, dimension(:), intent(in) :: r, rho
    double precision, intent(in) :: rdisk, Mgas
    double precision, dimension(size(r)) :: Sigma, height
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, parameter :: unit_conversion=FGSL_CONST_MKSA_SOLAR_MASS*1d3 &
                      / (FGSL_CONST_MKSA_PARSEC*1e5)**3 ! 1 Msun/kpc3 in g/cm3
    integer :: i

    Sigma = Mgas/2d0/pi*(rs_to_r50/rdisk)**2  * exp(- abs(r) * rs_to_r50/rdisk)
    height = Sigma / rho * unit_conversion ! result in kpc

    if ( any(height <= 0.0)) &
      print *, 'midplane_density_simple: warning, non-positive heights.'

    do i=1,size(height)
      if (height(i)>5)  &
        print *, 'r', r(i), 'h', height(i), rho(i)/Hmass
    end do

  end function scaleheight_simple


end module pressureEquilibrium
