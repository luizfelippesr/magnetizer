! Contains subroutines which compute the pressure and density at the midplane
! and disk scaleheight assuming (magneto)hydrostatic equilibrium and accounting 
! for the various pressure contributions.
module pressureEquilibrium
  use math_constants
  implicit none
  private

  procedure(press_dbl), pointer, public :: midplane_pressure
  procedure(dens_dbl), pointer, public :: midplane_density
  procedure(dens_dbl), pointer, public :: midplane_scaleheight
  
  abstract interface
    function press_dbl(r, rdisk, Mgas, Mstars) 
        import
        double precision, dimension(:), intent(in) :: r
        double precision, intent(in) ::rdisk, Mgas, Mstars
        double precision, dimension(size(r)) :: press_dbl
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
  
  public set_density_procedures
  
contains

  subroutine set_density_procedures(densityProcedure, pressureProcedure)
    ! Sets which function will be used to compute the midplane pressure
    ! and which function is used to compute the midplane density
    character(len=*), intent(in) :: densityProcedure, pressureProcedure
    
    select case (trim(pressureProcedure))
    case('simple'); midplane_pressure => midplane_pressure_simple
    case default; stop
    end select
    
    select case (trim(densityProcedure))
    case('simple'); midplane_density => midplane_density_simple
    case default; stop
    end select
  end subroutine set_density_procedures

  function midplane_pressure_simple(r, rdisk, Mgas, Mstars) result(P)
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
    double precision, dimension(size(r)) :: P
    double precision, dimension(size(r)) :: exp_m_r_over_rs
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, parameter :: convertPressureSItoGaussian=10
    double precision :: constant 
    
    ! G_SI * Msun_SI**2 / kpc_SI**4 (i.e. G and unit adjustments)
    constant = FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT &
              *FGSL_CONST_MKSA_SOLAR_MASS**2  &
              / FGSL_CONST_MKSA_PARSEC**4 / 1e3**4 * convertPressureSItoGaussian ! J/m^3 -> erg/cm^3
    
    P = constant/4d0*(rs_to_r50/rdisk)**2 * Mgas*(Mgas+Mstars) &
          * exp(- abs(r) * rs_to_r50/rdisk )
    
  end function midplane_pressure_simple

  function midplane_density_simple(r, P, B, cs, gamma, csi, vturb) result(rho)
    ! Computes the density in the midplane under the assumption of 
    ! c_s = v_{turb} = constant
    ! Input: r -> radii array (in kpc)
    !        P -> total pressure in the midplane, erg/cm^3
    !        B -> array containing the large scale magnetic field B(r), in Gauss
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
    
    if (present(vturb)) &
      print *, 'midplane_density_simple: warning, vturb will be ignored.'
      
    rho = (P - B**2/4d0/pi)/(cs*KMS_TO_CMS)**2/(csi+0.5d0+1d0/gamma)    
  end function midplane_density_simple
  
  
end module pressureEquilibrium
