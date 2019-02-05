module LoSintegrate_aux
  ! Auxiliary functions used by the program LoSintegrate
  implicit none
  private
  public :: Compute_RM, Compute_Stokes_I

  contains

  pure function Compute_RM(Bpara, ne, x_path, z_path)
    ! Computes the Faraday rotation measure for one specific line of sight
    !
    ! Input: Bpara -> 1d-array, B parallel to the LoS, in microgauss
    !        ne -> 1d-array, number of thermal electrons, in cm^-3
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: RM, in rad m^-2
    double precision, dimension(:), intent(in) :: Bpara, ne, x_path, z_path
    double precision, dimension(size(ne)) :: integrand
    double precision :: Compute_RM

    integrand = 0.812 * Bpara * ne * 1d3
    ! obs: last factor converts from kpc to pc

    Compute_RM = Integrator(integrand, x_path, z_path)
  end function Compute_RM

  pure function Compute_Stokes_I(Bperp, ncr, x_path, z_path, wavelength, alpha)
    ! Computes Stokes parameter I along a line of sight
    !
    ! Input: Bperp -> 1d-array, B parallel to the LoS, in microgauss
    !        ne -> 1d-array, number of thermal electrons, in cm^-3
    !        wavelength -> real, wavelength, in meters
    !        alpha -> real, spectral index of the CR electrons, default: 3
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: total synchrotron emissivity in arbitrary units
    double precision, dimension(:), intent(in) :: Bperp, ncr, x_path, z_path
    double precision, dimension(size(ncr)) :: integrand
    double precision, intent(in) :: wavelength
    double precision, optional, intent(in) :: alpha
    double precision :: Compute_Stokes_I

    if (present(alpha)) then
      integrand = emissivity(Bperp, ncr, wavelength, alpha)
    else
      integrand = emissivity(Bperp, ncr, wavelength)
    endif

    Compute_Stokes_I = Integrator(integrand, x_path, z_path)
  end function Compute_Stokes_I

  pure function emissivity(Bperp, ncr, wavelength, alpha)
    ! Synchrotron emissivity
    !
    ! Input: Bperp -> 1d-array, magnitude of B perpendicular to LoS, in microgauss
    !        ncr -> 1d-array, number of CR electrons, in cm^-3
    !        wavelength -> real, wavelength, in meters
    !        alpha -> real, spectral index of the CR electrons, default: 3
    ! Output: emissivity -> 1d-array, in arbitrary units
    !
    double precision, dimension(:), intent(in) :: Bperp, ncr
    double precision, intent(in) :: wavelength
    double precision, optional, intent(in) :: alpha
    double precision, dimension(size(ncr)) :: emissivity
    double precision :: aa

    if (present(alpha)) then
      aa = alpha
    else
      aa = 3.0
    endif

    emissivity = Bperp**((alpha+1.0)/2d0) * wavelength**((alpha-1.0)/2d0)
  end function emissivity


  pure function Integrator(integrand, x_path, z_path) result(integral)
    ! Basic trapezoidal rule integrator (used for LOS integration)
    double precision, dimension(:), intent(in) :: integrand, x_path, z_path
    double precision :: dr, integral
    integer :: i
    integral = 0

    do i=1, size(integrand)-1
        dr = sqrt((x_path(i+1)-x_path(i))**2d0 + (z_path(i+1)-z_path(i))**2d0)
        integral = integral + 0.5*(integrand(i)+integrand(i+1))*dr
    end do
  end function Integrator

end module LoSintegrate_aux
