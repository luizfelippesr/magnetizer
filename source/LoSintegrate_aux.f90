module LoSintegrate_aux
  ! Auxiliary functions used by the program LoSintegrate
  implicit none
  private
  public :: Compute_RM, Compute_Stokes_I,LoSintegrate
  public :: Galaxy_Properties, LoS_data
  public :: alloc_Galaxy_Properties

  type Galaxy_Properties
    double precision, allocatable, dimension(:,:) :: Br, Bp, Bz
    double precision, allocatable, dimension(:,:) :: Rcyl, h, n
    integer :: n_redshifts, n_grid
  end type

  type LoS_data
    double precision, allocatable, dimension(:) :: Stokes_I
    double precision, allocatable, dimension(:) :: RM
    double precision, allocatable, dimension(:) :: number_of_cells
    double precision :: wavelength = 20d-2 ! 20 cm or 1.49 GHz
    double precision :: alpha = 3d0
    double precision :: theta = 0d0
    logical :: B_scale_with_z = .false.
  end type

  contains

  subroutine LoSintegrate(props, impact_y, impact_z, data)
    type(Galaxy_Properties), intent(in) :: props
    type(LoS_data), intent(inout) :: data
    double precision, intent(in) :: impact_y ! kpc, Impact parameter in y-direction
    double precision, intent(in) :: impact_z  ! kpc, Impact parameter relative to z-direction

    double precision,dimension(props%n_redshifts,2*props%n_grid) :: Bpara_all, &
                                                      Bperp_all, xc, zc, ne_all
    logical, dimension(props%n_redshifts,2*props%n_grid) :: valid
    double precision, allocatable, dimension(:) :: x_path, z_path
    double precision, allocatable, dimension(:) :: Bpara, Bperp, ne
    double precision, dimension(props%n_redshifts) :: Bx, By, Bz, Bmag
    double precision, dimension(props%n_redshifts) :: angle_B_LoS, tmp
    integer, dimension(2*props%n_grid) :: js
    integer :: i, j, it

    ! Allocates everything
    valid = .false.

    ! Works on all redshifts simultaneously, looping over the radial grid-points
    ! i -> index for quantities in a cartesian box
    ! j -> index for axi-symmetric quantities
    do i=1,2*props%n_grid
      ! Constructs coordinates and auxiliary indices js
      if (i<props%n_grid+1) then
        ! Index for behind the x-z plane
        j=props%n_grid-i+1
        js(i) = j
        ! sign for x coordinate
        xc(:,i) = -1
      else
        ! Index ahead of the x-z plane
        j = i-props%n_grid
        js(i) = j
        ! sign for x coordinate
        xc(:,i) = 1
      end if

      ! The available values for x are x_i = sqrt(R_i^2-b^2)
      tmp = props%Rcyl(:,j)**2-impact_y**2 ! NB allocated on-the-fly: Fortran2003 feature
      where (tmp>=0)
        xc(:,i) = xc(:,i) * sqrt(tmp)
        ! Stores a mask to be used with the Fortran 2003 'pack' intrinsic function
        valid(:,i) = .true.
      endwhere

      ! Sets z-coord
      zc(:,i) = xc(:,i) / tan(data%theta) + impact_z

      ! Bx = Br * x/R - Bp * y/R  (NB allocated on-the-fly)
      Bx = props%Br(:,j) * xc(:,i)/props%Rcyl(:,j) - props%Bp(:,j)* impact_y/props%Rcyl(:,j)
      ! By = Br * y/R + Bp * x/R (NB allocated on-the-fly)
      By = props%Br(:,j) * impact_y/props%Rcyl(:,j) - props%Bp(:,j)* xc(:,i)/props%Rcyl(:,j)
      ! Bz = Bzmod * ? (NB allocated on-the-fly)
      Bz = props%Bz(:,j)

      ! Simple vector calculations
      ! |B|
      Bmag = sqrt(Bx**2 + By**2 + Bz**2)
      ! B_\parallel = dot(B,n), where n is the LoS direction
      ! [NB  n = (cos(theta),0,sin(theta))  ]
      Bpara_all(:,i) = Bx*cos(data%theta) + Bz*sin(data%theta)
      ! angle = arccos(Bpara/|B|)
      angle_B_LoS = acos(Bpara_all(:,i)/Bmag)
      ! B_\perp = |B|*sin(angle) -- magnitude of the perpendicular component
      Bperp_all(:,i) = Bmag*sin(angle_B_LoS)

      tmp = abs(zc(:,i))/props%h(:,j)

      ! Scales the density with z (coordinate)
      ne_all(:,i) = props%n(:,j) * exp(-tmp)

      if (data%B_scale_with_z) then
        ! Scale with z (coordinate)
        Bperp_all(:,i) = Bperp_all(:,i) * exp(-tmp)
        Bpara_all(:,i) = Bpara_all(:,i) * exp(-tmp)
      else
        ! Constant for |z|<h, zero otherwise
        where (tmp>1)
          Bperp_all(:,i) = 0d0
          Bpara_all(:,i) = 0d0
        endwhere
      endif
    enddo

    if (.not.allocated(data%RM)) then
      allocate(data%RM(props%n_redshifts))
      allocate(data%Stokes_I(props%n_redshifts))
      allocate(data%number_of_cells(props%n_redshifts))
    endif

    ! Now work is done for each redshift (as the valid section of each array may
    ! may be different)
    do it=1,props%n_redshifts
      ! Filters away invalid parts of the arrays
      Bpara = pack(Bpara_all(it,:),valid(it,:))
      Bperp = pack(Bperp_all(it,:),valid(it,:))
      ne = pack(ne_all(it,:),valid(it,:))
      x_path = pack(xc(it,:),valid(it,:))
      z_path = pack(zc(it,:),valid(it,:))

      data%number_of_cells(it) = size(x_path)
      ! Synchrotron emission
      ! NB Using the total density as a proxi for cosmic ray electron density
      data%Stokes_I(it) = Compute_Stokes_I(Bperp, ne, x_path, z_path, &
                                           data%wavelength, data%alpha)
      ! Faraday rotation (backlit)
      ! NB Using the total density as a proxi for thermal electron density
      data%RM(it) = Compute_RM(Bpara, ne, x_path, z_path)
    enddo
  end subroutine LoSintegrate


  pure function Compute_RM(Bpara, ne, x_path, z_path)
    ! Computes the Faraday rotation measure for one specific line of sight
    !
    ! Input: Bpara -> 1d-array, B parallel to the LoS, in microgauss
    !        ne -> 1d-array, number of thermal electrons, in cm^-3
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: RM, in rad m^-2
    !
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

  subroutine alloc_Galaxy_Properties(nz, nr, galprops)
    integer, intent(in) :: nz, nr
    type(Galaxy_Properties), intent(inout) :: galprops
    allocate( galprops%Br(nz,nr),   &
              galprops%Bp(nz,nr),   &
              galprops%Bz(nz,nr),   &
              galprops%Rcyl(nz,nr), &
              galprops%h(nz,nr)     )
    galprops%n_redshifts = nz
    galprops%n_grid = nr
  end subroutine alloc_Galaxy_Properties

end module LoSintegrate_aux
