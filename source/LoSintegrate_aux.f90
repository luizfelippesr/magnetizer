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
module LoSintegrate_aux
  ! Auxiliary functions used by the program LoSintegrate
  use math_constants
  use tools, only: linspace
  use FRB, only: get_FRB_LoS_z_position
  implicit none
  private
  public :: LoSintegrate
  public :: Galaxy_Properties, LoS_data
  public :: alloc_Galaxy_Properties
  public :: print_image, IntegrateImage

  ! 30% fixed ionisation fraction!
  double precision, parameter :: ionisation_fraction = 0.3

  type Galaxy_Properties
    double precision, allocatable, dimension(:,:) :: Br, Bp, Bz
    double precision, allocatable, dimension(:,:) :: Rcyl, h, n
    double precision, allocatable, dimension(:) :: z, Mgas_disk, Mstars_disk, r_disk
    double precision :: h_FRB = 0.1
    integer :: n_redshifts, n_grid, igal
    integer :: n_RMs = 10
  end type

  type LoS_data
    double precision :: Stokes_I
    double precision :: Stokes_Q
    double precision :: Stokes_U
    double precision :: test_Bx, test_By, test_Bz
    double precision, allocatable, dimension(:) :: RM, DM, SM
    double precision, allocatable, dimension(:) :: column_density
    double precision :: number_of_cells
    double precision, allocatable, dimension(:) :: psi0
    double precision :: wavelength = 20d-2 ! 20 cm or 1.49 GHz
    double precision :: alpha = 3d0
    double precision :: dust_alpha = 1d0
    double precision :: dust_p0 = 0.2
    ! Initial 1000 value is used as a tag to request random theta
    double precision :: theta = 1000d0
    logical :: B_scale_with_z = .false.
    logical :: ignore_small_scale_field = .false.
    integer :: nz_points = 500
  end type

  ! Global variables for the integration
  type(LoS_data) :: data_glb
  logical :: dust_glb
  type(Galaxy_Properties) :: props_glb
  integer :: iz_glb

  contains

  subroutine LoSintegrate(props, impacty_rmax, impactz_rmax, data, iz, &
                          FRB_mode, RM_out, I_out, Q_out, U_out, iRM, error, &
                          DM_out, SM_out, dust_mode, test)
    use input_constants
    use messages
    use global_input_parameters, only: p_molecularHeightToRadiusScale, p_ISM_turbulent_length
    use input_constants
    use FRB
    type(Galaxy_Properties), intent(in) :: props
    type(LoS_data), intent(inout) :: data
    ! Impact parameter in y- and z-directions in units of rmax
    double precision, intent(in) :: impactz_rmax, impacty_rmax
    logical, intent(in), optional :: I_out, Q_out, U_out, FRB_mode, test
    logical, intent(in), optional :: RM_out, DM_out, SM_out
    integer, optional, intent(in) :: iRM
    double precision, dimension(2*props%n_grid) :: xc, zc, ne_all, h_all, n_mol_all
    double precision, dimension(2*props%n_grid) :: Bx_all, By_all, Bz_all
    logical, dimension(2*props%n_grid) :: valid
    logical :: lI, lQ, lU, lRM, lDM, lSM, lI_dust, lQ_dust, lU_dust
    double precision, allocatable, dimension(:) :: x_path, z_path
    double precision, allocatable, dimension(:) :: psi, psi0, n_mol
    double precision, allocatable, dimension(:) :: Bpara, Bperp, Brnd, ne, h, Brnd_B, cos2gamma
    double precision, allocatable, dimension(:) :: Bx, By, Bz, B_perp_not_y
    double precision, allocatable, dimension(:) :: z_path_new, x_path_new
    double precision :: rmax, wavelength_gal, h_m
    double precision :: tmp, impact_y, impact_z
    double precision :: z_start, z_end, xmin, xmax
    integer, dimension(2*props%n_grid) :: js
    integer :: i, j, iz, this_RM, k
    logical, parameter :: ldense_z = .true.
    logical :: lFRB, lTest, ldust
    logical, optional, intent(inout) :: error
    double precision, parameter :: EPS=1e-6
    logical, intent(in), optional :: dust_mode


    lFRB = .false.
    if (present(FRB_mode)) lFRB = FRB_mode

    ! Allocates/initialises everything
    valid = .false.
    if (.not.allocated(data%RM)) then
      allocate(data%RM(props%n_RMs))
      allocate(data%column_density(props%n_RMs))
      data%RM = 0; data%number_of_cells = 0;
      data%Stokes_I = 0; data%Stokes_U = 0; data%Stokes_Q = 0
    endif
    ! Skips if outside the galaxy
    if (abs(impacty_rmax)>1) then
      if (present(error)) error = .true.
      return
    endif

    ! Sets optional arguments and their default values
    lRM=.false.; lI=.false.; lQ=.false.; lU=.false.
    lI_dust = .false.; lQ_dust = .false.; lU_dust = .false.
    if (present(RM_out)) lRM = RM_out

    if (present(I_out)) lI = I_out
    if (present(Q_out)) lQ = Q_out
    if (present(U_out)) lU = U_out

    ldust = .false.
    if (present(dust_mode)) then
      ldust = dust_mode
      if (ldust) then
        if (present(I_out)) lI_dust = I_out
        if (present(Q_out)) lQ_dust = Q_out
        if (present(U_out)) lU_dust = U_out
      endif
    endif

    lSM = .false.; lDM =.false.
    if (present(SM_out)) then
      lSM = SM_out
      if (.not.allocated(data%SM)) allocate(data%SM(props%n_RMs))
    endif
    if (present(DM_out)) then
      lDM = DM_out
      if (.not.allocated(data%DM)) allocate(data%DM(props%n_RMs))
    endif

    lTest = .false.
    if (present(test)) then
      lTest = test
    endif


    ! Computes maximum radius and impact parameters
    rmax = props%Rcyl(iz,props%n_grid)
    impact_y = impacty_rmax * rmax
    impact_z = impactz_rmax * rmax

    ! Computes wavelength at the frame of the galaxy
    wavelength_gal = data%wavelength/(1d0 + props%z(iz))

    ! Loops over the radial grid-points
    ! i -> index for quantities in a cartesian box
    ! j -> index for axi-symmetric quantities
    do i=1,2*props%n_grid
      ! Constructs coordinates and auxiliary indices js
      if (i<props%n_grid+1) then
        ! Index for behind the x-z plane
        j=props%n_grid-i+1
        js(i) = j
        ! sign for x coordinate
        xc(i) = -1
      else
        ! Index ahead of the x-z plane
        j = i-props%n_grid
        js(i) = j
        ! sign for x coordinate
        xc(i) = 1
      end if

      ! The available values for x are x_i = sqrt(R_i^2-b^2)
      tmp = props%Rcyl(iz,j)**2-impact_y**2
      if (tmp>=0) then
        xc(i) = xc(i) * sqrt(tmp)
        ! Stores a mask to be used with the Fortran 2003 'pack' intrinsic function
        valid(i) = .true.
      endif

      ! Sets z-coord
      zc(i) = xc(i) / tan(data%theta) + impact_z


      ! Bx = Br * x/R - Bp * y/R
      Bx_all(i) =   props%Br(iz,j) * xc(i)/props%Rcyl(iz,j)    &
                  - props%Bp(iz,j)*impact_y/props%Rcyl(iz,j)
      ! By = Br * y/R + Bp * x/R
      By_all(i) = props%Br(iz,j) * impact_y/props%Rcyl(iz,j)   &
                  + props%Bp(iz,j)* xc(i)/props%Rcyl(iz,j)
      ! Bz = Bzmod * sign(z)
      Bz_all(i) = sign(props%Bz(iz,j), zc(i))

      ! Stores density and scale height
      ! NB the scaling with of n with z is done further ahead!
      ne_all(i) = props%n(iz,j) * ionisation_fraction
      h_all(i) = props%h(iz,j)

      if (ldust) then
        ! At the moment, the following is NOT strictly speaking the molecular
        ! density, but is proportional to it
        ! NB the scaling with of n with z is done further ahead!
        n_mol_all(i) = compute_molecular_density(props%Rcyl(iz,j), zc(i), &
                                                 props%r_disk(iz),        &
                                                 props%Mgas_disk(iz),     &
                                                 props%Mstars_disk(iz) )
      endif
    enddo

    ! Now work is done over the whole grid
    ! Skips empty cells
    if (all(.not.valid(:)))  then
      if (present(error)) error = .true.
      return
    endif

    ! Filters away invalid parts of the arrays
    Bx = pack(Bx_all(:),valid(:))
    By = pack(By_all(:),valid(:))
    Bz = pack(Bz_all(:),valid(:))
    ne = pack(ne_all(:),valid(:))
    n_mol = pack(n_mol_all(:),valid(:))
    h = pack(h_all(:),valid(:))
    x_path = pack(xc(:),valid(:))
    z_path = pack(zc(:),valid(:))

    ! Requests at least 3 grid points
    if (size(z_path)<3) then
      if (present(error)) error = .true.
      return
    endif

    ! The following makes the grid denser, to allow meaningful z-integration
    ! This step is irrelevant for edge-on, but is important for face-on
    if (ldense_z) then

      if (z_path(1) < z_path(size(z_path))) then
        ! Sets a tentative z-minimum/maximum to -3.5/3.5 times the (maximum) scale-height
        z_start = -3.5d0*h(1)
        z_end = 3.5d0*h(size(h))

        ! If in dust mode, use the molecular gas as reference instead
        if (ldust) then
          z_start = max(z_start, -3d0*h_m)
          z_end = min(z_end, 3d0*h_m)
        endif

        ! Checks whether the line of sight is not too far away to be relevant
        if (z_end<z_path(1) .or. z_start>z_path(size(z_path))) then
          if (present(error)) error = .true.
          return
        endif

        ! Keeps the extremities within the interval
        z_start = max(z_start, z_path(1))
        z_end = min(z_end, z_path(size(z_path)))
        ! Adds an epsilon to avoid extrapolation
        z_end = z_end - EPS
        z_start = z_start + EPS
      else
        ! Sets a tentative z-minimum/maximum to -3.5/3.5 times the scale-height
        z_start = 3.5d0*h(1)
        z_end = -3.5d0*h(size(h))

        ! If in dust mode, use the molecular gas as reference instead
        if (ldust) then
          z_start = min(z_start, 3d0*h_m)
          z_end = max(z_end, -3d0*h_m)
        endif

        ! Checks whether the line of sight is not too far away to be relevant
        if (z_start<z_path(size(z_path)) .or. z_end>z_path(1) ) then !
          if (present(error)) error = .true.
          return
        endif

        ! Keeps the extremities within the interval
        z_end = max( z_end, z_path(size(z_path)))
        z_start = min( z_start, z_path(1))
        ! Adds an epsilon to avoid extrapolation
        z_end = z_end + EPS
        z_start = z_start - EPS
      endif

      ! Produces the dense grid
      z_path_new = linspace(z_start, z_end, data%nz_points)
      x_path_new = (z_path_new - impact_z)*tan(data%theta)

      ! If this is in FRB mode, uses the LoS path to find the FRB position
      if (lFRB) then
        ! Gets random position with probability proportional to molecular density
        z_start = get_FRB_LoS_z_position(x_path_new, impact_y, z_path_new, &
                                         props%h_FRB,                      &
                                         props%r_disk(iz),                 &
                                         props%Mgas_disk(iz),              &
                                         props%Mstars_disk(iz),            &
                                         rmax=props%Rcyl(iz,props%n_grid) )
        ! Traps cases where an FRB cannot be found in that LoS
        ! NB the extremely negative z value is a tag for invalid LoS
        if (z_start<-1d6) then
          if (present(error)) error = .true.
          return
        endif
        ! Produces a new dense grid starting at the *FRB position*
        z_path_new = linspace(z_start, z_end, data%nz_points)
        x_path_new = (z_path_new - impact_z) * tan(data%theta)
      endif

      call densify(Bx, z_path, z_path_new)
      call densify(By, z_path, z_path_new)
      call densify(Bz, z_path, z_path_new)
      call densify(ne, z_path, z_path_new)
      call densify(n_mol, z_path, z_path_new)
      call densify(h, z_path, z_path_new)

      z_path = z_path_new
      x_path = x_path_new

    endif


    ! Traps zero scaleheights
    where (h<1e-6) h=1e-6

    ! Adds the z dependence to the density
    ne = ne  * exp(-abs(z_path)/h)
    ! And molecular density
    h_m = props%r_disk(iz) * constDiskScaleToHalfMassRatio * p_molecularHeightToRadiusScale
    n_mol = n_mol * exp( -abs(z_path)/h_m )



    ! Adds the z dependence the magnetic field
    if (data%B_scale_with_z) then
      ! Scale with z (coordinate)
      Bx = Bx * exp(-abs(z_path)/h)
      By = By * exp(-abs(z_path)/h)
      Bz = Bz * exp(-abs(z_path)/h)
    else
      ! Constant for |z|<h, zero otherwise
      where (abs(z_path)/h>1)
        Bx = 0d0
        By = 0d0
        Bz = 0d0
      endwhere
    endif


    ! Allocates everything
    allocate(Bpara(size(Bx)))
    allocate(B_perp_not_y(size(Bx)))
    allocate(Bperp(size(Bx)))

    ! Simple vector calculations
    ! B_\parallel = dot(B,n), where n is the LoS direction
    ! i.e n = (sin(theta),0,cos(theta))
    Bpara = Bx*sin(data%theta) + Bz*cos(data%theta)
    ! B_\perp = dot(B,m), where m is perpendicular to the LoS direction and y
    ! i.e. m = (-cos(theta),0,sin(theta))
    B_perp_not_y = -Bx*cos(data%theta) + Bz*sin(data%theta)

    ! Magnitude of the total perpendicular field
    Bperp = sqrt(B_perp_not_y**2 + By**2)

    allocate(Brnd(size(Bpara)))
    allocate(Brnd_B(size(Bpara)))
    allocate(cos2gamma(size(Bpara)))

    ! Computes the random field in microgauss, assuming hydrogen mass and:
    !                      v0=1e5cm/s f=0.5  1d6\muG
    Brnd = sqrt(4*pi*ne*Hmass)*10d5 * 0.5d0 *1d6

    Brnd_B = Brnd/sqrt(Bperp**2 + Bpara**2 + 1e-40)

!     if (data%ignore_small_scale_field) then
!       Brnd = ne*0d0
!     else
!       ! Computes the random field in microgauss, assuming hydrogen mass and:
!       !                      v0=1e5cm/s f=0.5  1d6\muG
!       Brnd = sqrt(4*pi*ne/ionisation_fraction*Hmass)*10d5 * 0.5d0 *1d6
!
!       ! TODO change this to the values used in the simulation
!     endif

    data%number_of_cells = size(x_path)

    if (lTest) then
      data%test_Bx = Compute_Test('x', x_path, z_path, Bx, By, Bz)
      data%test_By = Compute_Test('y', x_path, z_path, Bx, By, Bz)
      data%test_Bz = Compute_Test('z', x_path, z_path, Bx, By, Bz)
      return
    endif

    ! ------------ Synchrotron ------------
    if (lI) then
      ! Synchrotron total emission
      ! NB Using the total density as a proxy for cosmic ray electron density
      data%Stokes_I = Compute_Stokes('I', Bperp, ne, x_path, z_path,&
                                     wavelength_gal, data%alpha, Brnd)
    endif
    if (lQ .or. lU) then
      ! Intrinsic polarization angle
      psi0 = pi/2d0 + atan2(By,B_perp_not_y)
      where (psi0>pi)
        psi0 = psi0-2d0*pi
      endwhere
      where (psi0<-pi)
        psi0 = psi0+2d0*pi
      endwhere
      allocate(psi(size(z_path)))
      do i=1, size(z_path)
        ! Integrates (i.e. computes RM) until the i-th layer
        psi(i) = psi0(i) &
                + (wavelength_gal)**2 * Compute_RM(Bpara(1:i),ne(1:i), &
                                                   x_path(1:i),z_path(1:i))
      enddo

      if (lQ) then
        data%Stokes_Q = Compute_Stokes('Q', Bperp, ne, x_path, z_path,&
                                       wavelength_gal, data%alpha, psi)
      endif
      if (lU) then
        data%Stokes_U = Compute_Stokes('U', Bperp, ne, x_path, z_path,&
                                       wavelength_gal, data%alpha, psi)
      endif
    endif

    ! ------------ Aligned dust grains emission ------------
    ! Intrinsic polarization angle
    if (lQ_dust .or. lU_dust) then
      psi0 = pi/2d0 + atan2(By,B_perp_not_y)
      where (psi0>pi)
        psi0 = psi0-2d0*pi
      endwhere
      where (psi0<-pi)
        psi0 = psi0+2d0*pi
      endwhere
    endif

    cos2gamma = Bperp**2/(Bpara**2 + Bperp**2 + 1d-40 ) + 1d-40

    if (lI_dust) then
      ! Synchrotron emission
      data%Stokes_I = Compute_Stokes_dust('I', n_mol, x_path, z_path,&
                                          data%dust_alpha, data%dust_p0, psi0, Brnd_B, cos2gamma)

    endif
    if (lQ_dust) then

      data%Stokes_Q = Compute_Stokes_dust('Q', n_mol, x_path, z_path,&
                                          data%dust_alpha, data%dust_p0, psi0, Brnd_B, cos2gamma)
    endif
    if (lU_dust) then
      data%Stokes_U = Compute_Stokes_dust('U', n_mol, x_path, z_path,&
                                          data%dust_alpha, data%dust_p0, psi0, Brnd_B, cos2gamma)
    endif

    ! This will be used by RM, DM and SM calculations
    if (present(iRM)) then
      this_RM = iRM
    else
      this_RM = 1
    endif

    ! ------------ RM (backlit), SM, DM, column density ------------
    if (lRM) then
      ! Faraday rotation (backlit) and (neutral warm gas) column density
      ! NB Using the total density as a proxy for thermal electron density
      data%RM(this_RM) = Compute_RM(Bpara, ne, x_path, z_path)
      ! Computes column density using neutral gas
      data%column_density(this_RM) = Compute_Column_Density( &
         ne*(1d0-ionisation_fraction)/ionisation_fraction, x_path, z_path)
    endif
    ! Dispersion measure
    if (lDM)  data%DM(this_RM) = Compute_DM(ne, x_path, z_path)
    ! Scattering measure
    if (lSM)  data%SM(this_RM) = Compute_SM(ne, h, x_path, z_path, p_ISM_turbulent_length)

  end subroutine LoSintegrate

  function compute_molecular_density(R, z, r_disk, Mgas_disk, Mstars_disk) result(n_mol)
    use global_input_parameters, only: p_molecularHeightToRadiusScale
    use surface_density, only: molecular_gas_surface_density
    use input_constants
    double precision, intent(in) :: R, z
    double precision, intent(in) :: r_disk, Mgas_disk, Mstars_disk
    double precision :: n_mol, h_m
    double precision, dimension(1) :: Sigma_m

    Sigma_m = molecular_gas_surface_density([R], r_disk, Mgas_disk, Mstars_disk)
    h_m = r_disk * constDiskScaleToHalfMassRatio * p_molecularHeightToRadiusScale
    n_mol = Sigma_m(1)
  end function compute_molecular_density

  function IntegrateImage(im_type, props,data,iz,number_of_calls,method,error,dust) result(res)
    use fgsl
    use, intrinsic :: iso_c_binding
    character(len=*), intent(in) :: im_type
    type(Galaxy_Properties), intent(in) :: props
    logical, optional, intent(in) :: dust
    type(LoS_data), intent(inout) :: data
    integer, intent(in) :: iz
    integer, intent(in), optional :: number_of_calls
    character(len=*), intent(in), optional :: method
    type(fgsl_monte_function) :: gfun
    type(fgsl_rng) :: r
    type(fgsl_rng_type) :: t
    type(fgsl_monte_plain_state) :: s
    type(fgsl_monte_miser_state) :: m
    type(fgsl_monte_vegas_state) :: v
    integer(fgsl_size_t) :: calls
    real(fgsl_double) :: xl(2), xu(2), res, err
    double precision, optional, intent(out) :: error
    type(c_ptr) :: ptr
    integer(fgsl_int) :: status
    character(len=10) :: mthd
    logical :: ldust

    dust_glb = .false.
    if (present(dust)) dust_glb = dust

    calls = 2500
    if (present(number_of_calls)) calls = number_of_calls
    mthd = 'MISER'
    if (present(method)) mthd = method

    data_glb = data
    props_glb = props
    iz_glb = iz

    xl = -1
    xu = 1

    ! Sets up FGSL stuff
    t = fgsl_rng_env_setup()
    r = fgsl_rng_alloc(t)
    if (im_type=='I') then
      gfun = fgsl_monte_function_init(IntegrandImage_I, 2_fgsl_size_t, ptr)
    elseif (im_type=='PI') then
      gfun = fgsl_monte_function_init(IntegrandImage_PI, 2_fgsl_size_t, ptr)
    elseif (im_type=='Q') then
      gfun = fgsl_monte_function_init(IntegrandImage_Q, 2_fgsl_size_t, ptr)
    elseif (im_type=='U') then
      gfun = fgsl_monte_function_init(IntegrandImage_U, 2_fgsl_size_t, ptr)
    else
      stop 'Error.'
    endif

    select case (trim(mthd))
      case('MISER')
        m = fgsl_monte_miser_alloc(2_fgsl_size_t)
        status = fgsl_monte_miser_integrate(gfun, xl, xu, 2_fgsl_size_t, &
                                            calls, r, m, res, err)
        call fgsl_monte_miser_free(m)
      case('VEGAS')
        v = fgsl_monte_vegas_alloc(2_fgsl_size_t)
        status = fgsl_monte_vegas_integrate(gfun, xl, xu, 2_fgsl_size_t, &
                                            calls, r, v, res, err)
        call fgsl_monte_vegas_free(v)
      case('plain')
        s = fgsl_monte_plain_alloc(2_fgsl_size_t)
        status = fgsl_monte_plain_integrate(gfun, xl, xu, 2_fgsl_size_t, &
                                            calls, r, s, res, err)
        call fgsl_monte_plain_free(s)
    end select

    if (present(error)) error=err
  end function IntegrateImage

  function IntegrandImage_I(v_c, n, params) bind(c)
    ! Wrapper to allow using LoSintegrate with FGSL
    use fgsl
    use, intrinsic :: iso_c_binding

    integer(c_size_t), value :: n
    type(c_ptr), value :: v_c, params
    real(c_double) :: IntegrandImage_I
    real(c_double), dimension(:), pointer :: v
    ! Reads the memory address
    call c_f_pointer(v_c, v, [n])
    call LoSintegrate(props_glb, v(1), v(2), data_glb, iz_glb,   &
                                   RM_out=.false., I_out=.true., &
                                   Q_out=.false., U_out=.false., &
                                   dust_mode=dust_glb)
    IntegrandImage_I = data_glb%Stokes_I
  end function IntegrandImage_I

  function IntegrandImage_PI(v_c, n, params) bind(c)
    ! Wrapper to allow using LoSintegrate with FGSL
    use fgsl
    use, intrinsic :: iso_c_binding

    integer(c_size_t), value :: n
    type(c_ptr), value :: v_c, params
    real(c_double) :: IntegrandImage_PI
    real(c_double), dimension(:), pointer :: v
    ! Reads the memory address
    call c_f_pointer(v_c, v, [n])
    call LoSintegrate(props_glb, v(1), v(2), data_glb, iz_glb,    &
                                   RM_out=.false., I_out=.false., &
                                   Q_out=.true., U_out=.true.,    &
                                   dust_mode=dust_glb)
    IntegrandImage_PI = sqrt(data_glb%Stokes_Q**2 + data_glb%Stokes_U**2)
  end function IntegrandImage_PI

  function IntegrandImage_Q(v_c, n, params) bind(c)
    ! Wrapper to allow using LoSintegrate with FGSL
    use fgsl
    use, intrinsic :: iso_c_binding

    integer(c_size_t), value :: n
    type(c_ptr), value :: v_c, params
    real(c_double) :: IntegrandImage_Q
    real(c_double), dimension(:), pointer :: v
    ! Reads the memory address
    call c_f_pointer(v_c, v, [n])
    call LoSintegrate(props_glb, v(1), v(2), data_glb, iz_glb,    &
                                   RM_out=.false., I_out=.false., &
                                   Q_out=.true., U_out=.false.,    &
                                   dust_mode=dust_glb)
    IntegrandImage_Q = data_glb%Stokes_Q
  end function IntegrandImage_Q

  function IntegrandImage_U(v_c, n, params) bind(c)
    ! Wrapper to allow using LoSintegrate with FGSL
    use fgsl
    use, intrinsic :: iso_c_binding

    integer(c_size_t), value :: n
    type(c_ptr), value :: v_c, params
    real(c_double) :: IntegrandImage_U
    real(c_double), dimension(:), pointer :: v
    ! Reads the memory address
    call c_f_pointer(v_c, v, [n])
    call LoSintegrate(props_glb, v(1), v(2), data_glb, iz_glb,    &
                                   RM_out=.false., I_out=.false., &
                                   Q_out=.false., U_out=.true.,    &
                                   dust_mode=dust_glb)
    IntegrandImage_U = data_glb%Stokes_U
  end function IntegrandImage_U

  function Compute_Column_Density(nn, x_path, z_path)
    ! Computes column density measure for one specific line of sight
    !
    ! Input: nn -> 1d-array, number of neutral atoms, in cm^-3
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: column density, cm^-2
    !
    double precision, dimension(:), intent(in) :: nn, x_path, z_path
    double precision :: Compute_Column_Density

    Compute_Column_Density = Integrator(nn, x_path*cm_kpc, z_path*cm_kpc)

  end function Compute_Column_Density

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

    integrand = 812 * Bpara * ne

    Compute_RM = Integrator(integrand, x_path, z_path)
  end function Compute_RM

  pure function Compute_DM(ne, x_path, z_path)
    ! Computes the Faraday rotation measure for one specific line of sight
    !
    ! Input: ne -> 1d-array, number of thermal electrons, in cm^-3
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: DM, in pc cm^-3
    !
    double precision, dimension(:), intent(in) :: ne, x_path, z_path
    double precision, dimension(size(ne)) :: integrand
    double precision :: Compute_DM

    integrand = ne

    Compute_DM = Integrator(integrand, x_path, z_path) *1d3 ! kpc to pc
  end function Compute_DM

  pure function Compute_SM(ne, h, x_path, z_path, l0)
    ! Computes the Faraday rotation measure for one specific line of sight
    !
    ! Input: ne -> 1d-array, number of thermal electrons, in cm^-3
    !        h -> 1d-array, scaleheight profile, in kpc
    !        l0 -> turbulent length in kpc
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: SM, in kpc m^-20/3
    !
    double precision, dimension(:), intent(in) :: ne, h, x_path, z_path
    double precision, intent(in) :: l0
    double precision, dimension(size(ne)) :: integrand, C2N, l
    double precision :: Compute_SM

    l = l0
    where (h>l0)
      l=l0
    elsewhere
      l=h
    end where

    integrand = 1.8d-3*(ne/0.01)**2*(l/1d-6)**(-2d0/3d0)

    Compute_SM = Integrator(integrand, x_path, z_path)
  end function Compute_SM


  pure function Compute_Stokes(S, Bperp, ncr, x_path, z_path, &
                                 wavelength, alpha, Brnd_or_psi)
    ! Computes Stokes parameter I, Q or U along a line of sight
    !
    ! Input: Bperp -> 1d-array, B parallel to the LoS, in microgauss
    !        ne -> 1d-array, number of thermal electrons, in cm^-3
    !        wavelength -> real, wavelength, in meters
    !        alpha -> real, spectral index of the CR electrons, default: 3
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    ! Output: total synchrotron emissivity in arbitrary units
    character(len=*), intent(in) :: S
    double precision, dimension(:), intent(in) :: Bperp, Brnd_or_psi, ncr, x_path, z_path
    double precision, intent(in) :: wavelength
    double precision, optional, intent(in) :: alpha
    double precision, dimension(size(ncr)) :: integrand, emissivity, psi, Brnd
    double precision :: p0
    double precision :: aa
    double precision :: Compute_Stokes

    if (present(alpha)) then
      p0 = (alpha+1d0)/(alpha+7d0/3d0)
      aa = alpha
    else
      p0 = 0.75
      aa = 3d0
    endif

    if (S=='I') then
      Brnd = Brnd_or_psi
      integrand = compute_emissivity(Bperp, ncr, wavelength, aa, Brnd)
    else
      psi = Brnd_or_psi
      ! Intrinsic polarization
      p0 = (aa+1d0)/(aa+7d0/3d0)
      ! The random field should not influence the polarized emission
      emissivity = compute_emissivity(Bperp, ncr, wavelength, aa)

      if (S=='Q') then
        integrand = p0*emissivity*cos(2d0*psi)
      else
        integrand = p0*emissivity*sin(2d0*psi)
      endif
    endif

    Compute_Stokes = Integrator(integrand, x_path, z_path)

  end function Compute_Stokes


  function Compute_Stokes_dust(S, n_mol, x_path, z_path, &
                                    alpha, p0, psi, Brnd_B, cos2gamma)
    ! Computes Stokes parameter I, Q or U along a line of sight
    !
    ! Input: n_mol -> 1d-array, density of H2 , in arbitrary units (future: cm^-3)
    !        alpha -> real, parameter of the dust emissivity
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    !        S -> character with the type of Stokes parameter to be computed
    ! Output: integrated emission in arbitrary units
    character(len=*), intent(in) :: S
    double precision, dimension(:), intent(in) :: psi, n_mol, x_path, z_path
    double precision, dimension(:), intent(in) :: Brnd_B, cos2gamma
    double precision, intent(in) :: alpha, p0
    double precision, dimension(size(x_path)) :: integrand, emissivity, p0_actual
    double precision :: Compute_Stokes_dust
    integer i
    emissivity = n_mol**alpha

    p0_actual = p0*dust_instrinsic_polarization(Brnd_B)

    if (S=='I') then
      integrand = emissivity
    else if (S=='Q') then
        integrand = p0_actual * emissivity * cos(2d0*psi) * cos2gamma
    else if (S=='U') then
        integrand = p0_actual * emissivity * sin(2d0*psi) * cos2gamma
    else
        stop 'Error in Compute_Stokes_dust'
    endif

    Compute_Stokes_dust = Integrator(integrand, x_path, z_path)

  end function Compute_Stokes_dust


  function Compute_Test(S, x_path, z_path, Bx, By, Bz)
    ! Computes Stokes parameter I, Q or U along a line of sight
    !
    ! Input: n_mol -> 1d-array, density of H2 , in arbitrary units (future: cm^-3)
    !        alpha -> real, parameter of the dust emissivity
    !        x_path,z_path -> 1d-array, positions along path, in kpc
    !        S -> character with the type of Stokes parameter to be computed
    ! Output: integrated emission in arbitrary units
    character(len=*), intent(in) :: S
    double precision, dimension(:), intent(in) :: x_path, z_path, Bx, By, Bz
    double precision, dimension(size(x_path)) :: integrand, B
    double precision :: Compute_Test, w
    double precision :: z_threshold

    if (S=='x') then
        B = Bx
    else if (S=='y') then
        B = By
    else if (S=='z') then
        B = Bz
    else
        stop 'Error...'
    endif

    z_threshold = minval(abs(z_path))

    integrand = 0d0
    where (abs(z_path)<=z_threshold)
      integrand = B
    end where

    Compute_Test = Integrator(integrand, x_path, z_path)

  end function Compute_Test



  pure function compute_emissivity(Bperp, ncr, wavelength, alpha, Brnd)
    ! Synchrotron emissivity
    !
    ! Input: Bperp -> 1d-array, magnitude of B perpendicular to LoS, in microgauss
    !        ncr -> 1d-array, number of CR electrons, in cm^-3
    !        wavelength -> real, rest frame wavelength, in meters
    !        alpha -> real, spectral index of the CR electrons, default: 3
    ! Output: emissivity -> 1d-array, in arbitrary units
    !
    double precision, dimension(:), intent(in) :: Bperp, ncr
    double precision, intent(in) :: wavelength, alpha
    double precision, dimension(:), optional, intent(in) :: Brnd
    double precision, dimension(size(ncr)) :: compute_emissivity, Bperp_tot2


    if (present(Brnd)) then
      Bperp_tot2 = Bperp**2 + Brnd**2
    else
      Bperp_tot2 = Bperp**2
    endif
    compute_emissivity = ncr * Bperp_tot2**((alpha+1.0)/4d0) * wavelength**((alpha-1.0)/2d0)

  end function compute_emissivity

  subroutine print_image(props, data, directory, ymax, zmax, nprint, isnap, &
                         dust, test)
      use messages
      use IO
      character(len=*) :: directory
      character(len=20) :: form
      integer, optional, intent(in) :: nprint, isnap
      logical, optional, intent(in) :: dust, test
      integer :: n, i, j, it, iz, it_max, it_min
      integer, dimension(7) :: unit
      double precision :: impact_y, impact_z, zmax, ymax, rmax
      real, dimension(:,:), allocatable :: RM_im, I_im, Q_im, U_im, N_im
      real, dimension(:),allocatable :: y, z
      type(Galaxy_Properties), intent(in) :: props
      type(LoS_data), intent(inout) :: data
      logical :: to_file, ldust, lTest


      n = 60; iz = 1 ! Default values
      if (present(nprint)) n = nprint
      if (present(isnap)) iz = isnap

      to_file = .true.; ldust = .false.; lTest = .false.
      if (trim(directory)=='-') to_file = .false.
      if (present(dust)) ldust = dust
      if (present(test)) lTest = test
      rmax = props%Rcyl(iz,props%n_grid)

      form ='('//str(n)//'E15.5)'

      if (to_file) then
        open(newunit=unit(1),file=trim(directory)//'I.dat', FORM='FORMATTED', status='replace')
        open(newunit=unit(2),file=trim(directory)//'Q.dat', FORM='FORMATTED', status='replace')
        open(newunit=unit(3),file=trim(directory)//'U.dat', FORM='FORMATTED', status='replace')
        open(newunit=unit(4),file=trim(directory)//'RM.dat', FORM='FORMATTED',status='replace')
        open(newunit=unit(5),file=trim(directory)//'cells.dat', FORM='FORMATTED',status='replace')
        open(newunit=unit(6),file=trim(directory)//'y.dat', FORM='FORMATTED', status='replace')
        open(newunit=unit(7),file=trim(directory)//'z.dat', FORM='FORMATTED',status='replace')
      endif


      allocate(I_im(n,n))
      allocate(Q_im(n,n))
      allocate(U_im(n,n))
      allocate(RM_im(n,n))
      allocate(N_im(n,n))
      allocate(y(n))
      allocate(z(n))

      if (to_file) then
        it_min=iz
        it_max=iz
      else
        it_min=1
        it_max=props%n_redshifts
      endif

      do it=it_min, it_max
        do j=1,n
          impact_z = -zmax + 2*zmax/dble(n)*j
          z = real(impact_z)
          do i=1,n
            impact_y =  -ymax + 2*ymax/dble(n)*i
            y(i) = real(impact_y)

            if (lTest) then
              call LoSintegrate(props, impact_y/rmax, impact_z/rmax, data, it, &
                                test=.true.)
              I_im(i,j) = real(data%test_Bz)
              data%Stokes_I = 0

              Q_im(i,j) = real(data%test_By)
              data%Stokes_Q = 0

              U_im(i,j) = real(data%test_Bx)
              data%Stokes_U = 0

            else
              if (.not.ldust) then
                  call LoSintegrate(props, impact_y/rmax, impact_z/rmax, data, it, &
                                    I_out=.true., U_out=.true., Q_out=.true., RM_out=.true.)
              else
                  call LoSintegrate(props, impact_y/rmax, impact_z/rmax, data, it, &
                                    I_out=.true., U_out=.true., Q_out=.true., &
                                    dust_mode=.true.)
              endif
              I_im(i,j) = real(data%Stokes_I)
              data%Stokes_I = 0

              Q_im(i,j) = real(data%Stokes_Q)
              data%Stokes_Q = 0

              U_im(i,j) = real(data%Stokes_U)
              data%Stokes_U = 0

              N_im(i,j) = real(data%number_of_cells)
              data%number_of_cells = 0

              RM_im(i,j) = real(data%RM(1))
              data%RM = 0
            endif
          enddo

          if (to_file) then
            write(unit(1), trim(form)) I_im(:,j)
            write(unit(2), trim(form)) Q_im(:,j)
            write(unit(3), trim(form)) U_im(:,j)
            write(unit(4), trim(form)) RM_im(:,j)
            write(unit(5), trim(form)) N_im(:,j)
            write(unit(6), trim(form)) y
            write(unit(7), trim(form)) z
          endif
        enddo
      enddo

      if (to_file) then
        do i=1,7
          close(unit(i))
        enddo
      endif

      deallocate(I_im); deallocate(Q_im); deallocate(U_im); deallocate(RM_im);
      deallocate(N_im); deallocate(y); deallocate(z)

  end subroutine print_image

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
    allocate( galprops%Br(nz,nr),      &
              galprops%Bp(nz,nr),      &
              galprops%Bz(nz,nr),      &
              galprops%Rcyl(nz,nr),    &
              galprops%h(nz,nr),       &
              galprops%z(nz),          &
              galprops%Mstars_disk(nz),&
              galprops%Mgas_disk(nz),  &
              galprops%r_disk(nz))
    galprops%n_redshifts = nz
    galprops%n_grid = nr
  end subroutine alloc_Galaxy_Properties

  subroutine densify(array, x0, x_dense)
    use interpolation
    double precision, allocatable, dimension(:), intent(in) :: x0, x_dense
    double precision, allocatable, dimension(:), intent(inout) :: array
    double precision, allocatable, dimension(:) :: new_array

    allocate(new_array(size(x_dense)))
      call interpolate(x0, array, x_dense, new_array)

    deallocate(array)
    call move_alloc(new_array, array)
  end subroutine

  function dust_instrinsic_polarization(b_B) result(p)
    use interpolation
    double precision, dimension(:), intent(in) :: b_B
    double precision, dimension(size(b_B)) :: p
    double precision, dimension(2,173) :: v
    double precision, dimension(173) :: p_ref, b_B_ref
    double precision, dimension(1) :: p_1 ! workaround
    integer u, i

    ! Reads the file
    open(newunit=u,file='data/depolar_random.dat', FORM='FORMATTED')
    read(u,*) v
    close(u)

    b_B_ref = v(1,:); p_ref = v(2,:)
    do i=1, size(b_B)
      if (b_B(i) >= maxval(b_B_ref)) then
        p(i) = 0d0
      else if (b_B(i) <= minval(b_B_ref)) then
        p(i) = 1d0
      else
        call interpolate(b_B_ref, p_ref, [b_B(i)], p_1)
        p(i) = p_1(1)
      endif
    enddo
  end function dust_instrinsic_polarization

end module LoSintegrate_aux
