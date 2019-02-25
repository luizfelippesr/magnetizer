module LoSintegrate_aux
  ! Auxiliary functions used by the program LoSintegrate
  implicit none
  private
  public :: Compute_RM, Compute_Stokes_I,LoSintegrate
  public :: Galaxy_Properties, LoS_data
  public :: alloc_Galaxy_Properties
  public :: print_image, IntegrateImage



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
    integer :: nz_points = 500
  end type

  ! Global variables for the integration
  type(LoS_data) :: data_glb
  type(Galaxy_Properties) :: props_glb
  integer :: iz_glb

  contains
  subroutine LoSintegrate(props, impacty_rmax, impactz_rmax, data, &
                                       iz, RM_out, I_out)
    type(Galaxy_Properties), intent(in) :: props
    type(LoS_data), intent(inout) :: data
    ! Impact parameter in y- and z-directions in units of rmax
    double precision, intent(in) :: impactz_rmax, impacty_rmax
    logical, intent(in), optional :: RM_out, I_out
    double precision,dimension(2*props%n_grid) :: Bpara_all, &
                                              Bperp_all, xc, zc, ne_all, h_all
    logical, dimension(2*props%n_grid) :: valid
    logical :: lRM, lI
    double precision, allocatable, dimension(:) :: x_path, z_path
    double precision, allocatable, dimension(:) :: Bpara, Bperp, ne, h
    double precision, allocatable, dimension(:) :: z_path_new, x_path_new
    double precision :: Bx, By, Bz, rmax
    double precision :: tmp, impact_y, impact_z
    double precision :: zmin, zmax, xmin, xmax
    integer, dimension(2*props%n_grid) :: js
    integer :: i, j, iz
    logical, parameter :: ldense_z = .true.

    ! Allocates/initialises everything
    valid = .false.
    if (.not.allocated(data%RM)) then
      allocate(data%RM(props%n_redshifts))
      allocate(data%Stokes_I(props%n_redshifts))
      allocate(data%number_of_cells(props%n_redshifts))
      data%RM = 0; data%number_of_cells = 0; data%Stokes_I = 0
    endif
    ! Skips if outside the galaxy
    if (abs(impacty_rmax)>1) return

    lRM = .false.; lI = .true. ! Default values
    if (present(RM_out)) lRM = RM_out
    if (present(I_out)) lI = I_out

    rmax = props%Rcyl(iz,props%n_grid)
    impact_y = impacty_rmax * rmax
    impact_z = impactz_rmax * rmax

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
      Bx = props%Br(iz,j) * xc(i)/props%Rcyl(iz,j) - props%Bp(iz,j)*impact_y/props%Rcyl(iz,j)
      ! By = Br * y/R + Bp * x/R
      By = props%Br(iz,j) * impact_y/props%Rcyl(iz,j) + props%Bp(iz,j)* xc(i)/props%Rcyl(iz,j)
      ! Bz = Bzmod * sign(z)
      Bz = sign(props%Bz(iz,j), zc(i))

      ! Stores density and scale height
      ! NB the scaling with of n with z will be done later!
      ne_all(i) = props%n(iz,j)
      h_all(i) = props%h(iz,j)

      ! Simple vector calculations

      ! B_\parallel = dot(B,n), where n is the LoS direction
      ! [NB  n = (sin(theta),0,cos(theta))  ]
      Bpara_all(i) = Bx*sin(data%theta) + Bz*cos(data%theta)
      ! B_\perp = \sqrt{ (B_x - B_\parallel\sin\theta)^2 + B_y^2 +
      !                  + (B_z -B_\parallel\cos\theta)^2 }
      Bperp_all(i) = ( Bx - Bpara_all(i)*sin(data%theta) )**2   &
                         + By**2 +                                  &
                       ( Bz - Bpara_all(i)*cos(data%theta) )**2
      Bperp_all(i) = sqrt(Bperp_all(i))
    enddo

    ! Now work is done over the whole grid
    ! Skips empty cells
    if (all(.not.valid(:))) return

    ! Filters away invalid parts of the arrays
    Bpara = pack(Bpara_all(:),valid(:))
    Bperp = pack(Bperp_all(:),valid(:))
    ne = pack(ne_all(:),valid(:))
    h = pack(h_all(:),valid(:))
    x_path = pack(xc(:),valid(:))
    z_path = pack(zc(:),valid(:))

    ! The following makes the grid denser, to allow meaningful z-integration
    ! This step is irrelevant for edge-on, but is important for face-on
    if (ldense_z) then
      ! Sets a tentative z-maximum to 3.5 times the scale-height
      zmin = -3.5d0*h(1); zmax = 3.5d0*h(size(h))

      ! Computes corresponding values for x
      xmin = (zmin - impact_z)*tan(data%theta)
      xmax = (zmax - impact_z)*tan(data%theta)
      ! Forces x coordinate to live within the original grid
      if (abs(xmin) > abs(x_path(1))) then
        xmin = x_path(1)
        zmin = z_path(1)
      endif
      if (abs(xmax) > abs(x_path(size(x_path)))) then
        xmax = x_path(size(x_path))
        zmax = z_path(size(x_path))
      endif

      z_path_new = linspace(zmin,zmax, data%nz_points)
      x_path_new = (z_path_new - impact_z)*tan(data%theta)

      call densify(Bperp, x_path, x_path_new)
      call densify(Bpara, x_path, x_path_new)
      call densify(ne, x_path, x_path_new)
      call densify(h, x_path, x_path_new)

      z_path = z_path_new
      x_path = x_path_new
    endif
    ! Adds the z dependence to the density
    ne = ne  * exp(-abs(z_path)/h)
    ! Adds the z dependence the magnetic field
    if (data%B_scale_with_z) then
      ! Scale with z (coordinate)
      Bpara = Bpara * exp(-abs(z_path)/h)
      Bperp = Bperp * exp(-abs(z_path)/h)
    else
      ! Constant for |z|<h, zero otherwise
      where (abs(z_path)/h>1)
        Bpara = 0d0
        Bpara = 0d0
        Bperp = 0d0
      endwhere
    endif

    data%number_of_cells(iz) = data%number_of_cells(iz) + size(x_path)
    if (lI) then
      ! Synchrotron emission
      ! NB Using the total density as a proxy for cosmic ray electron density
      ! NB2 Increments previous calculation
      data%Stokes_I(iz) = data%Stokes_I(iz)+Compute_Stokes_I(Bperp, ne,      &
                                                              x_path, z_path, &
                                                              data%wavelength,&
                                                              data%alpha)
    endif
    if (lRM) then
      ! Faraday rotation (backlit)
      ! NB Using the total density as a proxy for thermal electron density
      data%RM(iz) = Compute_RM(Bpara, ne, x_path, z_path)
    endif

  end subroutine LoSintegrate

  subroutine LoSintegrate_all_redshifts(props, impacty_rmax, impactz_rmax, data, RM_out, I_out)
    type(Galaxy_Properties), intent(in) :: props
    type(LoS_data), intent(inout) :: data
    ! Impact parameter in y- and z-directions in units of rmax
    double precision, intent(in) :: impactz_rmax, impacty_rmax
    logical, intent(in), optional :: RM_out, I_out
    double precision,dimension(props%n_redshifts,2*props%n_grid) :: Bpara_all, &
                                              Bperp_all, xc, zc, ne_all, h_all
    logical, dimension(props%n_redshifts,2*props%n_grid) :: valid
    logical :: lRM, lI
    double precision, allocatable, dimension(:) :: x_path, z_path
    double precision, allocatable, dimension(:) :: Bpara, Bperp, ne, h
    double precision, allocatable, dimension(:) :: z_path_new, x_path_new
    double precision, dimension(props%n_redshifts) :: Bx, By, Bz, rmax
    double precision, dimension(props%n_redshifts) :: tmp, impact_y, impact_z
    double precision :: zmin, zmax, xmin, xmax
    integer, dimension(2*props%n_grid) :: js
    integer :: i, j, it
    logical, parameter :: ldense_z = .true.

    ! Allocates/initialises everything
    valid = .false.
    if (.not.allocated(data%RM)) then
      allocate(data%RM(props%n_redshifts))
      allocate(data%Stokes_I(props%n_redshifts))
      allocate(data%number_of_cells(props%n_redshifts))
      data%RM = 0; data%number_of_cells = 0; data%Stokes_I = 0
    endif
    ! Skips if outside the galaxy
    if (abs(impacty_rmax)>1) return

    lRM = .false.; lI = .true. ! Default values
    if (present(RM_out)) lRM = RM_out
    if (present(I_out)) lI = I_out


    rmax = props%Rcyl(:,props%n_grid)
    impact_y = impacty_rmax * rmax
    impact_z = impactz_rmax * rmax


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
      tmp = props%Rcyl(:,j)**2-impact_y(:)**2 ! NB allocated on-the-fly: Fortran2003 feature
      where (tmp>=0)
        xc(:,i) = xc(:,i) * sqrt(tmp)
        ! Stores a mask to be used with the Fortran 2003 'pack' intrinsic function
        valid(:,i) = .true.
      endwhere

      ! Sets z-coord
      zc(:,i) = xc(:,i) / tan(data%theta) + impact_z(:)

      ! Bx = Br * x/R - Bp * y/R
      Bx = props%Br(:,j) * xc(:,i)/props%Rcyl(:,j) - props%Bp(:,j)*impact_y(:)/props%Rcyl(:,j)
      ! By = Br * y/R + Bp * x/R
      By = props%Br(:,j) * impact_y(:)/props%Rcyl(:,j) + props%Bp(:,j)* xc(:,i)/props%Rcyl(:,j)
      ! Bz = Bzmod * sign(z)
      Bz = sign(props%Bz(:,j), zc(:,i))

      ! Stores density and scale height
      ! NB the scaling with of n with z will be done later!
      ne_all(:,i) = props%n(:,j)
      h_all(:,i) = props%h(:,j)

      ! Simple vector calculations

      ! B_\parallel = dot(B,n), where n is the LoS direction
      ! [NB  n = (sin(theta),0,cos(theta))  ]
      Bpara_all(:,i) = Bx*sin(data%theta) + Bz*cos(data%theta)
      ! B_\perp = \sqrt{ (B_x - B_\parallel\sin\theta)^2 + B_y^2 +
      !                  + (B_z -B_\parallel\cos\theta)^2 }
      Bperp_all(:,i) = ( Bx - Bpara_all(:,i)*sin(data%theta) )**2   &
                         + By**2 +                                  &
                       ( Bz - Bpara_all(:,i)*cos(data%theta) )**2
      Bperp_all(:,i) = sqrt(Bperp_all(:,i))
    enddo

    ! Now work is done for each redshift (as the valid section of each array may
    ! may be different)
    do it=1,props%n_redshifts
      ! Skips empty cells
      if (all(.not.valid(it,:))) cycle

      ! Filters away invalid parts of the arrays
      Bpara = pack(Bpara_all(it,:),valid(it,:))
      Bperp = pack(Bperp_all(it,:),valid(it,:))
      ne = pack(ne_all(it,:),valid(it,:))
      h = pack(h_all(it,:),valid(it,:))
      x_path = pack(xc(it,:),valid(it,:))
      z_path = pack(zc(it,:),valid(it,:))

      ! Skips (bogus) empty cells
!       if (maxval(z_path)-minval(z_path)<1d-6) cycle

      ! The following makes the grid denser, to allow meaningful z-integration
      ! This step is irrelevant for edge-on, but is important for face-on
      if (ldense_z) then
        ! Sets a tentative z-maximum to 3.5 times the scale-height
        zmin = -3.5d0*h(1); zmax = 3.5d0*h(size(h))

        ! Computes corresponding values for x
        xmin = (zmin - impact_z(it))*tan(data%theta)
        xmax = (zmax - impact_z(it))*tan(data%theta)
        ! Forces x coordinate to live within the original grid
        if (abs(xmin) > abs(x_path(1))) then
          xmin = x_path(1)
          zmin = z_path(1)
        endif
        if (abs(xmax) > abs(x_path(size(x_path)))) then
          xmax = x_path(size(x_path))
          zmax = z_path(size(x_path))
        endif

        z_path_new = linspace(zmin,zmax, data%nz_points)
        x_path_new = (z_path_new - impact_z(it))*tan(data%theta)

        call densify(Bperp, x_path, x_path_new)
        call densify(Bpara, x_path, x_path_new)
        call densify(ne, x_path, x_path_new)
        call densify(h, x_path, x_path_new)

        z_path = z_path_new
        x_path = x_path_new
      endif
      ! Adds the z dependence to the density
      ne = ne  * exp(-abs(z_path)/h)
      ! Adds the z dependence the magnetic field
      if (data%B_scale_with_z) then
        ! Scale with z (coordinate)
        Bpara = Bpara * exp(-abs(z_path)/h)
        Bperp = Bperp * exp(-abs(z_path)/h)
      else
        ! Constant for |z|<h, zero otherwise
        where (abs(z_path)/h>1)
          Bpara = 0d0
          Bpara = 0d0
          Bperp = 0d0
        endwhere
      endif

      data%number_of_cells(it) = data%number_of_cells(it) + size(x_path)
      if (lI) then
        ! Synchrotron emission
        ! NB Using the total density as a proxy for cosmic ray electron density
        ! NB2 Increments previous calculation
        data%Stokes_I(it) = data%Stokes_I(it)+Compute_Stokes_I(Bperp, ne,      &
                                                               x_path, z_path, &
                                                               data%wavelength,&
                                                               data%alpha)
      endif
      if (lRM) then
        ! Faraday rotation (backlit)
        ! NB Using the total density as a proxy for thermal electron density
        data%RM(it) = Compute_RM(Bpara, ne, x_path, z_path)
      endif
    enddo
  end subroutine LoSintegrate_all_redshifts

  function IntegrateImage(props,data,iz,number_of_calls,method,error) result(res)
    use fgsl
    use, intrinsic :: iso_c_binding
    type(Galaxy_Properties), intent(in) :: props
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
    real(fgsl_double) :: chisq, xl(2), xu(2), res, err, y, yy, yyy
    double precision, optional, intent(out) :: error
    type(c_ptr) :: ptr
    integer(fgsl_int) :: status
    character(len=10) :: mthd

    calls = 5000
    if (present(number_of_calls)) calls = number_of_calls
    mthd = 'MISER'
    if (present(method)) mthd = method

    ! Global variables for the integration
    data_glb = data
    props_glb = props
    iz_glb = iz

    xl = -1
    xu = 1

    ! Sets up FGSL stuff
    t = fgsl_rng_env_setup()
    r = fgsl_rng_alloc(t)
    gfun = fgsl_monte_function_init(IntegrandImage, 2_fgsl_size_t, ptr)

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

  function IntegrandImage(v_c, n, params) bind(c)
    ! Wrapper to allow using LoSintegrate with FGSL
    use fgsl
    use, intrinsic :: iso_c_binding

    integer(c_size_t), value :: n
    type(c_ptr), value :: v_c, params
    real(c_double) :: IntegrandImage
    real(c_double), dimension(:), pointer :: v
    ! Reads the memory address
    call c_f_pointer(v_c, v, [n])
    data_glb%Stokes_I(iz_glb) = 0d0
    call LoSintegrate(props_glb, v(1), v(2), data_glb, iz_glb, &
                                   RM_out=.false., I_out=.true.)
    IntegrandImage = data_glb%Stokes_I(iz_glb)
  end function IntegrandImage


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

    Compute_Stokes_I = Integrator(integrand, x_path, z_path) * 1d6

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

    emissivity = ncr * Bperp**((alpha+1.0)/2d0) * wavelength**((alpha-1.0)/2d0)

  end function emissivity


  subroutine print_image(props, data, directory, ymax, zmax, nprint, isnap)
      use messages
      character(len=*) :: directory
      character(len=20) :: form
      integer, optional, intent(in) :: nprint, isnap
      integer :: n, i, j, it
      integer, dimension(5) :: unit
      double precision :: impact_y, impact_z, zmax, ymax, rmax
      real, dimension(:,:), allocatable :: RM_im, I_im, N_im
      real, dimension(:),allocatable :: y, z
      type(Galaxy_Properties), intent(in) :: props
      type(LoS_data), intent(inout) :: data

      n = 60; it = 1 ! Default values
      if (present(nprint)) n = nprint
      if (present(isnap)) it = isnap

      rmax = props%Rcyl(it,props%n_grid)

      form ='('//str(n)//'E15.5)'

      open(newunit=unit(1),file=directory//'I.dat', FORM='FORMATTED', status='replace')
      open(newunit=unit(2),file=directory//'RM.dat', FORM='FORMATTED',status='replace')
      open(newunit=unit(3),file=directory//'cells.dat', FORM='FORMATTED',status='replace')
      open(newunit=unit(4),file=directory//'y.dat', FORM='FORMATTED', status='replace')
      open(newunit=unit(5),file=directory//'z.dat', FORM='FORMATTED',status='replace')

      allocate(I_im(n,n))
      allocate(RM_im(n,n))
      allocate(N_im(n,n))
      allocate(y(n))
      allocate(z(n))



      do j=1,n
        impact_z = -zmax + 2*zmax/dble(n)*j
        z = impact_z
        do i=1,n
          impact_y =  -ymax + 2*ymax/dble(n)*i
          y(i) = impact_y
!           call LoSintegrate_all_redshifts(props, impact_y/rmax, impact_z/rmax, data, &
!                             I_out=.true., RM_out=.true.)
          call LoSintegrate(props, impact_y/rmax, impact_z/rmax, data, it, &
                            I_out=.true., RM_out=.true.)
          I_im(i,j) = data%Stokes_I(it)
          data%Stokes_I(it) = 0
          N_im(i,j) = data%number_of_cells(it)
          data%number_of_cells(it) = 0
          RM_im(i,j) = data%RM(it)
          data%RM(it) = 0
        enddo
        write(unit(1), trim(form)) I_im(:,j)
        write(unit(2), trim(form)) RM_im(:,j)
        write(unit(3), trim(form)) N_im(:,j)
        write(unit(4), trim(form)) y
        write(unit(5), trim(form)) z
      enddo

      do i=1,5
        close(unit(i))
      enddo
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
    allocate( galprops%Br(nz,nr),   &
              galprops%Bp(nz,nr),   &
              galprops%Bz(nz,nr),   &
              galprops%Rcyl(nz,nr), &
              galprops%h(nz,nr)     )
    galprops%n_redshifts = nz
    galprops%n_grid = nr
  end subroutine alloc_Galaxy_Properties

  function linspace(z_min, z_max, n) result(array)
    double precision, intent(in) :: z_min,z_max
    double precision, allocatable, dimension(:) :: array
    double precision :: step
    integer :: n,i
    step = (z_max-z_min)/dble(n-1)
    allocate(array(n))
    array = (/( (i-1)*step+z_min, i=1 , n)/)
  end function linspace

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

end module LoSintegrate_aux
