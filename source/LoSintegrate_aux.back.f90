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
  public :: write_I_textfile, init_write_I_textfile


  ! 30% fixed ionisation fraction!
  double precision, parameter :: ionisation_fraction = 0.3
  type Galaxy_Properties
    double precision, allocatable, dimension(:,:) :: Br, Bp, Bz
    double precision, allocatable, dimension(:,:) :: Rcyl, h, n
    double precision, allocatable, dimension(:) :: z, Mgas_disk, Mstars_disk, r_disk

!for small scale field    
    double precision, allocatable, dimension(:,:) :: beq, shear_galaxy
    double precision, allocatable, dimension(:) :: tau, Uz, l0, v0 
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
    integer :: inc_random_field = 1 !0 - no small scale field. 1-isotropic small scale field. 2-anisotropic small scale field 
    integer :: nz_points = 500
  end type

  !Global and public variables
  double precision, public,  dimension(300) :: write_IQU=0.0
  integer, public :: iw, iw_gal 
  
   !Global variables for Small scale field
  integer, parameter:: steps_rint=200, steps_pint=200
  double precision, dimension(steps_rint)::br2_muG2,sheargal, Sfac, Sfacsq 
  double precision :: tau0_myr, Uz_kms, Lu_kpc,const_factor_b 
  double precision :: bp2mid_muG2, bz2mid_muG2, br2mid_muG2, brbpmid_muG2, bpbzmid_muG2, brbzmid_muG2
  double precision :: Sfac_mid, Sfacsq_mid
  double precision :: bxp2_muG2, byp2_muG2
  double precision :: B2Lavg, b2savg
  double precision, parameter :: v0_ISM_cms = 10d5  !rms speed of turbulant ism 10 km/s = 10d5 cm/s 

  ! Global variables for the integration
  type(LoS_data) :: data_glb
  logical :: dust_glb
  type(Galaxy_Properties) :: props_glb
  integer :: iz_glb
  integer :: print_cj = 1, test_CJ=0 
  logical :: CR_B_equipartition = .true.
  contains

  subroutine LoSintegrate(props, impacty_rmax, impactz_rmax, data, iz, &
                          FRB_mode, RM_out, I_out, Q_out, U_out, iRM, error, &
                          DM_out, SM_out, dust_mode, test)
    use input_constants
    use messages
    use global_input_parameters, only: p_molecularHeightToRadiusScale,  p_ISM_turbulent_length
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
    double precision :: rmax, lambda_m, h_m
    double precision :: tmp, impact_y, impact_z
    double precision :: z_start, z_end 
    integer :: i, j, iz, this_RM
    logical, parameter :: ldense_z = .true.
    logical :: lFRB, lTest, ldust
    logical, optional, intent(inout) :: error
    double precision, parameter :: EPS=1e-6
    logical, intent(in), optional :: dust_mode

    integer :: rint 
    double precision, dimension(props%n_grid) :: rval,Bzrval
    real(fgsl_double) :: rlim(2)
    double precision:: drval


    lFRB = .false.
    if (present(FRB_mode)) lFRB = FRB_mode

    ! Allocates/initialises everything
    valid = .false.
    if (.not.allocated(data%RM)) then
              !props%n_RMs = 1 for observables_single, else it can be more than 1
      allocate(data%RM(props%n_RMs))
      allocate(data%column_density(props%n_RMs))
      data%RM = 0; data%number_of_cells = 0;
      data%Stokes_I = 0; data%Stokes_U = 0; data%Stokes_Q = 0
    endif
    ! Skips if outside the galaxy
    
    if (abs(impacty_rmax)>1) then
      if (print_cj ==1) then
         print *, 'mod(impacty_rmax) greater than 1', impacty_rmax
         print *, 'not computing RM'
      end if 
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

    if (test_CJ==1) then
      rlim(1) = 0.0001
      rlim(2) = 1.0
      rmax = rlim(2)
      drval  = ( rlim(2) - rlim(1) )/props%n_grid 
!creating  grids rval and pval with all r and phi points used in integration
      rval(1) =  rlim(1) +  (drval/2.0)
      do rint =2, props%n_grid, 1
         rval(rint) = rval(rint-1)  + drval
     end do
    impact_y = impacty_rmax * rmax
    impact_z = impactz_rmax * rmax
    end if 


    ! Computes wavelength at the frame of the galaxy
    lambda_m = data%wavelength/(1d0 + props%z(iz) )
    if (test_CJ == 1) then
       lambda_m = 1.0 
    end if 


    ! Loops over the radial grid-points
    ! i -> index for quantities in a cartesian box
    ! j -> index for axi-symmetric quantities (for radial grid points)
    do i=1,2*props%n_grid
      ! Constructs coordinates and auxiliary indices js
      if (i<props%n_grid+1) then
        ! Index for behind the x-z plane
        j=props%n_grid-i+1
        !js(i) = j
        ! sign for x coordinate
        xc(i) = -1
      else
        ! Index ahead of the x-z plane
        j = i-props%n_grid
        !js(i) = j
        ! sign for x coordinate
        xc(i) = 1
      end if

      ! The available values for x are x_i = sqrt(R_i^2-b^2)
      tmp = props%Rcyl(iz,j)**2-impact_y**2

      if (test_CJ==1) then
        tmp = rval(j)**2-impact_y**2
      end if

      if (tmp>=0) then
        xc(i) = xc(i) * sqrt(tmp)
        ! Stores a mask to be used with the Fortran 2003 'pack' intrinsic function
        valid(i) = .true.
      endif

      
      ! Sets z-coord
      zc(i)     = xc(i) / tan(data%theta) + impact_z
   
      ! Bx = Br * x/R - Bp * y/R
      Bx_all(i) =   props%Br(iz,j) * xc(i)/props%Rcyl(iz,j)    &
                  - props%Bp(iz,j)*impact_y/props%Rcyl(iz,j)
      ! By = Br * y/R + Bp * x/R
      By_all(i) = props%Br(iz,j) * impact_y/props%Rcyl(iz,j)   &
                  + props%Bp(iz,j)* xc(i)/props%Rcyl(iz,j)
      ! Bz = Bzmod * sign(z)
      Bz_all(i) = sign(props%Bz(iz,j), zc(i))
      !print *, zc(i), Bz_all(i)
      ! Stores density and scale height
      ! NB the scaling with of n with z is done further ahead!
      ne_all(i) =  props%n(iz,j) * ionisation_fraction
      h_all(i)  = props%h(iz,j)
     

    if (test_CJ == 1) then
      !beatit 
      Bx_all(i) =  0.0 * xc(i)/rval(j) - 0.0 * impact_y/rval(j)
      By_all(i) =  0.0 * impact_y/rval(j) + 0.0 * xc(i)/rval(j)
      ! Bz = Bzmod * sign(z)
      Bzrval(j) =  1.0 
      Bz_all(i) = sign(Bzrval(j), zc(i))
      ! NB the scaling with of n with z is done further ahead!
      ne_all(i) = 1.0 
 
      h_all(i) = 0.05 
   end if
 
     if (ldust) then
        ! At the moment, the following is NOT strictly speaking the molecular
        ! density, but is proportional to it
        ! NB the scaling with of n with z is done further ahead!
        ! r_disk, Mgas_disk, Mstars_disk etc has a single value
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
    ! ldense_z is always assigned to be .true. 
    if (ldense_z) then

      if (z_path(1) < z_path(size(z_path))) then
        ! Sets a tentative z-minimum/maximum to -3.5/3.5 times the (maximum) scale-height
        z_start = -3.5d0*h(1)
        z_end = 3.5d0*h(size(h))

        ! If in dust mode, use the molecular gas as reference instead
        if (ldust) then
          h_m = 0
          print *, 'initialize h_m properly. Currently h_m is set to', h_m
          z_start = max( z_start, -3d0*h_m )
          z_end = min( z_end, 3d0*h_m )
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
        z_end   = -3.5d0*h(size(h))

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
    if (test_CJ ==1) then
      ne = ne  * exp(abs(z_path)/h)
    end if
    ! And molecular density
    h_m = props%r_disk(iz) * constDiskScaleToHalfMassRatio * p_molecularHeightToRadiusScale
    n_mol = n_mol * exp( -abs(z_path)/h_m )
    if (test_CJ ==1) then
      n_mol = n_mol * exp( abs(z_path)/h_m )
    end if
    
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
    Brnd = 0 * sqrt(4*pi*ne*Hmass)*10d5 * 0.5d0 *1d6

    Brnd_B =  Brnd/sqrt(Bperp**2 + Bpara**2 + 1e-40)
    

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
                                     lambda_m, data%alpha, Brnd)
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
                + (lambda_m)**2 * Compute_RM(Bpara(1:i),ne(1:i), &
                                                   x_path(1:i),z_path(1:i))
      enddo

      if (lQ) then
        data%Stokes_Q = Compute_Stokes('Q', Bperp, ne, x_path, z_path,&
                                       lambda_m, data%alpha, psi)
      endif
      if (lU) then
        data%Stokes_U = Compute_Stokes('U', Bperp, ne, x_path, z_path,&
                                       lambda_m, data%alpha, psi)
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


!Function that Observables/Observables_single calls to calculate various intensities, RM etc   
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
   real(fgsl_double) :: xl(2), xu(2), integral, res, err
   double precision, optional, intent(out) :: error
   type(c_ptr) :: ptr
   integer(fgsl_int) :: status
   character(len=10) :: mthd

! Below varibles are used to find the integrals in brute force manner with method=IntLoS
   integer:: xint, zint, steps_xint, steps_zint
   double precision :: theta_rad, integrant_syn, xval, zval, dxval, dzval
   real :: time1, time2, time3
! These stores total and polarized intensities
   double precision:: res1, res2, res3, resq, resu
!steps sizes and grid sizes for los integration. 
   double precision:: d_los,hfac
   integer :: steps_xp, steps_yp, steps_los_orig, i_run_vol_int 
   logical ::run_vol_int 
   double precision :: resq_last, resu_last

!Spacify this to calculate total synchrotron intensity as 1D/2D integral   
   character(len=2) :: vol_int_dimn
!Below variables stores time taken for calculation as well as accuracy of result.
!If error_>= p, then program to calculate intensities are rerun with 
!finer integration grid to improve accuracy
   double precision :: time_calc, error_i, error_q, error_u
!tolerence for stokes parameters in percent. Also threshold (I_pol/I).  
   double precision :: tol_i, tol_qu, pol_th



    dust_glb = .false.
    if (present(dust)) dust_glb = dust

    calls = 2500
    if (present(number_of_calls)) calls = number_of_calls
    mthd = 'MISER'
    if (present(method)) mthd = method

    data_glb = data
    props_glb = props
    iz_glb = iz

    calls = 3700
    theta_rad = data_glb%theta

    xl(1) = -1
    xu(1) = 1
!    xl(2) = -1/ data%theta
!    xu(2) = 1/ data%theta

!    mthd = 'IntLoS' ! Not an LoS integral, but a volume integral
!    mthd = 'IntVol'
!    mthd = 'plain' ! 
    mthd = 'IntLoS' ! 
     

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
      stop 'Error'
    endif

   select case (trim(mthd))
      case('MISER')
        m = fgsl_monte_miser_alloc(2_fgsl_size_t)
    
        status = fgsl_monte_miser_integrate(gfun, xl, xu, 2_fgsl_size_t, &
                                            calls, r, m, res, err)
        call fgsl_monte_miser_free(m)
      case('VEGAS')
        v = fgsl_monte_vegas_alloc(2_fgsl_size_t)
        ! This function runs 'calls' number of times
        status = fgsl_monte_vegas_integrate(gfun, xl, xu, 2_fgsl_size_t, &
                                            calls, r, v, res, err)
        call fgsl_monte_vegas_free(v)
      case('plain')
        s = fgsl_monte_plain_alloc(2_fgsl_size_t)
        status = fgsl_monte_plain_integrate(gfun, xl, xu, 2_fgsl_size_t, &
                                            calls, r, s, res, err)
        call fgsl_monte_plain_free(s)

      !To Compute Total Intensity by 2D and ID integrals using volume integral method (CJ)
      case('IntVol') 
          vol_int_dimn = '1D'
          call Total_Synchrotron_Intensity_volume_integral(props_glb, data_glb, iz_glb, vol_int_dimn, res1)
          print "(A,  e12.4, 3x, A, I2)", 'Tot Intensity=', res1, 'Dimn of Volume integral=', vol_int_dimn  
          res = res1
      !To Compute total and polarized Intensities by doing an LoS integral along different LoSs (CJ)
      case('IntLoS') 
        CR_B_equipartition = .true. 
        ! Calls function to compute the synchrotron total intensity by doing volume integral to obtain I_vol (stored in res2) 
        ! I_vol will be compared with total synchrotron intensity, I_los, obtained through the LoS integration (stored in res3).
        ! We adjust grid points used for integration until I_los is within 3-5% of I_vol.
        vol_int_dimn = '2D'
        call Total_Synchrotron_Intensity_volume_integral(props_glb, data_glb, iz_glb, vol_int_dimn, res2)
        ! Calls functions to compute total and polarized intensities (stokes parameters) by doiing an LoS integral       
        ! The accuracy of calculation depends on  steps_xp^2, the number of grid points on x'y plane (see manual). 
        ! steps_xp is varied in steps until integrations to compute intensities converges or integral falls below some threshold. 
        ! The program will run atleast for 2 values of steps_xp because only then we can determine if Q and U are converged. 

        
        !Parameters affecting convergence of results including steps_xp
        steps_xp       = 50 !Should be even
        steps_los_orig = 70  !if this is bigger than 90, numerical issues
        d_los          = 2.e-2 !Minimum grid seperation along any given LoS 
        hfac           = 0.8   !Maximum Height from disc plane upto which integrations are done is hfac * max(scale height)
                              !If hfac = 1 produces good accuracy for exponentially decaying magnetic fields and gas densities
      
!tol_i and tol_qu -  percentage numerical convergence of integral for computing I,Q and U should be less than this 
!pol_th - If percentage of total Q or U is less than pol_th * I, then it won't computed.  
        tol_i  = 5.0
        tol_qu = 5.0
        pol_th = 0.01        

! resq_last, resu_last store  Q,U computed in previous integrations which will be compare with values computed 
! from latest integration (stored in resq, resu) until required tolerence is reached. 
! res3 is the total I comuted from LoS method and compared with res2 (computed above) until tolerence is reached. 
        resq_last = 0.0    
        resu_last = 0.0
        res3  = 0.0
        resq  = 0.0
        resu  = 0.0

!To avoid numerical convergence issues, intensities are computed multiple times by varying steps_xp until convergence is reached        
        run_vol_int   = .true. !Intensities are computed for different steps_xp as long as run_vol_int = true
        i_run_vol_int = 0
        if (res2 < 1.e-27) then !If I is below a threshold, there is no point in finding I,Q,U using LoS integral 
          run_vol_int= .false.
          print *, 'I(syn) less than threshold. Will not compute I, Q and U', props%igal, iz, res2
        end if
        time_calc = 0.0 !total time for obtaining results
        do while (run_vol_int) 
          steps_yp =  steps_xp 
     
          call cpu_time(time1)
                    !Stokes parameters I, Q, U are computed in one go along an LoS and stored in res3, resq and resu
          call Total_Synchrotron_Intensity_LoS_integral(props_glb, data_glb, iz_glb, steps_xp, steps_yp, & 
                                                           steps_los_orig, d_los, hfac,  res3, resq, resu)
          time_calc = time_calc + (time2-time1)

!To check accuracy of computed I, Q, U. The integration stops when specified accuray is reached or intensities fall below 
!a thresold or the specified grid size steps_xp is greater than some thresold.  
!When calculated intensity is zero, then also error is set to zero. 
          if (res2 .ne. 0.0) then 
            error_i= 100*abs( (res3-res2)/res2 )
          else 
            error_i = 0.0
          end if
          if (resq .ne. 0.0) then 
            error_q = 100*abs( (resq-resq_last)/resq)
          else 
            error_q = 0.0
          end if
          if (resu .ne. 0.0) then 
            error_u = 100*abs( (resu-resu_last)/resu)
          else 
            error_u = 0.0
          end if
          !To compare previous Q and U values with new Q and U values 
          resu_last = resu
          resq_last = resq
      
          ! Below decides when to stop the program
          ! If errors are below a threshold (or zero), then run_vol_int is set to False.
          ! since i_run_vol_int >=1, program is run atleast twice 
          if (i_run_vol_int >=1 .and. res3 > 0.0) then
            ! If all integrals are converged, then run_vol_int is set to False
            ! This also happens if all integrals give zero result
            if (error_i<= tol_i .and. error_u <= tol_qu .and. error_q <= tol_qu ) then
              run_vol_int= .False.
            end if 
            !If Q is below a fraction of I, then run_vol_int is set to False, provided I, U integrals converged.
            if (error_i<= tol_i .and. error_u<= tol_qu .and. resq/res3 <= pol_th ) then 
               run_vol_int= .False.
               resq = 0.0
            end if   
            !If U is below a fraction of I, then run_vol_int is set to False, provided I, Q integrals converged.
            if (error_i<= tol_i .and. error_q <= tol_qu .and. resu/res3 <= pol_th ) then 
               run_vol_int= .False.
               resu = 0.0
            end if  
            !If U and Q are below a fraction of I, then run_vol_int is set to False
            !computed atleast twice. 
            if ( error_i<= tol_i .and. resq/res3 <= pol_th .and. resu/res3 <= pol_th ) then 
               run_vol_int= .False.
               resu = 0.0
               resq = 0.0
            end if  
          end if  
        
        ! Priting all outputs .  
          print "(i5, 2f10.4, 3i5, 3f5.2, 4e12.4, 4f10.3)",props%igal,props%z(iz), data%theta, steps_xp, steps_yp, &
               steps_los_orig,  data%wavelength, d_los, hfac, res2, res3, resq, resu, error_i, error_q, error_u, time_calc

          open(33,file='./output_text_files/out_I.txt', FORM='FORMATTED', status='old', position="append")
             write(33,"(i5, 2f12.6, 3i5, 3f5.2, 4e13.5, 4f12.3)" ),props%igal,props%z(iz), data%theta, steps_xp, steps_yp, &
              steps_los_orig,  data%wavelength, d_los, hfac, res2, res3, resq, resu, error_i, error_q, error_u, time_calc 
          close(33)

          !steps_xp is varied from 50-250 in steps of 50
          steps_xp = steps_xp +50
          if ( steps_xp >= 270) then 
            run_vol_int= .False.
          end if    
          i_run_vol_int = i_run_vol_int + 1 
        end do
        ! Writing final (possibly) converged output to a file. If steps_xp >250, not sure if converged unless q and u are zero.  
        open(34,file='./output_text_files/out_converged_I.txt', FORM='FORMATTED', status='old', position="append")
          write(34,"(i5, 2f12.6, 3i5, 3f5.2, 4e13.5, 4f12.3)"),props%igal,props%z(iz), data%theta, steps_xp, steps_yp, &
           steps_los_orig,  data%wavelength, d_los, hfac, res2, res3, resq, resu, error_i, error_q, error_u, time_calc 
        close(34)
        res = res3
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
    double precision, dimension(size(ne)) :: integrand,  l
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
    double precision :: Compute_Test
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


  ! This routine computes the total synchrotron intensity (not polarized intensity) by doing a volume integral. 
  ! The z-integral can be done analytically if we assume a functional form for the z-dependence of magnetic fields 
  ! and gas density. This will reduce the 3D integral to a 2D integral. 
  ! Also in there is cylindrical symmetry to the problem, the integration over azimuthal angle can be also done 
  ! analytically. This further reduces the original 3D integral to  1D numerical integral. 

subroutine Total_Synchrotron_Intensity_volume_integral(props,  data, iz, vol_int_dimn,  res1)
    use input_constants
    use messages
    use global_input_parameters, only: p_molecularHeightToRadiusScale, p_ISM_turbulent_length
    use input_constants
    use FRB
    use interpolation
    type(Galaxy_Properties), intent(in) :: props
    type(LoS_data), intent(inout) :: data

    integer ::  j, iz, this_RM
    double precision, intent(inout) :: res1 

    integer :: xint, zint, steps_xint, steps_zint
    double precision ::   zval
    real(fgsl_double) :: rlim(2), plim(2)
    integer:: rint, pint
    double precision :: res, lambda_m
   
    character(len=2), intent(in) :: vol_int_dimn
    integer, parameter:: steps_rint=150, steps_pint=150
    double precision, dimension(steps_rint) :: rval, B2_rval_muG2, Br_muG,  Bp_muG, Bz_muG,  &
           n_cm3, ncr_cm3, h_kpc
    double precision, dimension(steps_pint) :: nr_LoS, np_LoS, pval 
    double precision, dimension(steps_rint, steps_pint) :: Bparal_up, Bparal_low, int_syn_grid  
    double precision, dimension(steps_rint, steps_pint) :: int_syn_full_grid 

    double precision :: nz_LoS, Bmag, Bmag2, drval, dpval
    double precision :: sint, cost, splus1by4, splus1by2, sminus1by2, s, factor_z, factor_r, int_syn    
    !temporary variable
    double precision :: Brmid, Bpmid, Bzmid, ncrmid, Ke_fac, hmid, int_syn1, int_syn2, rmid
    real time1, time2, time3
    double precision :: Bmax, Bavg, B2avg, havg

! Volume integral in (r, phi,z) to calculate the total synchrtron emission.  
! result = integral_V  epsilon;  episilon is the emissivity at a given r, phi, z
! The z-integral is always done analytically and r-integrals is always done numerically
! It  phi integral also numerically and analytically depending on model
! limits of r and phi variables 

    rlim(1) = props%Rcyl(iz,1)+1.e-10    !limits of r-variable
    rlim(2) = props%Rcyl(iz,props%n_grid)-1.e-10 
    plim(1) = 0.0 !limits of phi-variable 0<phi<2pi
    plim(2) = 6.28318530718
    drval  = ( rlim(2) - rlim(1) )/(steps_rint-1) !dr 
    dpval  = ( plim(2) - plim(1) )/(steps_pint-1)   !dphi

!creating  grids rval and pval with all r and phi points used in integration
    do rint =1, steps_rint, 1
         rval(rint) = rlim(1) + (rint-1) * drval
    end do
    do pint =1, steps_pint, 1
       pval(pint)  =  plim(1) + (pint-1) * dpval
    end do 

  !Below step is to avoid floating point errors. Occassionally you may get  ray_kpc(steps_rint)> rlim(2) 
    if (rval(steps_rint)> rlim(2) ) then
          rval(steps_rint) = rlim(2)
    end if
    if (pval(steps_pint)> plim(2) ) then
          pval(steps_pint) = plim(2)
    end if

! Wavelength of observation 
    lambda_m = data%wavelength/(1d0 + props%z(iz))

!Input parameters 
    call interpolate( props%Rcyl(iz,:), props%Br(iz,:), rval, Br_muG)
    call interpolate( props%Rcyl(iz,:), props%Bp(iz,:), rval, Bp_muG)
    call interpolate( props%Rcyl(iz,:), props%Bz(iz,:), rval, Bz_muG)
    call interpolate( props%Rcyl(iz,:), props%n(iz,:), rval, n_cm3)
    call interpolate( props%Rcyl(iz,:), props%h(iz,:), rval, h_kpc)
    ncr_cm3  = n_cm3  * ionisation_fraction


    !To find the maximum and average value of magnetic field
    B2_rval_muG2 = Br_muG**2.0 + Bp_muG**2.0 + Bz_muG**2.0
    Bmax  = sqrt(  MAXVAL(B2_rval_muG2) )
    Bavg  = 0.0
    B2avg = 0.0
    havg  = 0.0
    do rint = 1, steps_rint-1, 1
         Bavg  = Bavg  +   sqrt(B2_rval_muG2(rint)) * rval(rint) * ( rval(rint+1)-rval(rint) ) 
         B2avg = B2avg +   B2_rval_muG2(rint) * rval(rint) * ( rval(rint+1)-rval(rint) ) 
         havg  = havg  + h_kpc(rint) * rval(rint) * ( rval(rint+1)-rval(rint) )
    end do
    Bavg  = Bavg  / ( rlim(2)**2.0/2.0 ) !dividing by area = integral r dr = r^2/2
    B2avg = B2avg / ( rlim(2)**2.0/2.0 ) !dividing by area = integral r dr = r^2/2
    havg  = havg /  ( rlim(2)**2.0/2.0 )

    print "(A5,  e12.6, A5, e12.6, A8, e12.6, A8, e12.6, A8, e12.6, A8, e12.6)",  &
             'z=', props%z(iz), 'r=', rlim(2), 'Bmax=', Bmax, 'Bavg=', Bavg, 'B2avg=', B2avg, 'havg=', havg   

!factor to account for integration along z-direction
    factor_z = 1.0/3.0 !.  !h(r)*factor_z = int_0^infinity exp(-3z/h(r)) see below  
                           !
    if (test_CJ == 1) then ! test run where simple models are assumed 
       call Test_Model(rval, Br_muG, Bp_muG, Bz_muG, ncr_cm3, h_kpc, drval, &
            rlim,plim, factor_z, lambda_m, steps_rint, steps_pint, data%theta) 
    end if    

!Cosine of Sine of LoS angle
    cost = COS(data%theta) 
    sint = SIN(data%theta)

   s          = data%alpha
   splus1by4  = (s+1.0)/4.0
   splus1by2  = 2.0*splus1by4  ! power of B_perp
   sminus1by2 = (s - 1.0)/2.0

   !Computing intensity by doing a single numerical r-integral
   if (vol_int_dimn == '1D') then
     print *, 'cosmic ray - magnetic field equipartition can not be applied for 1D integration'
     int_syn1 = 0
     int_syn2 = 0
     do rint =1, steps_rint-1, 1
       Brmid  =     ( Br_muG(rint+1)+Br_muG(rint) )/2.0 
       Bpmid  =     ( Bp_muG(rint+1)+Bp_muG(rint) )/2.0  
       Bzmid  =      ( Bz_muG(rint+1)+Bz_muG(rint) )/2.0 
       ncrmid =     ( ncr_cm3(rint+1)+ncr_cm3(rint) )/2.0 
       hmid   =     ( h_kpc(rint+1)+h_kpc(rint) )/2.0 
       rmid   =     ( rval(rint+1)+rval(rint) )/2.0 

       int_syn1 = int_syn1 + rmid * ncrmid * hmid * (Brmid**2.0+Bpmid**2.0) * ( rval(rint+1) - rval(rint) )
       int_syn2 = int_syn2 + rmid * ncrmid * hmid *  Bzmid**2.0 * ( rval(rint+1) - rval(rint) )
     end do
     int_syn = (6.2830*factor_z) * ( (2.0-sint*sint)*int_syn1 + 2.0* (1.0-cost*cost)*int_syn2 )!6.2830 accounts for phi integral
     res1  = int_syn*lambda_m** sminus1by2 

 
    ! 2D numerical integration on the galactic disc plane over  r and phi grid points defined above.
    ! seperate terms to integrate above and below the galactic plane
   else if ( vol_int_dimn == '2D' ) then
         !Here if equiparition b/w B and Cosmic rays assumed then integrant scale as exp(-4z/h(r)).  
      if ( CR_B_equipartition ) then 
         factor_z = 1.0/4.0 
      end if

         ! Array to store Components of LoS unit vectors as a function of phi 
      do pint =1, steps_pint, 1
        nr_LoS(pint) = COS( pval(pint) ) 
        np_LoS(pint) = SIN( pval(pint) ) 
      end do
         !LoS unit vectors
      nr_LoS = sint * nr_LoS 
      np_LoS = -sint * np_LoS 
      nz_LoS = cost

      do rint =1, steps_rint-1, 1
                   !Evaluating functions at grid mid points
        Brmid =    ( Br_muG(rint+1)+Br_muG(rint) )/2.0 
        Bpmid =    ( Bp_muG(rint+1)+Bp_muG(rint) )/2.0  
        Bzmid =    ( Bz_muG(rint+1)+Bz_muG(rint) )/2.0 
        hmid  =    ( h_kpc(rint+1)+h_kpc(rint) )/2.0 
        rmid  =    ( rval(rint+1)+rval(rint) )/2.0 
        ncrmid  =  ( ncr_cm3(rint+1)+ncr_cm3(rint) )/2.0 
        Bmag2  =   Brmid**2.0 + Bpmid**2.0 + Bzmid**2.0 !B^2

        if ( CR_B_equipartition ) then 
           Ke_fac = Ke_equipart( Bmag2, s)  !ncrmid is essentially the factor KE
        else
           Ke_fac = ncrmid 
        end if

        factor_r =  Ke_fac * rmid * (hmid * factor_z)  !hmid*factor_z = int_0^infinity exp(-3z/h(r))  
        do pint =1, steps_pint-1, 1
           Bparal_up(rint, pint)    = Brmid*nr_LoS(pint) + Bpmid*np_LoS(pint) + Bzmid*nz_LoS ! B.n 
           Bparal_low(rint, pint)   = Brmid*nr_LoS(pint) + Bpmid*np_LoS(pint) - Bzmid*nz_LoS ! B.n 

           int_syn_grid(rint, pint) = ( 2.0*Bmag2 -  Bparal_up(rint, pint)**splus1by2 - & 
                                       Bparal_low(rint, pint)**splus1by2 ) *factor_r 
           int_syn_grid(rint, pint) = int_syn_grid(rint, pint) * ( pval(pint+1) -pval(pint) )
        end do
      int_syn_grid(rint, :) = int_syn_grid(rint, :) * ( rval(rint+1) -rval(rint) )
      end do
 
      int_syn= 0.0 !to sum over integrant for all (r,phi) grid points
      do rint =1, steps_rint-1, 1
        do pint =1, steps_pint-1, 1
          int_syn = int_syn +  int_syn_grid(rint, pint)
        end do
      end do

      res1 = int_syn*lambda_m**sminus1by2
      if ( CR_B_equipartition ) then 
        res1 =  ext_factor_luminosity(props%z(iz))*emissivity_syn_constant(s)*res1 !to multiply various res1 with numerical constants
      else 
        res1 =  ext_factor_luminosity(props%z(iz))*res1 !to multiply various res1 with numerical constants
      end if
  !3D numerical volumen integral
   else if  ( vol_int_dimn == '3D' ) then 
       print *, '3D volume integral not implemented yet'
       stop
   else 
       print *,  'Set vol_int_dimn = 1D, 2D or 3D'
       stop
   end if
end subroutine Total_Synchrotron_Intensity_volume_integral

! This is a test model that can be used to test the numerical computation of total and polarized intensities. 
! The test model assumes a perfectly cylindrical galaxy. The magnetic fields and gas density are z-independent 
! inside the cylinder, but zero outside it. 
! The mangetic field and gas density is also a constant inside the disc. The intensity of the galaxy in this case 
! can be computed analytically and compared with numerical results. 

subroutine Test_Model(rval, Br_muG, Bp_muG, Bz_muG, ncr_cm3, h_kpc, drval, &
                rlim,plim,  factor_z, lambda_m, steps_rint, steps_pint, theta)
 
    real(fgsl_double), intent(inout) :: rlim(2), plim(2)
    integer, intent(in):: steps_rint, steps_pint
    double precision, intent(inout), dimension(steps_rint) :: rval, Br_muG,  Bp_muG
    double precision, intent(inout), dimension(steps_rint) ::  Bz_muG, ncr_cm3, h_kpc
    double precision, intent(inout) :: factor_z, lambda_m, theta, drval
    integer :: rint

    print *, 'test_mode'  
            ! r changes from 0 to 1. Disc radius 1
    rlim(1) = 0.0
    rlim(2) = 1.0
    drval  = ( rlim(2) - rlim(1) )/(steps_rint-1) !dr 
    do rint =1, steps_rint, 1
         rval(rint) = rlim(1)  + (rint-1)*drval
            !magnetic field, density, scale height etc. assigned
         Br_muG(rint)  = 0.0
         Bp_muG(rint)  = 0.0
         Bz_muG(rint)  = 1.0 
         ncr_cm3(rint) = 1.0 
         h_kpc(rint)   = 0.05  !height of disc
    end do

  
    factor_z = 1.0 !in test run, z-dependence is absent. So factor_z=1.  
    lambda_m = 1.0
    theta = 2.0*ATAN(1.d0) 

end subroutine Test_Model

!Several constants appear in the calculation for total synchrotron intensity. They  are estimated here. 
!This will convert synchrotron intensity from arbitrary units to Physical units
function emissivity_syn_constant(s) result(const)
 double precision, intent(in) :: s
 double precision :: const 
 double precision :: c, m_e, e, kpc 
  c     =  2.99792458e10 !cm s^-1 units
  m_e   =  9.1094e-28 ! weight in gram. cgs units
  e     =  4.8032e-10 ! cm^3/2 g^1/2 s^-1 units
  kpc    =  3.086e+21 ! kparsec in cm
  const =  e**3.0/(m_e*c**2.0) * ( 3.0*e/( 12.5663706144 * m_e**3.0 * c**5.0) )**((s-1.0)/2.0) &
           * sqrt(3.0) / ( 12.5663706144 * (s+1.0) ) * gamma( (3.0*s-1.0)/12.0 ) * gamma( (3.0*s+19.0)/12.0 ) & 
           / (c/100.00) ** ((s-1.0)/2.0) & !This factor is included because emissivity is calculated using lambda and not nu. 
                                           !Also lambda is in meter whereas c is here in cm/s. So c/100 -> c in meter/s
           / 1.e-23   & !converting to Janski
           * kpc**3.0 & ! volume element in kpc**3.0
           * (1.e-6)**((s+1.0)/2.0) !magnetic field in units of muG. Converting it into Gauss
end function  emissivity_syn_constant    

!Determines  K_E  (related to cosmic ray particle number density) by assuming equipartition between local magnetic fields and cosmic rays
function Ke_equipart( Bsq,s) result(Ke)
 double precision, intent(in) :: Bsq, s
 double precision :: Ke
 double precision :: ratio_np, E1, E2
    ratio_np   = 0.01
    E1         = 5
    E2         = 100
    Ke   =  ratio_np * (s-2.0) * (1.e-12 * Bsq/25.1327412287)/( ((0.00160218)**(2-s)) * (E1**(2.0-s) - E2**(2.0-s)) )!25.1327412287 = 8 pi
                                                  !1.0e-12 -> Magnetic field from micro-gauss * Gauss 
                                                  ! (0.00160218)**(2.0-s) -> Energy in GeV to erg
                                                  ! Ke/Bsq ~ 3.35*10**-18
end function Ke_equipart        

! For computing the relation between emitted luminosity and observed luminosity, luminosity distance is calculated below
! Also accounts for the k-correction to account for the shrinking of frequency band due to expansion of the Universe
function ext_factor_luminosity(z1) result(ext_factor)
 double precision, intent(in) ::  z1
 double precision :: ext_factor 
 double precision :: var1, var2, z2
 double precision :: omega_lambda, omega_m, h, Mpc, zfinal
 omega_lambda = 0.7 
 omega_m      = 1-omega_lambda
 h            = 0.7
 Mpc    =  3.086e+24 ! megaparsec in cm

 if (z1 <= 0.001) then
  zfinal = 0.001
 else 
  zfinal = z1
 end if

 var1 = 0.0 ;
 z2   = 0.0
 do while ( z2 <= zfinal )
   z2   = z2 + 0.0002 ;
   var2 = 1.0 / sqrt( omega_lambda + omega_m * (1.0+z2)**3.0 ) 
   var1 = var1 + 0.0004 * var2 
   z2   = z2 + 0.0002 
 end do
 var1     = (1.0+z1) * (2997.92/h) * var1 !here var1 is luminosity distance. Varified by online calculator
 ext_factor  =  (1.0+z1) / ( 4*3.1415 * (var1*Mpc)**2.0 )  ! megaparsec in cm
         !The last (1.0+z1) factor accounts for K correction shinking of observed frequency band. 
end function ext_factor_luminosity


subroutine Total_Synchrotron_Intensity_LoS_integral(props,data,iz,steps_xp,steps_yp, steps_los_orig, d_los, hfac, resI, resQ, resU)
!This subroutine computes physical quantities along all possible LoSs through the galaxy. 
!The LoS is defined by an angle it makes with with normal to galaxy disc
!First a plane is defined perpenticular to LoS. This helps not miss any LoS passing through the galaxy. 
!Then the location where the LoS enters and leaves the galaxy is found. 
!Various quantities are then integrated along the LoS. 
    use input_constants
    use messages
    use global_input_parameters, only: number_of_redshifts
    use input_constants
    use FRB
    use interpolation
    type(Galaxy_Properties), intent(in) :: props
    type(LoS_data), intent(inout) :: data

    integer ::  j, iz, this_RM
    double precision, intent(inout) :: resI, resQ, resU 

    integer, intent(in) :: steps_xp, steps_yp,  steps_los_orig
    double precision, intent(in) ::  d_los, hfac
  
    real(fgsl_double) :: rlim(2), plim(2)
    integer:: rint, pint
    double precision :: res, lambda_m,lambdasq_m2
    character(len=20) :: form


 !Quantitites related to large scale field
    double precision, dimension(steps_rint) :: ray_kpc, Br_muG,  Bp_muG, Bz_muG, &
           B2_rval_muG2, ncr_cm3, n_cm3, ne_cm3, h_kpc, Brnd_muG, sigBr_muG
    double precision, dimension(steps_pint) :: nr_LoS, np_LoS, p_ay 
    double precision :: nz_LoS, Bmag2, drval, dpval
    double precision :: sint, cost, tant, signt, s, splus1by2, splus1by4, sminus1by2, factor_z, int_syn    
    double precision :: Brmid, Bpmid, Bzmid, ncrmid, nemid, hmid, rmid, sigBrmid, sigBrmidsq



    double precision, allocatable, dimension(:) ::  xp_ay, yp_ay  
    integer:: xint, yint, i_xp, i_yp, i_los, i_los_max, steps_los
    double precision :: rval, pval, xval, zval, yval, dxvalby2, dyval, dzvalby2, dlfac 
    double precision :: dxp, dyp
    double precision :: dlmin, dlmax,m_dl, c_dl, dlfacmin, dlfacmax 
    double precision :: xp_max, len_los, z_max, hr_max
    
    double precision, allocatable,  dimension(:,:) :: RM_ay, Stokes_I_ay, Stokes_U_ay, Stokes_Q_ay
    double precision, allocatable, dimension(:) :: RM_LoSay, Psi0_LoSay, IQU_integrant_LoSay  
    integer, allocatable,  dimension(:,:) :: int_LoS
    double precision, allocatable,  dimension(:,:,:) :: xz_entexit_ay 
   
    integer :: i_tmp, j_tmp 
    double precision :: mod_xlim, x_entry, x_exit, z_entry, z_exit, tanth  
    double precision :: Bperp_los_sq, Bpar_los, Bperp_los
    double precision :: r_max, lval, dlval
    double precision :: np, nr
    double precision :: var_tmp, z_tmp, x_tmp, Stokes_I_tmp, Stokes_Q_tmp, Stokes_U_tmp, RM_tmp, const,p0, Ke_fac 
    double precision :: tot_syn_I, psi0, psi
    double precision :: x1,rvalsq, z2, x3,r3, z4  
    double precision :: ratio_del_r, ratio_del_p 
    double precision :: By, Bxp

!    double precision :: f_Br_muG, f_Bp_muG, f_Bz_muG, f_ncr_cm3, f_h_kpc
!    double precision :: a_Br_muG, a_Bp_muG, a_Bz_muG, a_ncr_cm3, a_h_kpc

! limits of r and phi variables. The r-cordinate gives the radial distance to any point from center of galaxy.  
    rlim(1) = props%Rcyl(iz,1)    !The limit of r are ~0 to R_disk_of_galaxy
    rlim(2) = props%Rcyl(iz,props%n_grid) 
    plim(1) = 0.0 !limits of phi-variable 0<phi<2pi
    plim(2) = 6.28318530718
    drval  = ( rlim(2) - rlim(1) )/int(steps_rint-1) !dr - seperation between r-grid points
    dpval  = ( plim(2) - plim(1) )/int(steps_pint-1) !dphi - sepration between phi-grid points

!creating  grids  of r and phi 
    do rint =1, steps_rint, 1
         ray_kpc(rint) = rlim(1) + (rint-1) * drval
    end do
    do pint =1, steps_pint, 1
       p_ay(pint)  =  plim(1) + (pint-1) * dpval
    end do

!Occassionally you may get the last r value greater than disc radius (i.e. ray_kpc(steps_rint)> rlim(2))
!Then you have to set it back as ray_kpc(steps_rint) = rlim(2) to avoid interpolation errors
    if (ray_kpc(steps_rint)> rlim(2) ) then
          ray_kpc(steps_rint) = rlim(2)
    end if
    if (p_ay(steps_pint)> plim(2) ) then
          p_ay(steps_pint) = plim(2)
    end if

!On the r-phi grid points, magnetic field, scale height, and density are assigned from magnetizer output by interpolation
    call interpolate( props%Rcyl(iz,:), props%Br(iz,:), ray_kpc, Br_muG)
    call interpolate( props%Rcyl(iz,:), props%Bp(iz,:), ray_kpc, Bp_muG)
    call interpolate( props%Rcyl(iz,:), props%Bz(iz,:), ray_kpc, Bz_muG)
    call interpolate( props%Rcyl(iz,:), props%n(iz,:), ray_kpc, n_cm3)!number density of diffuse gas
    call interpolate( props%Rcyl(iz,:), props%h(iz,:), ray_kpc, h_kpc)
    ne_cm3   = n_cm3  * ionisation_fraction
    ncr_cm3  = n_cm3  * ionisation_fraction

    !initializing small scale field and galaxy angular velocity

    if ( data%inc_random_field == 1 ) then
       sigBr_muG = 0.5 * sqrt(4*pi * n_cm3 * Hmass) *v0_ISM_cms * mkG_G/sqrt(3.0) !mkG_G - gauss to micro gauss
    end if
    if ( data%inc_random_field == 2 ) then
       tau0_myr = 1.0
       Lu_kpc   = 1.0
       br2_muG2 = 0.0
       call interpolate( props%Rcyl(iz,:), props%shear_galaxy(iz,:), ray_kpc, sheargal)
       Uz_kms     = props%Uz(iz)
       Sfac       = ( 1.0 + ABS(sheargal)*tau0_myr)
       Sfacsq     = Sfac * Sfac
       const_factor_b = tau0_myr * Uz_kms / Lu_kpc
    end if

! Wavelength of observation and power law index of electron/CR energy distribution 
    lambda_m     = data%wavelength/(1d0 + props%z(iz))
    lambdasq_m2  = lambda_m*lambda_m

    s          = data%alpha
    splus1by4  = (s+1.0)/4.0
    splus1by2  = 2.0*splus1by4  ! power of B_perp
    sminus1by2 = (s - 1.0)/2.0
    p0       = (s+1.0)/(s+7.0/3.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! The Line of Sight (LoS)  !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Cosine, Sine and tangent of  angle that LoS makes with z-axis
    if ( abs(data%theta) <= 1.e-10) then !If LoS angle is <1.e-10, it is set to that value to avoid numerical issues.  
      data%theta = 1.e-10
    end if 
    cost = COS(data%theta) 
    sint = SIN(data%theta)
    tant = TAN(data%theta)
    signt = sign(1.d0, sint)

! We assume the xy plane coincides with galaxy disc plane. 
! The LoS is then assumed to be contained in xz plane so that its y component is always zero (nx, ny=0, nz) 
! The components of LoS unit vector in this case function of phi variable.    
! Below arrays are used to to store components (nr, np, nz) of LoS unit vectors as a function of phi 
    do pint =1, steps_pint, 1
      nr_LoS(pint) = COS( p_ay(pint) ) 
      np_LoS(pint) = SIN( p_ay(pint) ) 
    end do
!Components of LoS unit vectors (nr, np, nz)
    nr_LoS = sint * nr_LoS 
    np_LoS = -sint * np_LoS 
    nz_LoS = cost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!To integrate quantities along the LoS!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!First two variables define the grid
!Last two variables define integration along LoS
!    steps_xp = 150 !Should be even
!    steps_yp = 150 !Should be even
!    steps_los_orig = 100
!    d_los     = 1.e-2 

    
!We assume galaxy is a disc with radius r_max and height z_max. 
!limit of height of galactic disc. The z-integration is limited up few times maximum scale height. 
    hr_max=  MAXVAL(h_kpc) !maximum of scale height
    z_max   = hfac*hr_max
    r_max   = rlim(2)

!Below defines an X'Y' plane perpenticular to Los. Y' axis coincides with Y axis of galactic plane
!All LoS that pass through the galaxy will also pass through X'Y' plane. 
       !The maximum of X' value in so that all the LoSs through galaxy are counted
   xp_max = r_max * abs(cost) + z_max* abs(sint)
          ! Allocating X' and Y' grids with definite size. Each point on grid will contain an LoS
          ! If integral is not sufficienty accurate, then it will be done with finer grid.
   if (modulo(steps_xp,2) /= 0 .or. modulo(steps_yp,2) /= 0) then 
       print *, 'steps_xp and steps_yp should be even integers'
       stop
    end if  
    allocate( xp_ay (steps_xp ) )
    allocate( yp_ay (steps_yp ) )
            ! Assign values of X' and Y' on respective grids
    dxp = 2 * xp_max/(steps_xp-1)  !grid seperation for X' grid
    dyp = 2 * r_max/(steps_yp-1) !grid seperation for Y' grid. Maximum of Y'= maximum of Y = r_max
!creating the x'y' grid 
    do i_xp =1, steps_xp, 1
         xp_ay(i_xp) = -xp_max + dxp * (i_xp-1)
    end do
    do i_yp =1, steps_yp, 1
         yp_ay(i_yp) = -r_max + dyp * (i_yp-1)
    end do
!below steps to avoid floating point errors. Occassionally you will get  xp_ay(steps_xp) >xp_max
    if (xp_ay(steps_xp) >xp_max) then 
      xp_ay(steps_xp) = xp_max
    end if
    if (yp_ay(steps_yp) >r_max) then 
      yp_ay(steps_yp) = r_max
    end if

!Some LoS will not pass through the galaxy. Below we find those LoS and set int_LoS=0
!For LoSs that enter the galaxy, they enter either through:
           ! z=-z_max plane or
           ! r=r_max  plane  (the entry point depends of yval also) or 
           ! it will not pass through the galaxy 
!These LoSs will exit the galaxy either through:
           ! z=-z_max plane or
           ! r=r_max plane (the exit point depends of yval also) 
!We start with the LoS passing through negative values of X' and successively find all relevant LoS

!int_LoS(i,j) = 1 means the LoS passes through galaxy. Otherwise, it doesn't pass through it
! xz_entexit_ay stores x_entry, z_entry, x_exit, z_exit for each LoS passing through galaxy
! For index k=1 to 4, 1st index -> x_entry ; 2nd index -> z_entry; 3rd index -> x_exit; 4th index -> z_exit
  allocate( int_LoS (steps_xp, steps_yp ) )
  allocate( xz_entexit_ay (steps_xp, steps_yp, 4 ) )

!Looping over grid points to find which LoS enter the galaxy as well as their entry and exit points
!We need loop over only a quarter of grid points on X'Y' plane due to symmetry
  do i_xp =1, steps_xp/2, 1
    xval = xp_ay(i_xp) * cost !Corresponding x value on XZ plane 
    zval = -xp_ay(i_xp) * sint !Corresponding z value on XZ plane
    do i_yp =1, steps_yp/2, 1 
       yval = yp_ay(i_yp) !Fixed y value of LoS 
!To check if the LoS entering galaxy through z=-z_max plane. 
       x1   =  xval + (-z_max - zval )*tant !value of x for a given LoS at z=-z_max plane 
       rvalsq =  x1**2.0 + yval**2.0  !Distance of LoS from center at z=-z_max plane which is less than r_max 
!To check if LoS entering galaxy through r=r_max  
       mod_xlim = -signt*sqrt(r_max*r_max - yval*yval) !At this value of x LoS enters the galaxy (if at all it enters) 
       z2 =  zval + (mod_xlim - xval )/tant !z value fo an LoS at r=r_max plane 
!If LoS with (x',y) passes through the galaxy, finds x_ent & z_ent it as well as LoS with (x',-y) 
!Also find x_exit and z_exit for LoSs with (-x',y) and (-x',-y)
       if (rvalsq <= r_max*r_max ) then  ! If LoS entering galaxy through z=-z_max plane
         int_LoS(i_xp, i_yp)           = 1
         xz_entexit_ay(i_xp, i_yp, 1)  = x1 
         xz_entexit_ay(i_xp, i_yp, 2)  = -z_max 
       else if (z2 >= -z_max .and. z2 <= z_max ) then ! LoS entering galaxy through r=r_max and phi>pi plane
         int_LoS(i_xp, i_yp)           = 1
         xz_entexit_ay(i_xp, i_yp, 1)  = mod_xlim 
         xz_entexit_ay(i_xp, i_yp, 2)  = z2 
       else 
         int_LoS(i_xp, i_yp)           = 0
       end if         

       if (int_LoS(i_xp, i_yp) == 1 ) then 
         x3     =  xval + (z_max - zval )*tant !value of x for a given LoS at z=z_max plane 
         rvalsq =  x3**2.0 + yval**2.0  !Distance of LoS from center at z=-z_max plane. This should be less than r_max

!For LoS exiting galaxy through r=r_max and phi< pi plane.  
         mod_xlim = signt*sqrt(r_max*r_max - yval*yval) !At this value of x LoS exits the galaxy  
         z4       = zval + (mod_xlim - xval )/tant 
         if (rvalsq <= r_max*r_max ) then  ! LoS exits galaxy through z=z_max plane
            xz_entexit_ay(i_xp, i_yp, 3)      = x3 
            xz_entexit_ay(i_xp, i_yp, 4)      = z_max 
         else if (z4 <= z_max .and. z4 >= -z_max ) then ! LoS exiting galaxy through r=r_max and phi<pi plane
            xz_entexit_ay(i_xp, i_yp, 3)      = mod_xlim 
            xz_entexit_ay(i_xp, i_yp, 4)      = z4 
         else
            print *, 'Error: LoS entering a galaxy is not exiting'
           ! print *, x3, rvalsq, yval, mod_xlim, z4, z_max 
         end if
       end if    

         !If LoS with i_xp, i_yp passes through the galaxy, other three LoS also passes through it
         !Their entry and exit points are related to the original LoS passing through i_xp, i_yp
       int_LoS(i_xp, steps_yp-i_yp+1)             = int_LoS(i_xp, i_yp) 
       int_LoS( steps_xp-i_xp+1, i_yp)            = int_LoS(i_xp, i_yp) 
       int_LoS( steps_xp-i_xp+1, steps_yp-i_yp+1) = int_LoS(i_xp, i_yp)
       if (int_LoS(i_xp, i_yp) ==1 ) then 
          !los with i_yp ->  steps_yp-i_yp+1
         xz_entexit_ay(i_xp, steps_yp-i_yp+1, 1) =  xz_entexit_ay(i_xp, i_yp, 1)
         xz_entexit_ay(i_xp, steps_yp-i_yp+1, 2) =  xz_entexit_ay(i_xp, i_yp, 2)
         xz_entexit_ay(i_xp, steps_yp-i_yp+1, 3) =  xz_entexit_ay(i_xp, i_yp, 3)  
         xz_entexit_ay(i_xp, steps_yp-i_yp+1, 4) =  xz_entexit_ay(i_xp, i_yp, 4) 
          !los with i_xp -> steps_xp-i_xp+1
         xz_entexit_ay( steps_xp-i_xp+1, i_yp, 3)  = -xz_entexit_ay(i_xp, i_yp, 1)
         xz_entexit_ay( steps_xp-i_xp+1, i_yp, 4)  = -xz_entexit_ay(i_xp, i_yp, 2) 
         xz_entexit_ay( steps_xp-i_xp+1, i_yp, 1)  = -xz_entexit_ay(i_xp, i_yp, 3) 
         xz_entexit_ay( steps_xp-i_xp+1, i_yp, 2)  = -xz_entexit_ay(i_xp, i_yp, 4)  
          !los with i_yp ->  steps_yp-i_yp+1 and  i_xp -> steps_xp-i_xp+1
         xz_entexit_ay( steps_xp-i_xp+1, steps_yp-i_yp+1, 3) = -xz_entexit_ay(i_xp, i_yp, 1)
         xz_entexit_ay( steps_xp-i_xp+1, steps_yp-i_yp+1, 4) = -xz_entexit_ay(i_xp, i_yp, 2)
         xz_entexit_ay( steps_xp-i_xp+1, steps_yp-i_yp+1, 1) = -xz_entexit_ay(i_xp, i_yp, 3)   
         xz_entexit_ay( steps_xp-i_xp+1, steps_yp-i_yp+1, 2) = -xz_entexit_ay(i_xp, i_yp, 4) 
       end if 
    end do
  end do   

  
 !Array to store final integrated products along each LoS
  allocate( RM_ay (steps_xp, steps_yp ) )
  allocate( Stokes_I_ay (steps_xp, steps_yp ) )
  allocate( Stokes_Q_ay (steps_xp, steps_yp ) )
  allocate( Stokes_U_ay (steps_xp, steps_yp ) )
  RM_ay = 0.0
  Stokes_I_ay = 0.0
  Stokes_Q_ay = 0.0
  Stokes_U_ay = 0.0

! (x_entry, z_entry) and (x_exit, z_exit) stores (x, z) locations where LoS passing through X'Y' plane intersects the galaxy.
! An LoS can potentially intersect galaxy one two of the four possible faces. Below we finds those locations. 
   do i_xp =1, steps_xp, 1
     xval = xp_ay(i_xp) * cost !y value of LoS on the plane perpenticular to LoS
     zval = xp_ay(i_xp) * (-sint) !z value of LoS on the plane perpenticular to LoS
     do i_yp =1, steps_yp, 1
       yval = yp_ay(i_yp)
       if (int_LoS(  i_xp, i_yp) == 1) then   
         x_entry = xz_entexit_ay(i_xp, i_yp,1)
         z_entry = xz_entexit_ay(i_xp, i_yp,2)
         x_exit  = xz_entexit_ay(i_xp, i_yp,3)
         z_exit  = xz_entexit_ay(i_xp, i_yp,4)
         
         steps_los = steps_los_orig
                   !Total length of LoS along the galaxy
         len_los = sqrt( ( x_exit-x_entry )**2.0 + (z_exit-z_entry)**2.0 )
         dlval = MIN(d_los, len_los/steps_los)
         if (dlval <= 0.0) then 
            dlval = 1.e-5
         end if
            !reassigning steps_los 
         if (len_los <= 0.0) then
           steps_los = 0 
         else
           steps_los = INT(len_los/dlval+1) 
         end if

            !Arrays to store values of RM, angle and QU integrant along a given LoS 
         allocate( RM_LoSay (steps_los ) )
         allocate( Psi0_LoSay (steps_los ) )
         allocate( IQU_integrant_LoSay (steps_los ) )
         RM_LoSay            = 0.0
         Psi0_LoSay          = 0.0
         IQU_integrant_LoSay = 0.0

                !, Half of the change in dx and dz, for one step along LoS
         dxvalby2 = dlval * sint/2.0
         dzvalby2 = dlval * cost/2.0
          !values x and z at which LoS enters the galaxy
         xval  = x_entry
         zval  = z_entry
           !Initializing various parameters
         lval         = 0
         RM_tmp       = 0.0
         Stokes_I_tmp = 0.0

         i_los = 1
         do while( lval < len_los-dlval )
                !x,z, r, phi values along the LoS
          xval = xval + dxvalby2  !dxval = dlval*sint/2 
          zval = zval + dzvalby2   !dzval = dzval*cost/2 
          rval =  SQRT( xval**2.0 + yval**2.0 )
          
          !print *, 'lenth_los, rval, zval', len_los, xval, yval,  rval, zval
          if ( rval< r_max .and. abs(zval) <z_max  ) then 
             !rval is located between rint and rint+1 on input r-value array. 
             rint = INT( (rval-rlim(1)) / drval) + 1  !
             if (rint==steps_rint) then !to avoid errors below during interpolation
                rint = steps_rint -1
             end if
!Interpolating to find Br, Bphi, Bz, ncr, h etc on the above r and phi
             ratio_del_r =  (rval - ray_kpc(rint)) / ( ray_kpc(rint+1) - ray_kpc(rint) )
             Brmid  =  Br_muG(rint)  + ratio_del_r * (Br_muG(rint+1)  - Br_muG(rint))  
             Bpmid  =  Bp_muG(rint)  + ratio_del_r * (Bp_muG(rint+1)  - Bp_muG(rint))  
             Bzmid  =  Bz_muG(rint)  + ratio_del_r * (Bz_muG(rint+1)  - Bz_muG(rint))  
             hmid   =  h_kpc(rint)   + ratio_del_r * (h_kpc(rint+1)   - h_kpc(rint))
             nemid  =  ne_cm3(rint)  + ratio_del_r * ( ne_cm3(rint+1) - ne_cm3(rint) ) 
             ncrmid =  ncr_cm3(rint)  + ratio_del_r * ( ncr_cm3(rint+1) - ncr_cm3(rint) )
             Bmag2  =  Brmid**2.0 + Bpmid**2.0 + Bzmid**2.0 !B^2

!Find various averages of small scale field from br and omega_r
             if ( data%inc_random_field ==1 ) then
              sigBrmidsq =( sigBr_muG(rint) + ratio_del_r * (sigBr_muG(rint+1) - sigBr_muG(rint)) )**2.0
             end if

             if ( data%inc_random_field ==2 ) then
               br2mid_muG2   = br2_muG2(rint) + ratio_del_r * ( br2_muG2(rint+1) - br2_muG2(rint) )    
               Sfac_mid    = Sfac(rint) + ratio_del_r * ( Sfac(rint+1) - Sfac(rint) )
               Sfacsq_mid  = Sfacsq(rint) + ratio_del_r * ( Sfacsq(rint+1) - Sfacsq(rint) )
             
               bp2mid_muG2   = br2mid_muG2 * Sfacsq_mid 
               bz2mid_muG2   = br2mid_muG2 * const_factor_b * (1.0+Sfac_mid)**2.0
               brbpmid_muG2  = br2mid_muG2 * Sfac_mid
               brbzmid_muG2  = (br2mid_muG2 + brbpmid_muG2)  * const_factor_b 
               bpbzmid_muG2  = (bp2mid_muG2 + brbpmid_muG2)  * const_factor_b
             end if


             if (CR_B_equipartition) then 
                Ke_fac = Ke_equipart( Bmag2, s)  
             else 
                Ke_fac = ncrmid
             end if

             pval =  ATAN2( yval, xval )  !returns correct value of phi between 0 and 2 pi
             if (pval <0.0) then
                pval = pval + 2*pi
             end if
             !pval is located between pint and pint+1 on the phi-value array. 
             pint = INT(pval/dpval) + 1 
             if (pint==steps_pint) then  !to avoid errors below during interpolation
              pint = steps_pint -1
             end if
!Interpolating to find components of unit LoS vector along r, phi, z directions.
             ratio_del_p =  (pval - p_ay(pint)) / ( p_ay(pint+1) - p_ay(pint) )
             np = np_LoS(pint) + ratio_del_p * (np_LoS(pint+1)-np_LoS(pint))     
             nr = nr_LoS(pint) + ratio_del_p * (nr_LoS(pint+1)-nr_LoS(pint))    

!Parelle and perpenticular components of B wrt to LoS at r and phi 
             Bpar_los      = Brmid*nr + Bpmid*np + sign(Bzmid, zval) * nz_LoS ! B.n 
             Bperp_los_sq  = Bmag2 -   Bpar_los**2.0  !A new formulae for Bpar_los is used
             By       =  Brmid*sin(pval) + Bpmid*cos(pval)
             Bxp      = (Brmid*cos(pval) - Bpmid*sin(pval) )*cost -  sign(Bzmid, zval)*sint
                  !    Below exponential factors are z-dependence of density (which is proportional to sqaure of field)
                                 ! and square of field
             if (data%inc_random_field == 1) then
                 Bperp_los_sq = (Brmid**2.0 + sigBrmidsq)* (1 - ( sint*cos(pval) )**2.0 )  & 
                            + (Bpmid**2.0 + sigBrmidsq) * (1 - ( sint*sin(pval) )**2.0 )  &
                            + (Bzmid**2.0 + sigBrmidsq) * sint**2.0  &
                            + sint**2.0 * sin(2.0*pval) * Brmid * Bpmid &
                            - sign(Bzmid,zval) * 2*sint*cost * ( Brmid*cos(pval) - Bpmid*sin(pval) ) 
             end if 
!small scale field in xp, yp coordinates                                 
             if ( data%inc_random_field == 2  ) then
               bxp2_muG2 = br2mid_muG2 * cos(pval)**2.0 + bp2mid_muG2 * sin(pval) **2.0 + &
                                         brbpmid_muG2 * sin(2.0*pval) 
               byp2_muG2 =  br2mid_muG2 * (sin(pval)*cost)**2.0 + bp2mid_muG2 * (cost*cos(pval)) **2.0 &
                   + bz2mid_muG2*sint**2.0 &
                   + ( brbzmid_muG2 * sin(pval) + bpbzmid_muG2 * cos(pval) ) * sin(2*data%theta) &
                   +  brbpmid_muG2 * sin(2.0*pval) * cost**2.0
             end if

!end small scale field

             if (CR_B_equipartition) then 
              IQU_integrant_LoSay(i_los) = Bperp_los_sq**splus1by4 * Ke_fac * exp(-abs(4.0*zval/hmid))*dlval 
              RM_tmp                     = RM_tmp +  nemid * Bpar_los * exp(-abs(3.0*zval/hmid) ) * dlval
             else 
              IQU_integrant_LoSay(i_los) = Bperp_los_sq**splus1by4 * Ke_fac * exp(-abs(3.0*zval/hmid))*dlval   
              RM_tmp                     = RM_tmp +  nemid * Bpar_los * exp(-abs(2.0*zval/hmid) ) * dlval
             end if

             Stokes_I_tmp      = Stokes_I_tmp + IQU_integrant_LoSay(i_los)  
             RM_LoSay(i_los)   = RM_tmp
             Psi0_LoSay(i_los) =  pi/2d0 + atan2(By,Bxp)
             i_los             = i_los + 1

          end if !end of if the point along LoS is inside the galaxy statement
          xval = xval + dxvalby2  
          zval = zval + dzvalby2
          lval = lval + dlval
         end do  !end of integration along LoS array
         i_los_max         = i_los-1
         RM_ay(i_xp, i_yp)       = 812.0 * RM_tmp  * lambdasq_m2 
         RM_LoSay                = 812.0 * RM_LoSay * lambdasq_m2   
         IQU_integrant_LoSay     = p0 * IQU_integrant_LoSay !Now this is QU_integrantLoSay. Not used anymore to compute I 
                                                              !p0 ~ 0.75  
! Calculating stokes Q and U parameters. For this We need RM as a function of the path length. 
! This is be done by doing another LoS integral as below
         Stokes_Q_tmp = 0.0
         Stokes_U_tmp = 0.0
         do i_los =1, i_los_max, 1
            psi          = Psi0_LoSay(i_los) + ( RM_ay(i_xp, i_yp) - RM_LoSay(i_los) ) 
            Stokes_Q_tmp =  Stokes_Q_tmp + IQU_integrant_LoSay(i_los) * COS( 2.0 * psi )
            Stokes_U_tmp =  Stokes_U_tmp + IQU_integrant_LoSay(i_los) * SIN( 2.0 * psi )
         end do

!In the end assigning integrated intensities to grid 
         Stokes_I_ay(i_xp, i_yp) = Stokes_I_tmp  
         Stokes_Q_ay(i_xp, i_yp) = Stokes_Q_tmp 
         Stokes_U_ay(i_xp, i_yp) = Stokes_U_tmp 

         deallocate(RM_LoSay)
         deallocate(Psi0_LoSay) 
         deallocate(IQU_integrant_LoSay)
       end if ! if statement to check if given LoS passes through galaxy (int_LoS ==1) 
     end do !end of i_y loop
  end do !end of i_x loop

!Converting quantities to Physical units 
 if ( CR_B_equipartition ) then
  const  =  ext_factor_luminosity(props%z(iz)) * emissivity_syn_constant(s) &
                      * dxp*dyp *lambda_m** sminus1by2 
  Stokes_I_ay =  const * Stokes_I_ay 
  Stokes_Q_ay =  const * Stokes_Q_ay 
  Stokes_U_ay =  const * Stokes_U_ay 
 else 
  const  =  ext_factor_luminosity(props%z(iz)) * dxp*dyp *lambda_m** sminus1by2
  Stokes_I_ay =  const * Stokes_I_ay 
  Stokes_Q_ay =  const * Stokes_Q_ay 
  Stokes_U_ay =  const * Stokes_U_ay 
 end if 

!Integrating to find total stokes I, Q, U. Total Stokes I will be compared with analytic result
 resI = 0.0
 resQ = 0.0
 resU = 0.0
 do i_xp =1, steps_xp, 1
  do i_yp =1, steps_yp, 1
    resI =  resI + Stokes_I_ay(i_xp, i_yp) 
    resQ =  resQ + Stokes_Q_ay(i_xp, i_yp) 
    resU =  resU + Stokes_U_ay(i_xp, i_yp) 
  end do
 end do

!to write everything to text file.
 write_IQU(3*iz-2) = resI
 write_IQU(3*iz-1) = resQ
 write_IQU(3*iz)   = resU

end subroutine Total_Synchrotron_Intensity_LoS_integral


!write final output to text file. These can be called from Observables. 
!Also possible to write output at all z for a given galaxy to a single file
function write_I_textfile( igal, iz, n_z) result(write_I)
 integer, intent(in) :: igal, iz, n_z
 integer :: write_I
 open(37,file=trim('./output_text_files')//'/out.dat', FORM='FORMATTED', status='old', position="append")
     write( 37, * ), igal, write_IQU(1:3*n_z) 
 close(37)
 write_IQU = 0.0
 write_I   = 1
end function write_I_textfile


function init_write_I_textfile( ) result(write_I)
integer :: write_I
! Opening the file.  
  open(37,file=trim('./output_text_files')//'/out.dat', FORM='FORMATTED', status='replace')
     write( 37, * ), '#' 
  close(37)
  write_I = 1
end function init_write_I_textfile




end module LoSintegrate_aux
