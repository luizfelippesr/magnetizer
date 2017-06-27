module rotationCurves
  use, intrinsic :: iso_c_binding
  use fgsl
  implicit none
  private

  public :: disk_rotation_curve, bulge_rotation_curve, halo_rotation_curve

contains
  subroutine disk_rotation_curve(r, r_disk, rmin, v_disk, Omega, Shear)
    ! Computes the rotation curve associated with an exponential disk
    ! Input:  r -> the radii where the rotation curve will be computed
    !         r_disk ->  the half mass radius of the disk
    !         v_disk -> the circular velocity at r_disk
    !         r_min -> the minimum radius to be considered
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:   V^2 \propto y^2 [I0(y)K0(y)-I1(y)K1(y)] where y=r/r_s
    !         r_s --> the scale radius
    !         For the shear, see  http://is.gd/MnDKyS
    ! Ref:    Binney & Tremaine or  Mo, Bosch & White

    use Bessel_Functions
    use input_constants
    implicit none
    double precision, intent(in) :: r_disk, v_disk, rmin
    double precision, dimension(:), intent(in)  :: r
    double precision, dimension(size(r)) :: A
    double precision, dimension(size(r)),intent(out) :: Omega
    double precision, dimension(size(r)),intent(out), optional :: Shear
    double precision, parameter :: TOO_SMALL=2e-7 ! chosen empirically
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, dimension(size(r)) :: y
    double precision :: rs, y50
    integer :: i

    ! Scale radius
    rs = r_disk*constDiskScaleToHalfMassRatio
    ! (Half-mass radius) / (scale radius)
    y50 = 1d0/constDiskScaleToHalfMassRatio

    ! Traps disks of negligible size
    if (r_disk < rmin) then
      if (present(Shear)) then
        Shear = 0
      endif
      Omega = 0
      return
    end if

    ! Sets y=r/rs, trapping r~0 (which breaks the evaluation
    ! of the Bessel functions)
    where (abs(r) > TOO_SMALL)
      ! The absolute value is taken to deal with the ghost zone
      y = abs(r) / rs
    elsewhere
      y = TOO_SMALL / rs
    endwhere

    do i=1,size(r)
      A(i) = (  I0(y(i)) * K0(y(i))            &
              - I1(y(i)) * K1(y(i)) )          &
            / ( I0(y50) * K0(y50)  &
              - I1(y50) * K1(y50))
      A(i) = sqrt(A(i))
    end do
    Omega = A*v_disk/r_disk

    if (present(Shear)) then
      do i=1,size(r)
        Shear(i) =    I1(y(i)) * K0(y(i))                       &
                    - I0(y(i)) * K1(y(i))                       &
                    -0.5d0 * K1(y(i)) *( I0(y(i)) + I2(y(i)) )  &
                    +0.5d0 * I1(y(i)) *( K0(y(i)) + K2(y(i)) )
      end do

      Shear = Shear/A/2d0*v_disk/r_disk
    endif
    return
  end subroutine disk_rotation_curve

  subroutine bulge_rotation_curve(r, r_bulge, v_bulge, Omega, Shear)
    ! Computes the rotation curve associated with an Hernquist profile
    ! Input:  r -> the radii where the rotation curve will be computed
    !         r_bulge ->  the half mass radius of the spheroid
    !         v_bulge -> the circular velocity at r_bulge
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:   V^2 \propto y/(y+1)^2  where y=r/a
    ! Ref:    http://adsabs.harvard.edu/abs/1990ApJ...356..359H

    implicit none
    double precision, intent(in) :: r_bulge, v_bulge
    double precision, dimension(:), intent(in)  :: r
    double precision, dimension(size(r)) :: v, dvdr
    double precision, dimension(size(r)),intent(out) :: Omega
    double precision, dimension(size(r)),intent(out), optional :: Shear
    double precision, parameter :: a_to_r50 = (sqrt(2.0d0)-1.0d0)
    double precision, dimension(size(r)) :: y, a
    double precision :: Constant, y50

    ! Traps and deals with bulge-less galaxies
    if (v_bulge<1e-3) then
      Omega = 0d0
      if (present(Shear)) Shear = 0d0
      return
    endif

    ! Computes characteristic radius
    a = a_to_r50*r_bulge
    y = abs(r)/a
    y50 = 1d0/a_to_r50 ! r_bulge/a

    Constant = v_bulge * (1d0 + y50) / sqrt(y50)

    ! v =  (y/a_to_r50) * (a_to_r50 +1d0)**2 * (y +1d0)**(-2)
    ! v = v_bulge * sqrt(v)
    v = Constant * sqrt(y)/(y+1d0)

    Omega = v/abs(r)

    if (present(Shear)) then
      dvdr = Constant/r_bulge * ( 0.5d0/sqrt(y)/(y+1d0) - sqrt(y)/(y+1d0)**2 )
      Shear = dvdr - Omega
    endif
  end subroutine bulge_rotation_curve

  subroutine halo_rotation_curve(r, baryon_fraction, r_halo, v_halo, nfw_cs1, &
                        r_bulge, v_bulge, r_disk, v_disk, rmin, Omega, contract)
    ! Computes the rotation curve associated with an NFW halo
    ! Warning: This ignores effects of adiabatic contraction!
    !          It needs to be accounted for later
    ! Input: r -> the radii where the rotation curve will be computed
    !        r_halo -> the virial radius of the halo (r_200)
    !        v_halo -> the circular velocity at r_halo (V_200)
    !        nfw_cs1 -> 1/c_s -> inverse of the NFW concentration parameter
    !                        i.e. (NFW scale radius)/(virial radius)
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:  V^2 = V_200^2 * {ln(1+cy) - cy/(1+cy)} /
    !                    {ln(1+c) - 1/(1+c)} / y
    ! Ref: NFW profile
    use messages
    implicit none
    double precision, intent(in) :: baryon_fraction, r_halo, v_halo, nfw_cs1
    double precision, intent(in) :: r_bulge, v_bulge, r_disk, v_disk, rmin
    double precision, dimension(:), intent(in)  :: r
    double precision, dimension(size(r)) :: v
    double precision, dimension(size(r)),intent(out) :: Omega
    logical, intent(in), optional :: contract
    logical :: contract_actual
    integer :: i

    if (present(contract)) then
      contract_actual = contract
      print *, 'contract',contract
    else
      contract_actual = .false.
    endif

    if (.not.contract_actual) then
      do i=1,size(r)
        v(i) = sqrt(1d0-baryon_fraction)*halo_velocity(r(i), r_halo, v_halo, nfw_cs1)
      end do
    else
      do i=1,size(r)
        v(i) = halo_velocity_contracted(abs(r(i)), baryon_fraction, r_halo, v_halo,&
                              nfw_cs1, r_bulge, v_bulge, r_disk, v_disk, rmin)
      end do
    endif

    Omega = v/r

  end subroutine halo_rotation_curve



  pure function halo_velocity(r, r_halo, v_halo, nfw_cs1) result(V)
    ! Computes the rotation curve associated with an NFW halo
    ! Warning: This ignores effects of adiabatic contraction!
    !          It needs to be accounted for later
    ! Input: r -> the radii where the rotation curve will be computed
    !        r_halo -> the virial radius of the halo (r_200)
    !        v_halo -> the circular velocity at r_halo (V_200)
    !        nfw_cs1 -> 1/c_s -> inverse of the NFW concentration parameter
    !                        i.e. (NFW scale radius)/(virial radius)
    !        baryon_fraction -> fraction of the mass of the halo in barions
    ! Output: velocity at r
    ! Info:  V^2 = V_200^2 * {ln(1+cy) - cy/(1+cy)} /
    !                    {ln(1+c) - 1/(1+c)} / y
    ! Ref: NFW profile
    use messages
    implicit none
    double precision, intent(in) :: r, r_halo, v_halo, nfw_cs1
    double precision :: V, y, A, B, v2, c
    double precision, parameter :: v2_tol = -1d0 ! (km/s)^2

    y = abs(r) / r_halo
    c = 1d0/nfw_cs1

    A = v_halo**2 / (log(1d0+c)-c/(1d0+c))
    B = (log(1d0+c*y) - c*y/(1d0+c*y))/y
    v2 = A*B

!     if (any(v2<v2_tol)) then
!       ! This tolerance was introduced to deal with (otherwise harmless)
!       ! numerical errors in the centre of the halo.
!       call message('halo_rotation_curve: error, halo properties lead to' &
!                    //' invalid rotation velocities! v2min=', minval(v2), &
!                    info=0)
!     endif

    V = sqrt(abs(v2))

  end function halo_velocity

  function halo_velocity_contracted(r, baryon_fraction,r_halo,v_halo,nfw_cs1, &
                              r_bulge, v_bulge, r_disk, v_disk, rmin) result(V)
    ! Computes the rotation curve associated with an NFW halo accounting for
    ! adiabatic contraction
    !          It needs to be accounted for later
    ! Input: r -> the radii where the rotation curve will be computed
    !        r_halo -> the virial radius of the halo (r_200)
    !        v_halo -> the circular velocity at r_halo (V_200)
    !        nfw_cs1 -> 1/c_s -> inverse of the NFW concentration parameter
    !                        i.e. (NFW scale radius)/(virial radius)
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:  V^2 = V_200^2 * {ln(1+cy) - cy/(1+cy)} /
    !                    {ln(1+c) - 1/(1+c)} / y
    ! Ref: NFW profile
    use messages
    use root_finder, only: FindRoot
    implicit none
    double precision, intent(in) :: r_halo, v_halo, nfw_cs1, baryon_fraction
    double precision, intent(in) :: r_bulge, v_bulge, r_disk, v_disk, rmin
    double precision, intent(in)  :: r
    double precision :: v, y, B, v2, A, c, r0
    real(fgsl_double), target, dimension(10) :: parameters_array
    double precision, parameter :: v2_tol = -1d0 ! (km/s)^2
    double precision :: fun_min
    type(c_ptr) :: params_ptr
    integer, parameter :: MAX_INTERVAL_EXPONENT=8
    integer :: i

    ! Interfacing gymnastics...
    parameters_array(1)  = r_halo
    parameters_array(2)  = v_halo
    parameters_array(3)  = nfw_cs1
    parameters_array(4)  = r_bulge
    parameters_array(5)  = v_bulge
    parameters_array(6)  = r_disk
    parameters_array(7)  = v_disk
    parameters_array(8)  = rmin
    parameters_array(9)  = baryon_fraction
    parameters_array(10) = r
    ! Prepares a pointer to parameters_array
    params_ptr = c_loc(parameters_array)

    fun_min = halo_velocity_contracted_fun(r, params_ptr)
    if (fun_min>0d0) then
      ! Avoids shell crossing!
      V = halo_velocity(r0, r_halo, v_halo, nfw_cs1)

    else
      do i=1,MAX_INTERVAL_EXPONENT
        ! Tries to find an interval which contains a root
        if (halo_velocity_contracted_fun(r*1.5d0**i, params_ptr)/fun_min <0) exit
      enddo

      r0 = FindRoot(halo_velocity_contracted_fun, params_ptr, [r, r*1.5d0**i])
      V = halo_velocity(r0, r_halo, v_halo, nfw_cs1)
    endif
  end function halo_velocity_contracted

  function halo_velocity_contracted_fun(r0, params) bind(c)
    ! f(x) = r0^2*Vh^2(r0) - r^2(f_b*Vh^2(r)*(r0/r)+Vd^2(r)+Vb^2(r))
    ! root finder GSL routine
    real(c_double), value :: r0
    type(c_ptr), value :: params
    real(c_double) :: halo_velocity_contracted_fun
    real(fgsl_double) :: r, r_halo, v_halo, nfw_cs1, rmin
    real(fgsl_double) :: r_bulge, v_bulge, r_disk, v_disk, f_b
    real(fgsl_double) :: V0, Vr, Vd2, Vb2
    real(fgsl_double), dimension(1) :: Omega_tmp
!     real(fgsl_double), target, dimension(10) :: p
    real(fgsl_double), pointer, dimension(:) :: p

    call c_f_pointer(params, p, [10])
    r_halo = p(1)
    v_halo = p(2)
    nfw_cs1 = p(3)
    r_bulge = p(4)
    v_bulge = p(5)
    r_disk = p(6)
    v_disk = p(7)
    rmin = p(8)
    f_b = p(9)
    r = p(10)

    V0 = halo_velocity(r0, r_halo, v_halo, nfw_cs1)
    Vr = halo_velocity(r,  r_halo, v_halo, nfw_cs1)
    call bulge_rotation_curve([r], r_bulge, v_bulge, Omega_tmp)
    Vb2 = (Omega_tmp(1)*r)**2
    call disk_rotation_curve([r], r_disk, rmin, v_disk, Omega_tmp)
    Vd2 = (Omega_tmp(1)*r)**2

    halo_velocity_contracted_fun = r0**2*V0**2 - r**2*((1d0-f_b)*Vr**2+Vd2+Vb2)

  end function halo_velocity_contracted_fun

end module rotationCurves
