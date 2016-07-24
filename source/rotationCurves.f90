module rotationCurves
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
    double precision, dimension(size(r)),intent(out) :: Omega, Shear
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
      Shear = 0
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
              ! Check this!
            / ( I0(y50) * K0(y50)  &
              - I1(y50) * K1(y50))
      A(i) = sqrt(A(i))
    end do
    Omega = A*v_disk/r_disk

    do i=1,size(r)
      Shear(i) =    I1(y(i)) * K0(y(i))                       &
                  - I0(y(i)) * K1(y(i))                       &
                  -0.5d0 * K1(y(i)) *( I0(y(i)) + I2(y(i)) )  &
                  +0.5d0 * I1(y(i)) *( K0(y(i)) + K2(y(i)) )
    end do

    Shear = Shear/A/2d0*v_disk/r_disk
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
    double precision, dimension(size(r)),intent(out) :: Omega, Shear
    double precision, parameter :: a_to_r50 = (sqrt(2.0d0)-1.0d0)
    double precision, dimension(size(r)) :: y, a
    double precision :: Constant, y50

    ! Traps and deals with bulge-less galaxies
    if (v_bulge<1e-3) then
      Omega = 0d0
      Shear = 0d0
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

    dvdr = Constant/r_bulge * ( 0.5d0/sqrt(y)/(y+1d0) - sqrt(y)/(y+1d0)**2 )
!     dvdr = Constant/2d0/((y+1d0)*sqrt(y)) - Constant*sqrt(y)/(1d0+y)**2

    Shear = dvdr - Omega

  end subroutine bulge_rotation_curve


  subroutine halo_rotation_curve(rx, r_halo, v_halo, nfw_cs1, Omega, Shear)
    ! Computes the rotation curve associated with an NFW halo
    ! Warning: This ignores effects of adiabatic contraction!
    !          It need to be accounted for later
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
    use deriv
    implicit none
    double precision, intent(in) :: r_halo, v_halo, nfw_cs1
    double precision, dimension(:), intent(in)  :: rx
    double precision, dimension(size(rx)), intent(out) :: Omega, Shear
    double precision, dimension(size(rx)) :: v, dvdr, dv2dy
    double precision, dimension(size(rx)) :: y, B
    double precision :: A, c
    double precision, parameter :: B_TOL = -1e-7

    y = abs(rx) / r_halo
    c = 1d0/nfw_cs1

    A = v_halo**2 / (log(1d0+c)-c/(1d0+c))
    B = (log(1d0+c*y) - c*y/(1d0+c*y))/y

    if (any(B<B_TOL)) then
      stop 'halo_rotation_curve: error, halo properties lead to negative '&
            // 'rotation curves!'
    endif

    v = sqrt(abs(A*B))

    dv2dy = A*(-B +  (1d0-c)/(1d0+c*y) +c**2*y/(1d0+c*y)**2 )/y
    dvdr = dv2dy/2d0/v/r_halo

    Omega = v/rx
    Shear = dvdr - Omega
  end subroutine halo_rotation_curve

end module rotationCurves
