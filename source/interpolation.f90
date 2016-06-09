! Contains a module which implements the rescaling of the f-array
! using interpolation routines from the FGSL
module interpolation
  implicit none
  private
  public :: interpolate,

  contains

  subroutine interpolate(x, y, xi, yi, method)
    ! Interpolates using FGSL
    ! Input: x, y -> 1d-arrays containing the known values of the function
    !        xi -> 1d-array with the point one want to interpolate into
    !      optional, method -> either 'linear', 'polynomial' or 'cubic_spline'
    ! Output: yi -> 1d-array containing the interpolated values
    use FGSL
    real(fgsl_double), dimension(:), intent(in) :: x, y
    real(fgsl_double), dimension(:), intent(in) :: xi
    real(fgsl_double), dimension(size(xi)), intent(out) :: yi
    character(len=*), intent(in), optional :: method
    integer(fgsl_size_t) :: n, i
    integer(fgsl_int) :: status
    type(fgsl_interp_accel) :: acc
    type(fgsl_spline) :: spline
    type(fgsl_interp_type) :: interp_type

    if (.not.present(method)) then
      interp_type = fgsl_interp_linear
    else
      select case (trim(method))
        case('linear')
          interp_type = fgsl_interp_linear
        case('polynomial')
          interp_type = fgsl_interp_polynomial
        case('cubic_spline')
          interp_type = fgsl_interp_cspline
        case default
          print *, 'interpolate: Unrecognized option ',method
      end select
    endif

    ! Prepares the interpolation/spline routine
    n = size(x)
    spline =  fgsl_spline_alloc(interp_type, n)
    status = fgsl_spline_init(spline, x, y, n)
    ! Stores the interpolated values
    do i=1, size(xi)
      yi(i) = fgsl_spline_eval(spline, xi(i), acc)
    enddo
    ! Frees FGSL stuff
    call fgsl_spline_free (spline)
    call fgsl_interp_accel_free (acc)
  end subroutine interpolate

end module interpolation
