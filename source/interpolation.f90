! Contains a module which implements the rescaling of the f-array
! using interpolation routines from the FGSL
module interpolation
  implicit none
!   private
!   public :: interpolate,

  contains

  subroutine interpolate(x, y, xi, yi, method)
    ! Interpolates using FGSL
    ! Input: x, y -> 1d-arrays containing the known values of the function
    !        xi -> 1d-array with the point one want to interpolate into
    !        optional, method -> either 'linear', 'polynomial', 'cubic_spline'
    !                            or 'akima'
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
        case('akima')
          interp_type = fgsl_interp_akima
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

  function rescale_array(y, n_new) result(y_new)
    ! Interpolates the n-array y into an n_new array, with the same endpoints
    double precision, dimension(:), intent(in) :: y
    double precision, dimension(size(y)) :: x
    integer, intent(in) :: n_new
    integer :: n, i
    double precision, dimension(n_new) :: y_new
    double precision, dimension(n_new) :: x_new

    n = size(y)
    ! If the old and new grid sizes are the same, do nothing
    if (n==n_new) then
      y_new = y
      return
    endif

    ! Generates x, and x_new
    ! I.e two 0 to 1 uniform arrays, with the correct number of points
    x_new = 0.0d0
    x = 0.0d0
    do i=1, n
      x(i) = dble(i-1)/dble(n-1)
    enddo
    do i=1, n_new
      x_new(i) = dble(i-1)/dble(n_new-1)
    enddo
    ! Calls the interpolation routine
    call interpolate(x, y, x_new, y_new)

  end function

  subroutine rescale_array_in_place(y, n)
    ! Interpolates the n-array into the f_rescaled array
    double precision, dimension(:), allocatable, intent(inout) :: y
    double precision, dimension(:), allocatable :: ytmp
    integer, intent(in) :: n
    ! Backups the array contents
    call move_alloc(y,ytmp)
    ! Allocates the new array
    allocate(y(n))
    ! Interpolates
    y = rescale_array(ytmp, n)
  end subroutine rescale_array_in_place



  subroutine rescale_f_array(f, n)
    ! Interpolates the f-array into the f_rescaled array
    double precision, dimension(:,:), allocatable, intent(inout) :: f
    double precision, dimension(:,:), allocatable :: ftmp
    integer, intent(in) :: n

    integer, dimension(2) :: shape_f
    integer :: nvar, nx, i

    shape_f = shape(f)
    nx   = shape_f(1)
    nvar = shape_f(2)

    call move_alloc(f,ftmp)
    allocate(f(n,nvar))

    do i=1,nvar
      f(:,i) = rescale_array(ftmp(:,i), n)
    enddo

  end subroutine rescale_f_array

end module interpolation
