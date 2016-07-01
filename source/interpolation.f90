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

  function rescale_array_old(y, n_new, method) result(y_new)
    ! Interpolates the n-array y into an n_new array, with the same endpoints
    double precision, dimension(:), intent(in) :: y
    double precision, dimension(size(y)) :: x
    integer, intent(in) :: n_new
    character(len=*), intent(in), optional :: method
    character(len=1000) :: method_actual
    integer :: n, i
    double precision, dimension(n_new) :: y_new
    double precision, dimension(n_new) :: x_new
    n = size(y)
    ! If the old and new grid sizes are the same, do nothing
    if (n==n_new) then
      y_new = y
      return
    endif

    print *, 'oi', n, n_new

    if (.not.present(method)) then
      method_actual = 'polynomial'
    else
      method_actual = method
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
    call interpolate(x, y, x_new, y_new, method=trim(method_actual))

  end function rescale_array_old


  subroutine rescale_array(y, y_new, method)
    ! Interpolates the n-array y into an n_new array, with the same endpoints
    double precision, dimension(:), intent(in) :: y
    double precision, dimension(:), intent(out) :: y_new
    double precision, dimension(size(y)) :: x
    double precision, dimension(size(y_new)) :: x_new
    character(len=*), intent(in), optional :: method
    character(len=100) :: method_actual
    integer :: n, n_new, i
    n = size(y)
    n_new = size(y_new)

    ! If the old and new grid sizes are the same, do nothing
    if (n==n_new) then
      y_new = y
      return
    endif

    if (.not.present(method)) then
      method_actual = 'linear'
    else
      method_actual = method
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
    call interpolate(x, y, x_new, y_new, method=trim(method_actual))

  end subroutine rescale_array


  subroutine rescale_f_array(f, nx_new, nghost)
    ! Interpolates the f-array into the f_rescaled array
    ! NB the ghost zone of the new f-array has to be corrected afterwards
    double precision, dimension(:,:), allocatable, intent(inout) :: f
    double precision, dimension(:,:), allocatable :: f_old
    integer, intent(in) :: nx_new, nghost
    integer, dimension(2) :: shape_f
    integer :: nvar, nx, i

    shape_f = shape(f)
    nx   = shape_f(1)
    nvar = shape_f(2)

    ! Copies f -> tmp
    ! move_alloc avoids copying data unnecessarily however
    ! there seems to be a fortran bug... which is very sad
    !call move_alloc(f,f_old)
    allocate(f_old(nx, nvar))
    f_old = f
    ! Allocates new f-array
    deallocate(f)
    allocate(f(nx_new,nvar))

    do i=1,nvar
      ! Rescales, avoiding the ghost zone
      call rescale_array(f_old(nghost+1:,i), f(nghost+1:,i))
    enddo

  end subroutine rescale_f_array

end module interpolation
