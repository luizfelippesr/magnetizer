program testInterpolation
  use interpolation
  implicit none
  double precision, dimension(30,4) :: f
  double precision, dimension(13,4) :: f_old
  double precision, dimension(30) :: y
  double precision, dimension(20) :: y_new
  integer :: i
  double precision, dimension(100) :: xi1, yi1
  double precision :: x1(4) = [0.0D0, 0.1D0, 0.27D0, 0.3D0], &
       y1(4) = [0.15D0, 0.7D0, -0.1D0, 0.15D0]

  ! Test the interpolate subroutine (using the FGSL example)
  do i=1, 100
     xi1(i) = (1.D0 - dble(i-1)/100.D0) * x1(1) + dble(i-1)/100.D0 * x1(4)
  end do
  call interpolate(x1, y1, xi1, yi1, method='cubic_spline')
  do i=1,100
    print *, xi1(i), yi1(i)
  enddo
  print *, '#----------'

  ! Test rescaling an array to another grid
  y = [ 0.00060015, -0.00465305, -0.01957277, -0.04545156, -0.07912795, &
        -0.11393337, -0.14165434, -0.15459364, -0.14721517, -0.11715889,&
        -0.06559698,  0.00299473, 0.0815216, 0.16125866,  0.23298575,   &
        0.28812186,  0.31973156,  0.32330388,  0.29723662,  0.24299389, &
        0.16493829,  0.06987072, -0.03366414, -0.13622948, -0.22841741, &
        -0.3017493,  -0.3494712,  -0.36717296, -0.35317848, -0.30867717]
  y_new = rescale_array(y,size(y_new))

  print *, '#----------'

  do i=1,size(y_new)
    print *, y_new(i)
  enddo

  do i=1,3
    f(:,i) = y
  enddo

  print *, '#----------'
  call rescale_f_array(f,f_old)
  do i=1,13
    print *, f_old(i,2)
  enddo

end program testInterpolation
