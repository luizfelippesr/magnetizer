program testInterpolation
  ! Small program for testing the interpolation routines
  use interpolation
  implicit none
  double precision, dimension(30) :: x_ref,y_ref
  double precision, dimension(20) :: x_20,y_20
  double precision, dimension(40) :: x_40, y_40
  integer :: i
  character(len=100), parameter :: method ='akima'
  integer, parameter :: u= 17

  ! Produces a uniform grid from 0 to 100 (for each variable)
  x_ref = [ (100*dble(i-1)/dble(size(x_ref)-1), i=1,size(x_ref),1) ]
  x_20 = [ (100*dble(i-1)/dble(size(x_20)-1), i=1,size(x_20),1) ]
  x_40 = [ (100*dble(i-1)/dble(size(x_40)-1), i=1,size(x_40),1) ]

  ! Computes the reference values (something complicated...)
  y_ref = x_ref**2*sin(x_ref/10)/100+cos(x_ref/4.)*3

  call interpolate(x_ref, y_ref, x_40, y_40, method=method)
  call interpolate(x_ref, y_ref, x_20, y_20, method=method)

  open(unit=u, file='/tmp/interp_test_ref.dat')
  do i=1,size(x_ref)
      write(u, '(2(F9.2))') x_ref(i), y_ref(i)
  end do
  close(u)

  open(unit=u, file='/tmp/interp_test_40.dat')
  do i=1,size(x_40)
      write(u, '(2(F9.2))') x_40(i), y_40(i)
  end do
  close(u)

  open(unit=u, file='/tmp/interp_test_20.dat')
  do i=1,size(x_20)
      write(u, '(2(F9.2))') x_20(i), y_20(i)
  end do
  close(u)

end program testInterpolation
