program testRoot
  use root_finder
  implicit none
  double precision :: x, a3, a2, a1, a0, guess

  a3 = 1d0
  a2 = -1d0
  a1 = -2d0
  a0 = -1d0
  guess = 20

  x = CubicRootClose(a3, a2, a1, a0, guess)
  print *, '(a3,a2,a1,a0)=',a3, a2, a1, a0
  print *, 'Found x=',x
  print *, 'Value at solution:', a3*x**3+a2*x**2+a1*x+a0

  a3 = 99.7d0
  a2 = 0.78787878d0
  a1 = 45d1
  a0 = -17d2
  guess = 10000

  x = CubicRootClose(a3, a2, a1, a0, guess)
  print *, '(a3,a2,a1,a0)=',a3, a2, a1, a0
  print *, 'Found x=',x
  print *, 'Value at solution:', a3*x**3+a2*x**2+a1*x+a0

end program testRoot
