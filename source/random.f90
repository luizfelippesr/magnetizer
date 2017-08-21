module random
! Subroutines and functions related to random number generation
  implicit none
  private

  public :: set_random_seed, random_normal

contains

  subroutine set_random_seed(igal, p_random_seed, use_array)
    ! Initializes the random seed based on igal and p_random_seed
    ! (actually, any two integers!)
    ! The optional argument use_array uses a more complex seed
    ! (possibly making the random sequence more random) but
    ! a the cost of making the RNG platform dependent.
    integer, intent(in) :: igal, p_random_seed
    logical, optional, intent(in) :: use_array
    integer, allocatable, dimension(:) :: seed
    integer :: i, n
    double precision :: p

    ! Discovers the size of the seed array and allocates it
    call random_seed(size = n)
    allocate(seed(n))

    ! Does *arbitrary* calculation to combine p_random_seed parameter and igal
    seed = p_random_seed**2 - igal**3

    ! Ideally, the seed should be an array. *However*, n may vary from
    ! platform to platform! Thus, this is inactive by _default_.
    if (present(use_array)) then
      if (use_array) then
        do i=1, n
          ! Arbitrary calculation with array index
          seed(i) = seed(i) + i*igal - i*p_random_seed
        enddo
      endif
    endif

    ! Zeros values for the seed lead to poor random number generation
    ! (according to gfortran documentation)
    where (seed==0)
      ! Just substitute it for anything else...
      seed = p_random_seed**igal
    end where

    ! Changes the seed!
    call random_seed(put=seed)
  end subroutine set_random_seed


  double precision function random_normal()
    ! Adapted from the code by Alan Miller (alan @ mel.dms.csiro.au) obtained
    ! at http://wp.csiro.au/alanmiller/random.html
    !
    ! His code was itself based on the Fortran 77 code
    !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
    !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
    !
    !  The function random_normal() returns a normally distributed pseudo-random
    !  number with zero mean and unit variance.   This version uses the default
    !  uniform random number generator which is in your fortran library.
    !
    !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
    !  and J.F. Monahan augmented with quadratic bounding curves.
    double precision, parameter :: s = 0.449871, t = -0.386595
    double precision, parameter :: a = 0.19600, b = 0.25472
    double precision, parameter :: half = 0.5, r1 = 0.27597, r2 = 0.27846
    double precision :: u, v, x, y, q

    ! Generate p = (u,v) uniform in rectangle enclosing acceptance region
    do
      call random_number(u)
      call random_number(v)
      v = 1.7156 * (v - half)

      ! Evaluate the quadratic form
      x = u - s
      y = abs(v) - t
      q = x**2 + y*(a*y - b*x)

      ! Accept p if inside inner ellipse
      if (q < r1) exit
      ! Reject p if outside outer ellipse
      if (q > r2) cycle
      ! Reject p if outside acceptance region
      if (v**2 < -4.0*log(u)*u**2) exit
    end do

    ! Return ratio of P's coordinates as the normal deviate
    random_normal = v/u
    return
  end function random_normal

end module random
