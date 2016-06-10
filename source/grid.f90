! Contains grid parameters and subroutine for constructing grid
module grid
  use math_constants
  use units
!
  implicit none
!
  integer, parameter :: nxphys= 51!201!61!91  !Resolution in r (excluding ghost zones) (for convenience should be N*r_max_kpc+1)
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in r
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in r
  double precision, parameter :: len_phi= 2*pi !phi domain
  double precision, parameter :: r_in=1.0d-6 ! min radius of disk
                                             ! NB the max radius of disk is 1.0 in code units
  double precision :: dx
  double precision, dimension(nx) :: x
  double precision, dimension(nx) :: r
  double precision, dimension(nx) :: r_kpc
!
  contains


  subroutine construct_grid

    integer :: i

      dx= (1.0 - r_in)/(nxphys-1)  !x corresponds to r coordinate
                                   !NB the max radius of disk is 1.0 in code units
      do i=1,nx
        x(i)=r_in -nxghost*dx +(dble(i)-1.)*dx
      enddo
      r=x !use explicit name
    endsubroutine construct_grid

  subroutine check_allocate(v, nv)
    ! Helper subroutine that checks whether 1d-array has the correct format,
    ! allocating it when necessary
    ! v -> the array that will be checked and possibly reallocated
    ! otpional, nv -> the expected length of the array
    !                 (if absent, use current value of nx)
    double precision, dimension(:), allocatable, intent(inout) :: v
    integer, intent(in), optional :: nv
    integer :: nv_actual

    if (present(nv)) then
      nv_actual = nv
    else
      nv_actual = nx
    endif

    if (allocated(v)) then
      if (size(v)==nv_actual) then
        ! If the v was allocated and has the correct shape
        ! then there is nothing to do..
        return
      else
        ! If the shape is wrong, deallocates
        deallocate(v)
      endif
    endif
    ! Allocates v with the correct shape
    allocate(v(nv_actual))
    return
  end subroutine check_allocate

end module grid
