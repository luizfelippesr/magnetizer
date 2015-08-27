! Contains grid parameters and subroutine for constructing grid
module grid
  use math_constants
  use units
!
  implicit none
!
  integer, parameter :: nxphys= 201!61!91  !Resolution in r (excluding ghost zones) (for convenience should be N*r_disk_kpc+1)
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in r
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in r
  double precision, parameter :: len_phi= 2*pi, r_in=0.0d0  !phi domain, min radius of disk, max radius of disk
  double precision :: dx
  double precision, dimension(nx) :: x
  double precision, dimension(nx) :: r
!
  contains
    subroutine construct_grid

    integer :: i

      dx= (r_disk-r_in)/(nxphys-1)  !x corresponds to r coordinate
      do i=1,nx
        x(i)=r_in -nxghost*dx +(dble(i)-1.)*dx
      enddo
      r=x !use explicit name
    endsubroutine construct_grid
end module grid
