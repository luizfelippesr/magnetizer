! Contains grid parameters and subroutine for constructing grid
module grid
  use math_constants
  use units
!
  implicit none
!
  integer, parameter :: nxphys= 201!61!91  !Resolution in r (excluding ghost zones) (for convenience should be N*r_max_kpc+1)
  integer, parameter :: nxghost= 3  !Number of ghost cells at each end in r
  integer, parameter :: nx= nxphys +2*nxghost  !Resolution in r
  double precision, parameter :: len_phi= 2*pi !phi domain
  double precision, parameter :: r_in=1.0d-6 ! min radius of disk
                                             ! NB the max radius of disk is 1.0 in code units
  double precision :: dx
  double precision, dimension(nx) :: x
  double precision, dimension(nx) :: r
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
end module grid
