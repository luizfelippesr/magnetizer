! Contains grid parameters and subroutine for constructing grid
module grid
  use math_constants
  use units
  use global_input_parameters

  implicit none

  integer :: nxphys= -1 !Resolution in r (excluding ghost zones) (for convenience should be N*r_max_kpc+1)
  integer, parameter :: nxghost= 3 !Number of ghost cells at each end in r
  integer, protected :: nx=-1  !Resolution in r
  double precision, parameter :: len_phi= 2*pi !phi domain
  double precision, parameter :: r_in=1.0d-6 ! min radius of disk in code units
                                             ! NB the max radius of disk is 1.0 in code units
  double precision, protected :: dx
  double precision, dimension(:), allocatable, protected :: x
  ! x and r have actually the exact same meaning
  ! one of them can be removed later
  double precision, dimension(:), allocatable, protected :: r
  double precision, dimension(:), allocatable, protected :: r_kpc
  double precision, protected :: dr_kpc

  contains


  subroutine construct_grid(r_max_kpc)
    double precision, intent(in) :: r_max_kpc
    integer :: i

    ! Initializes nxphys using the global input parameter value
    nxphys = p_nx_ref
    ! Includes the ghost zones
    nx= nxphys +2*nxghost
    ! First, prepares a grid from r_in to 1
    dx= (1.0 - r_in)/(nxphys-1)  !x corresponds to r coordinate
                                 !NB the max radius of disk is 1.0 in code units
    call check_allocate(x, nx)
    call check_allocate(r, nx)
    call check_allocate(r_kpc, nx)

    do i=1,nx
      x(i)=r_in -nxghost*dx +(dble(i)-1.)*dx
    enddo
    r=x !use explicit name

    ! Save the dimensional r grid
    r_kpc = r*r_max_kpc

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

  subroutine adjust_grid(f, previous_f, r_disk, previous_r_disk)
    double precision, dimension(:,:), intent(inout) :: f
    double precision, dimension(:,:), intent(inout) :: previous_f
    double precision, intent(in) :: r_disk
    double precision, intent(in) :: previous_r_disk

    if (p_use_fixed_physical_grid) then
      f = previous_f
    endif
!     else if (r_disk > previous_r_disk) then
!       if (p_rescale_field_for_expanding_disks) then
!         f = previous_f
!       else`
!     p_rescale_field_for_shrinking_disks,
!     p_rescale_field_for_expanding_disks,
!     p_nx_ref,
!     p_nx_MAX,
!
!

  end subroutine adjust_grid

end module grid
