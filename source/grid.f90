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
  double precision, parameter :: r_in=1.0d-7 ! min radius of disk in code units
                                             ! NB the max radius of disk is 1.0 in code units
                                             ! Should we make this zero?
  double precision, protected :: dx
  double precision, dimension(:), allocatable, target, protected :: x
  ! x and r have actually the exact same meaning
  ! one of them can be removed later
  double precision, dimension(:), pointer, protected :: r
  double precision, dimension(:), allocatable, protected :: r_kpc

  double precision, protected :: dr_kpc
  double precision, protected :: r_max_kpc

  double precision, private :: previous_r_max_kpc = -1

  contains

  subroutine construct_grid(r_disk, r_max_kpc_history)
    double precision, intent(in) :: r_disk, r_max_kpc_history
    integer :: i

    ! Initializes nxphys using the global input parameter value
    nxphys = p_nx_ref
    ! Includes the ghost zones
    nx= nxphys +2*nxghost
    ! First, prepares a grid from r_in to 1
    dx= (1d0 - r_in)/dble(nxphys-1)  !x corresponds to r coordinate
                                 !NB the max radius of disk is 1.0 in code units
    call check_allocate(x, nx)
    call check_allocate(r_kpc, nx)

    do i=1,nx
      x(i)=r_in -nxghost*dx +(dble(i)-1d0)*dx
    enddo

    ! For convenience. May be removed later
    r => x

    if (p_use_fixed_physical_grid) then
      r_max_kpc = r_max_kpc_history
    else
      r_max_kpc = p_rmax_over_rdisk*r_disk
    endif
    ! Save the dimensional r grid
    r_kpc = r*r_max_kpc
    dr_kpc = dx*r_max_kpc
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

  subroutine check_allocate_f_array(f, nvar)
    ! Helper subroutine that checks whether the f-array has the correct format,
    ! allocating it when necessary
    double precision, dimension(:,:), allocatable, intent(inout) :: f
    integer, intent(in) :: nvar
    integer, dimension(2) :: shape_f

    if (allocated(f)) then
      shape_f = shape(f)
      if (shape_f(1)==nx) then
        ! If the v was allocated and has the correct shape
        ! then there is nothing to do..
        return
      else
        ! If the shape is wrong, deallocates
        deallocate(f)
      endif
    endif
    ! Allocates f with the correct shape
    allocate(f(nx,nvar))
    f = 0
    return
  end subroutine check_allocate_f_array

  subroutine adjust_grid(f, dfdt, r_disk)
    ! Rescales the grid acording to the global input parameters
    ! The grid can be extended, keeping the resolution dr_kpc or the grid
    ! can be kept, and the r_kpc changed to accomodate the new disk size.
    ! The f-array is updated accordingly.
    use messages
    use interpolation
    double precision, dimension(:,:), allocatable, intent(inout) :: f, dfdt
    double precision, intent(in) :: r_disk
    double precision, dimension(:), allocatable :: tmp
    double precision, dimension(:,:), allocatable :: ftmp
    double precision :: r_max_candidate, relative_variation
    integer :: nx_new, i, nvar
    integer, dimension(2) :: shape_f
    double precision, parameter :: GRID_ADJ_THRES = 0.05

    if (p_use_fixed_physical_grid) then
      ! Do nothing!
      return
    endif

    shape_f = shape(f)
    nvar = shape_f(2)

    if (p_scale_back_f_array) then
      ! Restores the grid to the default resolution
      nx = p_nx_ref
      ! Rescales the f_array
      call rescale_f_array(f, nx)
      ! Rescales the dimensional grid
      call rescale_array_in_place(r_kpc, nx)
      ! Rescales the dimensionless grid
      call check_allocate(x, nx)
      x = r_kpc/r_max_kpc
      nullify(r); r => x
      ! Adjusts the shape of dfdt
      call check_allocate_f_array(dfdt, nx)
      ! Updates the resolution
      dr_kpc = r_kpc(2)-r_kpc(1)
    endif

    ! Saves previous state and computes targed r_max_kpc
    previous_r_max_kpc = r_max_kpc
    r_max_kpc = p_rmax_over_rdisk*r_disk

    ! Does nothing ignore if the change in disk size is negligible
    relative_variation = abs(r_max_kpc-previous_r_max_kpc)/previous_r_max_kpc
    if (relative_variation < GRID_ADJ_THRES) then
      ! Restors r_max_kpc and exits
      r_max_kpc = previous_r_max_kpc
      return
    endif

    if (r_max_kpc >= previous_r_max_kpc) then
      if (p_rescale_field_for_expanding_disks) then
        ! Keeps the f-array untouched
        ! Rescales the dimensional r grid
        r_kpc = r_kpc * r_max_kpc/previous_r_max_kpc
        dr_kpc = dr_kpc * r_max_kpc/previous_r_max_kpc
      else
        nx_new = nx
        r_max_candidate = previous_r_max_kpc
        do while (r_max_candidate < r_max_kpc)
          ! Keeps the resolution and increments 1 grid point
          nx_new = nx_new + 1
          r_max_candidate = (r_in/dx -nxghost +(dble(nx_new-nxghost)-1.))*dr_kpc
        enddo
        ! Fails if grid is too big
        if (nx_new > p_nx_MAX) then
          call message('adjust_grid: nx =', val_int=nx_new, &
                       msg_end='is greater than the limit set in p_nx_MAX.')
          stop
        endif

        ! Some log messages
        call message('adjust_grid: Temporarily rescaling the grid to nx=', &
                     val_int=nx_new, info=2)
        call message('adjust_grid: Previous r_max_kpc=', info=2, &
                     val=previous_r_max_kpc)
        call message('adjust_grid: New r_max_kpc=', info=2, &
                     val=r_max_candidate)

        ! Saves new maximum size
        r_max_kpc = r_max_candidate

        ! Adjusts the shape of the grid itself to the new nx:
        ! - First makes a copy of the dimensional grid
        call move_alloc(r_kpc, tmp)
        ! - Allocates the new dimensional grid
        call check_allocate(r_kpc, nv=nx_new)
        ! - Adds previous grid data
        r_kpc(:nx)=tmp
        deallocate(tmp)
        ! - Extends the grid
        do i=nx,nx_new
          r_kpc(i) = r_kpc(i-1)+dr_kpc
        enddo
        ! - Adjusts and sets the dimensionless grid
        call check_allocate(x, nv=nx_new)
        x = r_kpc/r_max_kpc
        nullify(r); r => x
        ! - Updates dx
        dx = x(2)-x(1)

        ! Prepares allocates new f-array
        call move_alloc(f, ftmp)
        allocate(f(nx_new,nvar))
        ! Moves the content from the previous to the new f-array (avoiding the
        ! boundary)
        f(:nx-nxghost-1,:) = ftmp(:nx-nxghost-1,:)
        ! Fills up the new parts
        do i=nx-nxghost,nx_new
          ! The value in the point before the boundary is used, scaled by 1/r
          f(i,:) = ftmp(nx-nxghost-1,:) * r_kpc(nx-nxghost-1)/r_kpc(i)
        enddo
        deallocate(ftmp)

        ! Checks the allocation of the dfdt array
        call check_allocate_f_array(dfdt, nx_new)

        ! Saves the new number of grid points
        nx = nx_new
      endif
    else
      if (p_rescale_field_for_shrinking_disks) then
        ! Keeps the f-array untouched
        ! Updates the dimensional r grid
        r_kpc = r_kpc * r_max_kpc/previous_r_max_kpc
        dr_kpc = dr_kpc * r_max_kpc/previous_r_max_kpc
      else
        call message('adjust_grid: p_rescale_field_for_shrinking_disks=False'&
                     //', option not yet implemented')
        stop
      endif
    endif

  end subroutine adjust_grid

end module grid
