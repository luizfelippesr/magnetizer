module tsDataObj
  implicit none
  private

  double precision, parameter :: INVALID = -99999d0
  integer, parameter :: NAMESIZE = 15

  type, public :: ts_array
    ! Package containing a time series array, the quantity name and
    ! the whether it should be in the output or not.
    character(len=NAMESIZE) :: name
    logical :: output
    logical :: scalar
    double precision, pointer, dimension(:,:) :: value
  end type ts_array

  type, public :: tsData
    ! Defines the tsData class
    type(ts_array), dimension(:), pointer :: data
    logical :: initialized = .false.
  contains
    procedure :: reset => reset_ts_arrays
    procedure :: get => get_ts_values
    procedure :: get_ts => get_ts_array
    procedure :: is_scalar => is_ts_scalar
    procedure :: set => set_ts_full
    procedure :: set_value => set_ts_single
    procedure :: set_value_scalar => set_ts_single_scalar
  end type tsData

  interface tsData
    module procedure new_tsData
  end interface

contains
  function new_tsData(scalars, profiles, nz, nx)
    implicit none
    character(len=NAMESIZE), dimension(:), intent(in) :: scalars, profiles
    type(tsData), target :: new_tsData
    integer, intent(in) :: nz, nx
    integer :: i, j, nscalars, nprofiles

    nscalars = size(scalars)
    nprofiles = size(profiles)

    allocate(new_tsData%data(nscalars+nprofiles))
    do i=1,nscalars
      new_tsData%data(i)%scalar = .true.
      new_tsData%data(i)%name = scalars(i)
      allocate(new_tsData%data(i)%value(nz,1))
      new_tsData%data(i)%value = INVALID
    enddo

    do i=1,nprofiles
      j = i + nscalars
      new_tsData%data(j)%scalar = .false.
      new_tsData%data(j)%name = profiles(i)
      allocate(new_tsData%data(j)%value(nz, nx))
      new_tsData%data(j)%value = INVALID
    enddo

    new_tsData%initialized = .true.
  end function new_tsData

  function get_ts_array(self, name) result(ts)
    class(tsData), intent(in) :: self
    type(ts_array), pointer :: ts
    character(len=*), intent(in) :: name
    integer :: i

    do i=1, size(self%data)
      if (trim(name) == trim(self%data(i)%name)) then
        ts => self%data(i)
        return
      endif
    enddo
  end function get_ts_array


  function get_ts_values(self, name) result(values)
    class(tsData), intent(in) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:,:), allocatable :: values
    type(ts_array) :: ts
    integer, dimension(2) :: dshape

    ts = self%get_ts(name)
    dshape = shape(ts%value)
    allocate(values(dshape(1),dshape(2)))
    values = ts%value
  end function get_ts_values


  subroutine set_ts_full(self, name, vals)
    class(tsData), intent(inout) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:,:) :: vals
    type(ts_array), pointer :: ts

    ts => self%get_ts(name)
    ts%value = vals
  end subroutine set_ts_full


  subroutine set_ts_single(self, name, it, vals)
    class(tsData), intent(inout) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:) :: vals
    type(ts_array), pointer :: ts
    integer :: it

    ts => self%get_ts(name)
    ts%value(it,:) = vals(:)
  end subroutine set_ts_single

  subroutine set_ts_single_scalar(self, name, it, vals)
    class(tsData), intent(inout) :: self
    character(len=*), intent(in) :: name
    double precision :: vals
    type(ts_array), pointer :: ts
    integer :: it

    ts => self%get_ts(name)
    ts%value(it,1) = vals
  end subroutine set_ts_single_scalar


  logical function is_ts_scalar(self, name)
    class(tsData), intent(in) :: self
    character(len=*), intent(in) :: name
    type(ts_array) :: ts

    ts = self%get_ts(name)
    is_ts_scalar = ts%scalar
  end function is_ts_scalar

  subroutine reset_ts_arrays(self)
    implicit none
    class(tsData), intent(inout) :: self
    integer :: i
    ! Deallocates the time series arrays
    do i=1, size(self%data)
      self%data(i)%value = INVALID
    end do
  end subroutine reset_ts_arrays

end module tsDataObj
