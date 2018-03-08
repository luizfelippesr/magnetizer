module tsDataObj
  implicit none
  private

  double precision, parameter :: INVALID = -99999d0 ! Invalid value tag
  integer, parameter :: NAMESIZE = 15 ! Maximum name size

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
    procedure :: get_scalar => get_ts_values_scalar
    procedure :: get_it => get_ts_single_value
    procedure :: get_ts => get_ts_array
    procedure :: is_scalar => is_ts_scalar
    procedure :: output => is_ts_in_output
    procedure :: set_full => set_ts_full
    procedure :: set => set_ts_single
    procedure :: set_scalar => set_ts_single_scalar
  end type tsData

  interface tsData
    module procedure new_tsData
  end interface

contains
  function new_tsData(scalars, profiles, nz, nx, output)
    implicit none
    character(len=*), dimension(:), intent(in) :: scalars, profiles
    character(len=*), dimension(:), intent(in), optional :: output
    type(tsData), target :: new_tsData
    integer, intent(in) :: nz, nx
    integer :: i, j, nscalars, nprofiles

    nscalars = size(scalars)
    nprofiles = size(profiles)

    allocate(new_tsData%data(nscalars+nprofiles))
    do i=1,nscalars
      new_tsData%data(i)%scalar = .true.
      if (present(output)) new_tsData%data(i)%output = is_in_list(scalars(i),output)
      new_tsData%data(i)%name = scalars(i)
      allocate(new_tsData%data(i)%value(nz,1))
      new_tsData%data(i)%value = INVALID
    enddo

    do i=1,nprofiles
      j = i + nscalars
      new_tsData%data(j)%scalar = .false.
      new_tsData%data(j)%name = profiles(i)
      if (present(output)) new_tsData%data(j)%output = is_in_list(profiles(i),output)
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
    ! If nothing is found, initializes the pointer to null
    nullify(ts)
  end function get_ts_array


  function get_ts_values(self, name) result(values)
    class(tsData), intent(in) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:,:), allocatable :: values
    type(ts_array), pointer :: ts
    integer, dimension(2) :: dshape

    ts => self%get_ts(name)
    dshape = shape(ts%value)
    allocate(values(dshape(1),dshape(2)))
    values = ts%value
  end function get_ts_values

  function get_ts_values_scalar(self, name) result(values)
    class(tsData), intent(in) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:), allocatable :: values
    type(ts_array), pointer :: ts
    integer, dimension(2) :: dshape

    ts => self%get_ts(name)
    dshape = shape(ts%value)
    allocate(values(dshape(1)))
    values = ts%value(:,1)
  end function get_ts_values_scalar

  function get_ts_single_value(self, name, it) result(values)
    class(tsData), intent(in) :: self
    character(len=*), intent(in) :: name
    type(ts_array), pointer :: ts
    double precision, dimension(:), allocatable :: values
    integer :: it
    integer, dimension(2) :: dshape


    ts => self%get_ts(name)
    dshape = shape(ts%value)
    allocate(values(dshape(2)))
    values = ts%value(it,:)
  end function get_ts_single_value

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
    type(ts_array), pointer :: ts

    ts => self%get_ts(name)
    is_ts_scalar = ts%scalar
  end function is_ts_scalar

  logical function is_ts_in_output(self, name)
    class(tsData), intent(in) :: self
    character(len=*), intent(in) :: name
    type(ts_array), pointer :: ts

    ts => self%get_ts(name)
    is_ts_in_output = ts%output
  end function is_ts_in_output

  subroutine reset_ts_arrays(self)
    implicit none
    class(tsData), intent(inout) :: self
    integer :: i
    ! Deallocates the time series arrays
    do i=1, size(self%data)
      self%data(i)%value = INVALID
    end do
  end subroutine reset_ts_arrays

  logical function is_in_list(quantity, list)
    character(len=*), dimension(:), intent(in) :: list
    character(len=*), intent(in) :: quantity
    integer :: i

    do i=1, size(list)
      if (trim(list(i))==trim(quantity)) then
        is_in_list = .true.
        return
      endif
    enddo
    is_in_list = .false.
  end function is_in_list
end module tsDataObj
