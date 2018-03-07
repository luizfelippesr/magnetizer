module tsDataObj
  implicit none
  private

!   double precision, parameter :: INVALID = -99999d0
  double precision, parameter :: INVALID = -17d0
  integer, parameter :: NAMESIZE = 15

  type, public :: ts_array
    ! Package containing a time series array, the quantity name and
    ! the whether it should be in the output or not.
    character(len=NAMESIZE) :: name
    logical :: output
    logical :: scalar
    double precision, pointer, dimension(:,:) :: value
  end type ts_array

  type, public :: ts_data
    ! Defines the ts_data class
    type(ts_array), dimension(:), pointer :: data
    logical :: initialized
  contains
    procedure :: initialize => initialize_ts_arrays
    procedure :: reset => reset_ts_arrays
    procedure :: get => get_ts_values
    procedure :: get_ts => get_ts_array
    procedure :: is_scalar => is_ts_scalar
    procedure :: set => set_ts_values
  end type ts_data

  interface ts_data
    module procedure new_ts_data
  end interface

contains
  function new_ts_data(globals, profiles)
    implicit none
    character(len=NAMESIZE), dimension(:), intent(in) :: globals, profiles
    type(ts_data), target :: new_ts_data
    integer :: i, j, nglobals, nprofiles

    new_ts_data%initialized = .false.
    nglobals = size(globals)
    nprofiles = size(profiles)

    allocate(new_ts_data%data(nglobals+nprofiles))
    do i=1,nglobals
      new_ts_data%data(i)%scalar = .true.
      new_ts_data%data(i)%name = globals(i)
    enddo

    do i=1,nprofiles
      j = i + nglobals
      new_ts_data%data(j)%scalar = .false.
      new_ts_data%data(j)%name = profiles(i)
    enddo
  end function new_ts_data


  subroutine initialize_ts_arrays(self,nz, nx)
    implicit none
    class(ts_data), intent(inout) :: self
    integer, intent(in) :: nz, nx
    integer :: i

    do i=1,size(self%data)
      if (self%data(i)%scalar) then
        allocate(self%data(i)%value(nz,1))
      else
        allocate(self%data(i)%value(nz, nx))
      endif
      self%data(i)%value = INVALID
    enddo
    self%initialized = .true.
  end subroutine initialize_ts_arrays


  function get_ts_array(self, name) result(ts)
    class(ts_data), intent(in) :: self
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
    class(ts_data), intent(in) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:,:), allocatable :: values
    type(ts_array) :: ts
    integer, dimension(2) :: dshape

    ts = self%get_ts(name)
    dshape = shape(ts%value)
    allocate(values(dshape(1),dshape(2)))
    values = ts%value
  end function get_ts_values


  subroutine set_ts_values(self, name, vals)
    class(ts_data), intent(inout) :: self
    character(len=*), intent(in) :: name
    double precision, dimension(:,:) :: vals
    type(ts_array), pointer :: ts

    ts => self%get_ts(name)
    ts%value = vals
  end subroutine set_ts_values


  logical function is_ts_scalar(self, name)
    class(ts_data), intent(in) :: self
    character(len=*), intent(in) :: name
    type(ts_array) :: ts

    ts = self%get_ts(name)
    is_ts_scalar = ts%scalar
  end function is_ts_scalar

  subroutine reset_ts_arrays(self)
    implicit none
    class(ts_data), intent(inout) :: self
    integer :: i
    ! Deallocates the time series arrays
    do i=1, size(self%data)
      self%data(i)%value = INVALID
    end do
  end subroutine reset_ts_arrays

end module tsDataObj
