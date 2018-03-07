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
    double precision, allocatable, dimension(:,:) :: value
  end type ts_array

  type, public :: ts_data
    ! Defines the ts_data class
    type(ts_array), dimension(:), allocatable :: data
  contains
    procedure :: reset => reset_ts_arrays
    procedure :: get => get_ts_values
    procedure :: get_ts => get_ts_array
    procedure :: is_scalar => is_ts_scalar
!     procedure :: set => set_ts_values
  end type ts_data

  interface ts_data
    module procedure new_ts_data
  end interface

contains
  function new_ts_data(globals, profiles, nz, nx)
    implicit none
    character(len=NAMESIZE), dimension(:), intent(in) :: globals, profiles
    type(ts_data) new_ts_data
    integer, intent(in) :: nz, nx
    integer :: i, j, nglobals, nprofiles

    nglobals = size(globals)
    nprofiles = size(profiles)

    allocate(new_ts_data%data(nglobals+nprofiles))
    do i=1,nglobals
      new_ts_data%data(i)%scalar = .true.
      new_ts_data%data(i)%name = globals(i)
      allocate(new_ts_data%data(i)%value(nz,1))
      new_ts_data%data(i)%value = INVALID
    enddo

    do i=1,nprofiles
      j = i + nglobals
      new_ts_data%data(j)%scalar = .false.
      new_ts_data%data(j)%name = profiles(i)
      allocate(new_ts_data%data(j)%value(nz, nx))
      new_ts_data%data(j)%value = INVALID
    enddo

  end function new_ts_data

  function get_ts_array(self, name) result(ts)
    class(ts_data), intent(in) :: self
    type(ts_array) :: ts
    character(len=*), intent(in) :: name
    integer :: i

    do i=1, size(self%data)
      if (trim(name) == trim(self%data(i)%name)) then
        ts = self%data(i)
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

program test_tsDataObj
  use tsDataObj
  implicit none
    integer, parameter :: NAMESIZE = 15
  character(len=NAMESIZE), dimension(3) :: profiles
  character(len=NAMESIZE), dimension(2) :: globals
  type(ts_data) :: ts_teste
  type(ts_array) :: ts
  double precision, dimension(3,2) :: teste
  double precision, dimension(3,1) :: testes

  globals = ['a              ','b              ']
  profiles = ['    AAAAA      ','    BBBBBB     ','teste          ']
  ts_teste = ts_data(globals, profiles,3,2)
!   call ts_teste%initialize(3,6)
!   call ts_teste%reset
!   call ts_teste%initialize(3,2)

!   ts = ts_teste%get('scalar')
!   print *, ts%name
!   print *, ts%scalar
!   print *, ts%value
!   print *, shape(ts%value)
!   print *, '--------------------'
  testes = ts_teste%get('a')
  print *, 'a'
  print *, testes
  print *, shape(testes)
  print *, ts_teste%is_scalar('a')
  print *, '--------------------'
  teste = ts_teste%get('teste')
  print *, 'teste'
  print *, teste
  print *, shape(teste)
  print *, ts_teste%is_scalar('a')
  print *, '--------------------'


end program test_tsDataObj