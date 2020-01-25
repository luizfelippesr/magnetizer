!# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
!#
!# This file is part of Magnetizer.
!#
!# Magnetizer is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# Magnetizer is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
!#
module data_transfer
  use mpi
  use IO
  implicit none
  private

  ! Maximum sizes
  integer, parameter :: NAMESIZE = 15
  integer, parameter :: DESCRIPTIONSIZE = 100
  integer, parameter :: UNITSSIZE = 100
  integer, parameter :: master_rank = 0
  interface send_or_write_dataset
    module procedure send_dataset_scalar
    module procedure send_dataset_vector
    module procedure send_dataset_code
  end interface send_or_write_dataset

  public send_or_write_dataset, send_end_message, receive_and_save

  contains


    subroutine check_allocate(a, a_shape)
      double precision, dimension(:,:), allocatable, intent(inout) :: a
      integer, intent(in), dimension(2) :: a_shape
      integer, dimension(2) :: present_shape
      if (allocated(a)) then
        present_shape = shape(a)
        if (present_shape(1) == a_shape(1) .and. &
            present_shape(2) == a_shape(2)) then
          ! If the a was allocated and has the correct shape
          ! then there is nothing to do..
          return
        else
          ! If the shape is wrong, deallocates
          deallocate(a)
        endif
      endif
      ! Allocates a with the correct shape
      allocate(a(a_shape(1),a_shape(2)))
      return
  end subroutine check_allocate

    subroutine check_allocate_code(a, length)
      character, dimension(:), allocatable, intent(inout) :: a
      integer, intent(in) :: length
      if (allocated(a)) then
        if (len(a) == length) then
          ! If the a was allocated and has the correct shape
          ! then there is nothing to do..
          return
        else
          ! If the shape is wrong, deallocates
          deallocate(a)
        endif
      endif
      ! Allocates a with the correct shape
      allocate(a(length))
      return
  end subroutine check_allocate_code


  subroutine receive_and_save(igal, iproc)
    integer, intent(in) :: igal, iproc
    character(len=DESCRIPTIONSIZE) :: description
    character(len=NAMESIZE) :: dataset_name, group
    character(len=UNITSSIZE) :: units
    double precision, dimension(:,:), allocatable :: data
    character, dimension(:), allocatable :: data_code
    integer, dimension(2) :: data_shape
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: ierr

    do
      call MPI_Recv(dataset_name, NAMESIZE, MPI_CHARACTER, iproc, 1, &
                    MPI_COMM_WORLD, status, ierr)
      if (trim(dataset_name) == 'done') exit

      call MPI_Recv(description, DESCRIPTIONSIZE, MPI_CHARACTER, iproc, 2, &
                    MPI_COMM_WORLD, status, ierr)

      call MPI_Recv(units, UNITSSIZE, MPI_CHARACTER, iproc, 3, &
                    MPI_COMM_WORLD, status, ierr)
      call MPI_Recv(group, NAMESIZE, MPI_CHARACTER, iproc, 4, &
               MPI_COMM_WORLD, status, ierr)
      call MPI_Recv(data_shape, 2, MPI_INTEGER, iproc, 5, &
                MPI_COMM_WORLD, status, ierr)

      if (data_shape(1)==0) then
        ! data_shape(1)=0 is a marker, indicating that the data is a character
        ! (i.e. error code) array
        call check_allocate_code(data_code, data_shape(2))
        call MPI_Recv(data_code, data_shape(2), MPI_CHARACTER, &
                 iproc, 6, MPI_COMM_WORLD, status, ierr)
        call IO_write_dataset(trim(dataset_name), igal, data_code, &
                              units=trim(units), description=trim(description),&
                              group=trim(group))
      else if (data_shape(1)==1) then
        ! data_shape(1)=1 indicates scalar quantity
        call check_allocate(data, data_shape)
        call MPI_Recv(data(1,:), data_shape(2), MPI_DOUBLE_PRECISION, &
                 iproc, 6, MPI_COMM_WORLD, status, ierr)
        call IO_write_dataset(trim(dataset_name), igal, data(1,:), &
                              units=trim(units), description=trim(description),&
                              group=trim(group))
      else
        call check_allocate(data, data_shape)
        call MPI_Recv(data, data_shape(1)*data_shape(2), MPI_DOUBLE_PRECISION, &
                 iproc, 6, MPI_COMM_WORLD, status, ierr)
        call IO_write_dataset(trim(dataset_name), igal, data, &
                              units=trim(units), description=trim(description))
      endif
    enddo
  end subroutine receive_and_save

  subroutine send_end_message
    character(len=NAMESIZE) :: name_a
    integer :: ierr
    if (mpirank==0) return

    name_a = 'done'
    call MPI_Send(name_a, NAMESIZE, MPI_CHARACTER, master_rank, 1, &
                  MPI_COMM_WORLD, ierr)
  end subroutine send_end_message

  subroutine send_dataset_vector(dataset_name, gal_id, data, &
                                 units, description, group)
    ! Writes a dataset to disk - vector version
    character(len=*), intent(in) :: dataset_name
    character(len=*), optional, intent(in) :: group
    integer, intent(in) :: gal_id
    double precision, dimension(:,:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    character(len=DESCRIPTIONSIZE) :: description_a
    character(len=NAMESIZE) :: name_a, group_a
    character(len=UNITSSIZE) :: units_a
    integer, dimension(2) :: data_shape
    integer :: ierr

    group_a = 'Output'; description_a = ''; units_a = ''
    if (present(group))  group_a = group
    if (present(description))  description_a = description
    if (present(units))  units_a = units
    name_a = dataset_name
    ! If this is the master, simply write things to disk
    if (mpirank==0) then
      call IO_write_dataset(dataset_name, gal_id, data, units_a, description_a)
      return
    endif

    call MPI_Send(name_a, NAMESIZE, MPI_CHARACTER, master_rank, 1, &
                  MPI_COMM_WORLD, ierr)

    call MPI_Send(description_a, DESCRIPTIONSIZE, MPI_CHARACTER, master_rank, 2, &
                  MPI_COMM_WORLD, ierr)

    call MPI_Send(units_a, UNITSSIZE, MPI_CHARACTER, master_rank, 3, &
                  MPI_COMM_WORLD, ierr)

    call MPI_Send(group_a, NAMESIZE, MPI_CHARACTER, master_rank, 4, &
                  MPI_COMM_WORLD, ierr)

    data_shape = shape(data)
    call MPI_Send(data_shape, 2, MPI_INTEGER, master_rank, 5, &
             MPI_COMM_WORLD, ierr)

    call MPI_Send(data, data_shape(1)*data_shape(2), MPI_DOUBLE_PRECISION, &
             master_rank, 6, MPI_COMM_WORLD, ierr)

    return
  end subroutine send_dataset_vector

  subroutine send_dataset_code(dataset_name, gal_id, data, &
                                     units, description, group)
    ! Writes a dataset to disk - vector version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    character(len=*), optional, intent(in) :: group
    character, dimension(:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    character(len=DESCRIPTIONSIZE) :: description_a
    character(len=NAMESIZE) :: name_a, group_a
    character(len=UNITSSIZE) :: units_a
    integer, dimension(2) :: data_shape
    integer :: ierr

    group_a = 'Output'; description_a = ''; units_a = ''
    if (present(group))  group_a = group
    if (present(description))  description_a = description
    if (present(units))  units_a = units
    name_a = dataset_name
    ! If this is the master, simply write things to disk
    if (mpirank==0) then
      call IO_write_dataset(dataset_name, gal_id, data, units_a, description_a, group_a)
      return
    endif

    call MPI_Send(name_a, NAMESIZE, MPI_CHARACTER, master_rank, 1, &
             MPI_COMM_WORLD, ierr)

    call MPI_Send(description_a, DESCRIPTIONSIZE, MPI_CHARACTER, master_rank, 2, &
             MPI_COMM_WORLD, ierr)

    call MPI_Send(units_a, UNITSSIZE, MPI_CHARACTER, master_rank, 3, &
             MPI_COMM_WORLD, ierr)

    call MPI_Send(group_a, NAMESIZE, MPI_CHARACTER, master_rank, 4, &
             MPI_COMM_WORLD, ierr)

    data_shape(1) = 0
    data_shape(2) = len(data)
    call MPI_Send(data_shape, 2, MPI_INTEGER, master_rank, 5, &
             MPI_COMM_WORLD, ierr)

    call MPI_Send(data, data_shape(2), MPI_CHARACTER, &
             master_rank, 6, MPI_COMM_WORLD, ierr)

  end subroutine send_dataset_code

subroutine send_dataset_scalar(dataset_name, gal_id, data, &
                               units, description, group)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    double precision, dimension(:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    character(len=*), optional, intent(in) :: group
    character(len=DESCRIPTIONSIZE) :: description_a
    character(len=NAMESIZE) :: name_a, group_a
    character(len=UNITSSIZE) :: units_a
    integer, dimension(2) :: data_shape
    integer :: ierr

    group_a = 'Output'; description_a = ''; units_a = ''
    if (present(group))  group_a = group
    if (present(description))  description_a = description
    if (present(units))  units_a = units
    name_a = dataset_name
    ! If this is the master, simply write things to disk
    if (mpirank==0) then
      call IO_write_dataset(dataset_name, gal_id, data, units_a, description_a, group_a)
      return
    endif

    call MPI_Send(name_a, NAMESIZE, MPI_CHARACTER, master_rank, 1, &
                  MPI_COMM_WORLD, ierr)
    call MPI_Send(description_a, DESCRIPTIONSIZE, MPI_CHARACTER, master_rank, 2, &
                  MPI_COMM_WORLD, ierr)

    call MPI_Send(units_a, UNITSSIZE, MPI_CHARACTER, master_rank, 3, &
             MPI_COMM_WORLD, ierr)

    call MPI_Send(group_a, NAMESIZE, MPI_CHARACTER, master_rank, 4, &
             MPI_COMM_WORLD, ierr)

    data_shape(1) = 1
    data_shape(2) = size(data)
    call MPI_Send(data_shape, 2, MPI_INTEGER, master_rank, 5, &
             MPI_COMM_WORLD, ierr)
    call MPI_Send(data, data_shape(2), MPI_DOUBLE_PRECISION, &
             master_rank, 6, MPI_COMM_WORLD, ierr)

  end subroutine send_dataset_scalar

end module
