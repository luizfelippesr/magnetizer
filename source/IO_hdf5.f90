!# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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
! Contains a module that implements hdf5 IO
module IO
  use hdf5
  use messages
  implicit none

  interface IO_write_dataset ! Overload the IO_write_dataset function
    module procedure IO_write_dataset_scalar
    module procedure IO_write_dataset_vector
    module procedure IO_write_dataset_code
  end interface

  private

  ! Number of galaxies
  integer :: gals_number, grid_points

  ! Is this a new run or are we continuing a previous one?
  logical :: lresuming_run


  ! Parameter setting the maximum possible number of datasets
  integer, parameter :: max_number_of_datasets=50
  ! Maximum length of a text attribute
  integer, parameter :: ATTRIB_MAX_LEN = 200

  integer(hid_t) :: file_id          ! hdf5 file identifier
  integer(hid_t) :: file_id_out      ! hdf5 file identifier
  integer(hid_t) :: output_group_id  ! hdf5 group identifier
  integer(hid_t) :: input_group_id   ! hdf5 group identifier
  integer(hid_t) :: log_group_id

  ! The following store dataset names and ids
  character(len=15), dimension(max_number_of_datasets) :: dset_names = ''
  integer(hid_t), dimension(max_number_of_datasets) :: dset_ids
  integer(hid_t), dimension(max_number_of_datasets) :: dataspace_ids
  integer(hid_t), dimension(max_number_of_datasets) :: memspace_ids

  integer :: ndsets = 0 ! number of datasets
  logical :: Initialized = .false. ! Initialisation flag
  logical :: lchunking, lcompression, lseparate_output
  integer :: chunksize, compression_level
  public IO_start, IO_start_galaxy, IO_finish_galaxy, IO_write_dataset, IO_end
  public IO_read_dataset_scalar, IO_read_dataset_single, IO_flush

contains

  subroutine IO_start(mpi_comm, mpi_info, resuming_run, time_str)
    use grid
    use global_input_parameters
    ! Initializes IO
    integer, intent(in) :: mpi_comm, mpi_info
    logical, intent(in) :: resuming_run
    integer(hid_t) :: plist_id      ! property list identifier
    integer :: error
    integer, parameter :: MAX_NMLSTR=20000
    character(len=MAX_NMLSTR) :: namelist_str
    character(len=*), intent(in) :: time_str

    ! Reads properties from the global parameters module
    grid_points = nx
    lchunking = p_IO_chunking
    lcompression = p_IO_compression
    compression_level = p_IO_compression_level
    chunksize = p_IO_number_of_galaxies_in_chunks
    lseparate_output = p_IO_separate_output

    lresuming_run = resuming_run

    if (Initialized) then
      call error_message('IO_start', &
                        'Fatal Error: start_IO trying to initialize already '&
                        //' initialized IO', abort=.true.)
    endif

    ! Initializes predefined datatypes
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, error)
    ! Opens the file collectively.
    call message('Opening file'//trim(input_file_name), info=3)
    call h5fopen_f(trim(input_file_name), H5F_ACC_RDWR_F, &
                    file_id, error, access_prp=plist_id)
    call check(error)

    ! Reads some properties from the global attributes
    gals_number = read_int_attribute(file_id,'Number of galaxies')
    ngals = gals_number
    number_of_redshifts = read_int_attribute(file_id,'Number of snapshots')

    if (lseparate_output) then
      ! Opens if file exists, creates, otherwise
      if (lresuming_run) then
        call message('Opening file'//trim(output_file_name), info=3)
        call h5fopen_f(trim(output_file_name), H5F_ACC_RDWR_F, &
                        file_id_out, error, access_prp=plist_id)
        call check(error)
      else
        call message('Creating file'//trim(output_file_name), info=3)
        call h5fcreate_f(trim(output_file_name), H5F_ACC_TRUNC_F, &
                        file_id_out, error, access_prp=plist_id)
        call check(error)
      endif
    else
      lresuming_run = .false. ! NB at the moment, resuming a run with
                              ! unified output/input files is NOT
                              ! supported! This will be added later.
      file_id_out = file_id
    endif

    ! Opens input group
    call h5gopen_f(file_id, 'Input', input_group_id, error)
    call check(error)

    if (.not.lresuming_run) then
      ! Saves global parameters file to the HDF5 file
      call add_text_attribute(file_id_out, 'Model name', model_name)

      write(namelist_str, nml=run_parameters)
      call add_text_attribute(file_id_out,'run_parameters',trim(namelist_str))

      write(namelist_str, nml=io_parameters)
      call add_text_attribute(file_id_out,'io_parameters',trim(namelist_str))

      write(namelist_str, nml=grid_parameters)
      call add_text_attribute(file_id_out,'grid_parameters',trim(namelist_str))

      write(namelist_str, nml=dynamo_parameters)
      call add_text_attribute(file_id_out,'dynamo_parameters',trim(namelist_str))

      write(namelist_str, nml=ISM_and_disk_parameters)
      call add_text_attribute(file_id_out,'ISM_and_disk_parameters',trim(namelist_str))

      write(namelist_str, nml=outflow_parameters)
      call add_text_attribute(file_id_out,'outflow_parameters',trim(namelist_str))

      call add_text_attribute(file_id_out,'File creation',trim(time_str))

      ! Creates log group
      call message('Creating groups', info=3)
      call h5gcreate_f(file_id_out, 'Log', log_group_id, error)
      call check(error)
      ! Creates output group
      call h5gcreate_f(file_id_out, 'Output', output_group_id, error)
      call check(error)
    else
      ! Opens output group
      call h5gopen_f(file_id_out, 'Output', output_group_id, error)
      call check(error)
      ! Opens log group
      call h5gopen_f(file_id_out, 'Log', log_group_id, error)
      call check(error)
    endif
    ! Closes property list
    call h5pclose_f(plist_id, error)
    call check(error)

    Initialized = .true.
    return

  end subroutine IO_start


  logical function IO_start_galaxy(gal_id)
    ! Returns true if the galaxy does NOT have a completion flag
    ! Returns false, otherwise.
    integer, intent(in) :: gal_id
    double precision, dimension(1) :: flag

    if (.not.Initialized) then
      call error_message('IO_start_galaxy','Fatal Error: IO not initialized',&
                         abort=.true.)
    endif

    IO_start_galaxy = .true.

    if (lresuming_run) then
      ! If continuing a previous Magnetizer run, checks whether the galaxy
      ! was finished earlier.
      call IO_read_dataset_scalar('completed', gal_id, flag, is_log=.true.)
      if (flag(1)>0.5d0) then
        IO_start_galaxy = .false.
        call message('galaxy previously run. Skipping.', gal_id=gal_id,info=2)
      endif
    endif

  end function IO_start_galaxy

  subroutine IO_write_dataset_scalar(dataset_name, gal_id, data, &
                                     units, description, is_log)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    double precision, dimension(:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    logical, intent(in), optional :: is_log
    integer(hid_t) :: group_id
    integer ::  idx, error
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset
    integer, dimension(1) :: data_shape

    integer(hsize_t), dimension(2) :: dimsf_sca
    integer(hsize_t), dimension(2) :: dimsf_sca_1gal

    group_id = output_group_id
    data_shape = shape(data)

    ! Sets dataset dimensions.
    dimsf_sca = [data_shape(1),gals_number]
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_sca_1gal = [data_shape(1),1]

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      ! Flag 'is_log' can be used to choose the group 'Log' instead of 'Output'
      if (present(is_log)) then
          if (is_log) group_id = log_group_id
      endif
      if (lresuming_run) then
        idx = open_dset(dataset_name, group_id=group_id)
      else
        idx = create_dset(dataset_name, scalar=.true., dimsf_sca=dimsf_sca, &
                          group_id=group_id)
        ! Also writes the attributes (if needed)
        if (present(units)) &
          call add_text_attribute(dset_ids(idx), 'Units', units)
        if (present(description)) &
          call add_text_attribute(dset_ids(idx), 'Description', description)
      endif

      call H5Screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      call H5Dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if

    ! Selects hyperslab in the file.
    offset = [0,gal_id-1] ! Means: in second dimension, start from gal_id-1
    call H5Sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_sca_1gal, error)
    ! At last, writes the dataset
    call H5Dwrite_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error,&
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))

    if (error/=0) then
      call error_message('IO_write_dataset_scalar','Error. Exiting..', &
                        gal_id=gal_id, info=0, abort=.true.)
    endif

    return

  end subroutine IO_write_dataset_scalar

subroutine IO_read_dataset_scalar(dataset_name, gal_id, data, nrows, is_log, &
                                  group_id)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    integer, intent(in), optional :: nrows
    logical, intent(in), optional :: is_log
    integer(hid_t), optional :: group_id

    logical :: full
    double precision, dimension(:), intent(out) :: data
    integer ::  idx, error, nrows_actual
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset
    integer, dimension(1) :: data_shape
    integer(hsize_t), dimension(2) :: dimsf_sca
    integer(hsize_t), dimension(2) :: dimsf_sca_1gal
    logical :: is_log_actual

    if (present(is_log)) then
      is_log_actual = is_log
    else
      is_log_actual = .false.
    endif

    if (present(nrows)) then
        nrows_actual = nrows
        full = .true.
    else
        nrows_actual = gals_number
        full = .false.
    endif
    data_shape = shape(data)
    ! Sets dataset dimensions.
    dimsf_sca = [data_shape(1),gals_number]
    ! Sets the dimensions associated with reading a single galaxy
    dimsf_sca_1gal = [data_shape(1),1]

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      if (is_log_actual) then
        idx = open_dset(dataset_name, log_group_id)
      elseif (present(group_id)
        idx = open_dset(dataset_name, group_id)
      else
        idx = open_dset(dataset_name) !Defaults to the input_group_id
      endif

      if (full) then
        dimsf_sca = [nrows_actual,1]
        call H5Screate_simple_f(rank, dimsf_sca, memspace_ids(idx), error)
      else
        call H5Screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      endif
      call H5Dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if

    ! Selects hyperslab in the file.
    offset = [0,gal_id-1] ! Means: in second dimension, start from gal_id-1
    if (.not.full) &
      call H5Sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                                 offset, dimsf_sca_1gal, error)
    ! At last, reads the dataset
    if (full) then
      call H5Dread_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error)
    else
      call H5Dread_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error,&
                      file_space_id=dataspace_ids(idx), &
                      mem_space_id=memspace_ids(idx))
    endif
    return

  end subroutine IO_read_dataset_scalar

  subroutine IO_read_dataset_single(dataset_name, gal_id, data)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    double precision, intent(out) :: data
    integer ::  idx, error
    integer, parameter :: rank = 1
    integer(hssize_t), dimension(1) :: offset

    integer(hsize_t), dimension(1) :: dimsf_sca
    integer(hsize_t), dimension(1) :: dimsf_sca_1gal

    ! Sets dataset dimensions.
    dimsf_sca = [gals_number]
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_sca_1gal = [1]

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)

    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      idx = open_dset(dataset_name)
      call H5Screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      call H5Dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if

    ! Selects hyperslab in the file.
    offset = [gal_id-1] ! Means: in second dimension, start from gal_id-1
    call H5Sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_sca_1gal, error)
    ! At last, reads the dataset
    call H5Dread_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error,&
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))
    return

  end subroutine IO_read_dataset_single


  subroutine IO_write_dataset_vector(dataset_name, gal_id, data, &
                                     units, description)
    ! Writes a dataset to disk - vector version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    double precision, dimension(:,:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    integer ::  idx, error
    integer, parameter :: rank = 3
    integer(hssize_t), dimension(3) :: offset
    integer, dimension(2) :: data_shape
    integer(hsize_t), dimension(3) :: dimsf_vec
    integer(hsize_t), dimension(3) :: dimsf_vec_1gal

    data_shape = shape(data)

    ! Sets dataset dimensions.
    dimsf_vec = [data_shape(1),data_shape(2),gals_number]
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_vec_1gal = [data_shape(1),data_shape(2),1]

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      if (lresuming_run) then
        idx = open_dset(dataset_name, group_id=output_group_id)
      else
        idx = create_dset(dataset_name, dimsf_vec=dimsf_vec)
        ! Also writes the attributes (if needed)
        if (present(units)) &
          call add_text_attribute(dset_ids(idx), 'Units', units)
        if (present(description)) &
          call add_text_attribute(dset_ids(idx), 'Description', description)
      endif

      call H5Screate_simple_f(rank, dimsf_vec_1gal, memspace_ids(idx), error)
      call H5Dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if
    ! Selects hyperslab in the file.
    offset = [0,0,gal_id-1] ! Means: in third dimension, start from gal_id-1
    call H5Sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_vec_1gal, error)
    ! At last, writes the dataset
    call H5Dwrite_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_vec, error, &
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))

    if (error/=0) then
      call error_message('IO_write_dataset_vector','Fatal error. Exiting..', &
                         gal_id=gal_id, info=0, abort=.true.)
    endif

    return

  end subroutine IO_write_dataset_vector

  subroutine IO_write_dataset_code(dataset_name, gal_id, data, &
                                     units, description, is_log)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id
    character, dimension(:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    logical, intent(in), optional :: is_log
    integer(hid_t) :: group_id
    integer ::  idx, error
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset
    integer, dimension(1) :: data_shape

    integer(hsize_t), dimension(2) :: dimsf_sca
    integer(hsize_t), dimension(2) :: dimsf_sca_1gal


    group_id = output_group_id
    if (present(is_log)) then
        if (is_log) group_id = log_group_id
    endif
    data_shape = shape(data)

    ! Sets dataset dimensions.
    dimsf_sca = [data_shape(1),gals_number]
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_sca_1gal = [data_shape(1),1]

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      if (lresuming_run) then
        idx = open_dset(dataset_name, group_id=group_id)
      else
        ! Flag 'is_log' can be used to choose the group 'Log' instead of 'Output'
        idx = create_dset(dataset_name, scalar=.true., dimsf_sca=dimsf_sca, &
                            group_id=group_id, datatype=H5T_NATIVE_CHARACTER)
        ! Also writes the attributes (if needed)
        if (present(units)) &
          call add_text_attribute(dset_ids(idx), 'Units', units)
        if (present(description)) &
          call add_text_attribute(dset_ids(idx), 'Description', description)
      endif
      call H5Screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      call H5Dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if

    ! Selects hyperslab in the file.
    offset = [0,gal_id-1] ! Means: in second dimension, start from gal_id-1
    call H5Sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_sca_1gal, error)
    ! At last, writes the dataset
    call H5Dwrite_f(dset_ids(idx), H5T_NATIVE_CHARACTER, data, dimsf_sca, error,&
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))
    if (error/=0) then
      call error_message('IO','Error in IO_write_dataset_code. Exiting..', &
                         gal_id=gal_id, info=0, abort=.true.)
    endif

    return

  end subroutine IO_write_dataset_code


  subroutine IO_finish_galaxy(gal_id)
    ! Finishes IO
    ! At the moment, this is a dummy
    integer, intent(in) :: gal_id

    call message('Finished writing datasets',gal_id=gal_id, info=2)

  end subroutine IO_finish_galaxy


  subroutine IO_end(date)
    ! Finishes IO (closing everything)
    integer :: i, error
    character(len=*), intent(in) :: date
    character(len=10) :: date_formatted

!     date_formatted = date(1:4)//'-'//date(5:6)//'-'//date(7:8)
!     call add_text_attribute(file_id_out, 'Run date', date_formatted)

    ! Closes all dataspaces, namespaces and datasets
    ! (this may be moved to IO_finish_galaxy, if necessary)
    do i=1,ndsets
      call message('Preparing to close dataset '//dset_names(i), info=3)
      call message('Closing dataspace ', val_int=i, info=4)
      call H5Sclose_f(dataspace_ids(i), error)
      call check(error)
      call message('Closing memspace ', val_int=i, info=4)
      call H5Sclose_f(memspace_ids(i), error)
      call check(error)
      call message('Closing dataset '//dset_names(i), info=4)
      call H5Dclose_f(dset_ids(i), error)
      call check(error)
    end do
    ndsets = 0
    ! Closes output group
    call message("Closing Output group", info=3)
    call h5gclose_f(output_group_id, error)
    call check(error)
    ! Closes input group
    call message("Closing Input group", info=3)
    call h5gclose_f(input_group_id, error)
    call check(error)
    ! Closes log group
    call message("Closing log group", info=3)
    call h5gclose_f(log_group_id, error)
    call check(error)
    ! Closes file
    call message("Closing input file", info=3)
    call h5fclose_f(file_id, error)
    call check(error)
    ! Closes file
    if (lseparate_output) then
      call message("Closing Output file", info=3)
      call h5fclose_f(file_id_out, error)
      call check(error)
    endif
    ! Resets initialisation flag
    Initialized =.false.

  end subroutine IO_end

  subroutine IO_flush()
    integer :: error
    call message('Flushing HDF5 file(s).', info=3)

    call H5Fflush_f(file_id, H5F_SCOPE_GLOBAL_F , error)
    call check(error)
    if (lseparate_output) then
      call H5Fflush_f(file_id_out, H5F_SCOPE_GLOBAL_F , error)
      call check(error)
    endif
  end subroutine IO_flush

  subroutine add_int_attribute(dset_id, attribute_name, attribute)
    ! Adds or modifies int attributes to an existing dataset
    !
    ! Input: dset_id -> id of an opened dataset
    !        attribute_name -> string with the name
    !        attribute -> strig with the value
    character(len=*), intent(in) ::  attribute_name ! attribute name
    character(len=*), intent(in) ::  attribute ! attribute data
    integer(hid_t) :: dset_id       ! Dataset identifier
    integer(hid_t) :: attr_id       ! attribute identifier
    integer(hid_t) :: aspace_id     ! attribute dataspace identifier
    integer(hid_t) :: atype_id      ! attribute dataspace identifier
    integer(hsize_t), dimension(1) :: adims = [1] ! attribute dimension
    integer     ::   arank = 1                      ! attribute rank
    integer(size_t) :: attrlen    ! length of the attribute string
    integer     ::   error ! error flag
    integer(hsize_t), dimension(1) :: data_dims
    integer :: old_attr_value
    attrlen = len(attribute)
    data_dims(1) = 1
    ! Creates scalar data space for the attribute.
    call H5Screate_simple_f(arank, adims, aspace_id, error)
    call check(error)

    ! Creates datatype for the attribute.
    call H5Tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
    call check(error)
    call H5Tset_size_f(atype_id, attrlen, error)
    call check(error)

    ! Creates dataset attribute.
    call H5Acreate_f(dset_id, attribute_name, atype_id, aspace_id, attr_id, &
                     error)
    if (error/=0) then
      ! If the attribute already exists (which will be an error),
      ! opens it to overwritting
      call H5Aopen_f(dset_id, attribute_name, attr_id, error)
      call check(error)
    endif
    ! Writes the attribute data.
    call H5Awrite_f(attr_id, atype_id, attribute, data_dims, error)
    call check(error)
    ! Closes the attribute.
    call H5Aclose_f(attr_id, error)
    call check(error)
    ! Terminates access to the data space.
    call H5Sclose_f(aspace_id, error)
    call check(error)
  end subroutine add_int_attribute

  subroutine add_text_attribute(dset_id, attribute_name, attribute)
    ! Adds new text attributes to an existing dataset
    ! If the attribute already exists, checks whether it has the same value
    ! and fails if the values differ.
    !
    ! Input: dset_id -> id of an opened dataset
    !        attribute_name -> string with the name
    !        attribute -> strig with the value
    character(len=*), intent(in) ::  attribute_name ! attribute name
    character(len=*), intent(in) ::  attribute ! attribute data
    integer(hid_t) :: dset_id       ! Dataset identifier
    integer(hid_t) :: attr_id       ! attribute identifier
    integer(hid_t) :: aspace_id     ! attribute dataspace identifier
    integer(hid_t) :: atype_id      ! attribute dataspace identifier
    integer(hsize_t), dimension(1) :: adims = [1] ! attribute dimension
    integer     ::   arank = 1                      ! attribute rank
    integer(size_t) :: attrlen    ! length of the attribute string
    integer     ::   error ! error flag
    integer(hsize_t), dimension(1) :: data_dims
    character(len=ATTRIB_MAX_LEN) :: old_attr_value
    attrlen = len(attribute)
    data_dims(1) = 1
    ! Creates scalar data space for the attribute.
    call H5Screate_simple_f(arank, adims, aspace_id, error)
    call check(error)

    ! Creates datatype for the attribute.
    call H5Tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
    call check(error)
    call H5Tset_size_f(atype_id, attrlen, error)
    call check(error)

    ! Creates dataset attribute.
    call H5Acreate_f(dset_id, attribute_name, atype_id, aspace_id, attr_id, &
                     error)
    if (error/=0) then
      ! If the attribute already exists (which will be an error), reads
      ! its value and compares with the one that one is trying to write
      old_attr_value = read_text_attribute(dset_id, attribute_name)
      if (trim(old_attr_value)/=trim(attribute)) then
        call error_message("IO","Attribute alread exists.", info=0, abort=.true.)
      endif
    else
      ! Writes the attribute data.
      call H5Awrite_f(attr_id, atype_id, attribute, data_dims, error)
      call check(error)
      ! Closes the attribute.
      call H5Aclose_f(attr_id, error)
      call check(error)
    endif

    ! Terminates access to the data space.
    call H5Sclose_f(aspace_id, error)
    call check(error)

  end subroutine add_text_attribute


  function read_int_attribute(dset_id, attribute_name)
    integer(hid_t) :: attr_id, atype_id       ! attribute identifier
    integer(hid_t) :: dset_id   ! Dataset identifier
    integer     ::   error ! error flag
    integer(hsize_t), dimension(1) :: dimsf_sca
    integer :: read_int_attribute
    character(len=*), intent(in) ::  attribute_name ! attribute name

    dimsf_sca = 0

    call H5Tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
    call check(error)

    call H5Aopen_f(dset_id, attribute_name, attr_id, error)
    call check(error)

    call H5Aread_f(attr_id, atype_id, read_int_attribute, dimsf_sca,error)
    call check(error)

    ! Closes the attribute.
    call H5Aclose_f(attr_id, error)
    call check(error)

  end function read_int_attribute

    function read_text_attribute(dset_id, attribute_name)
    integer(hid_t) :: attr_id, atype_id       ! attribute identifier
    integer(hid_t) :: dset_id ! Dataset identifier
    integer     ::   error ! error flag
    integer(hsize_t), dimension(1) :: dimsf_sca
    character(len=ATTRIB_MAX_LEN) :: read_text_attribute
    character(len=*), intent(in) ::  attribute_name ! attribute name

    dimsf_sca = 0
    call H5Tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
    call check(error)

    call H5Aopen_f(dset_id, attribute_name, attr_id, error)
    call check(error)

    call H5Aread_f(attr_id, atype_id, read_text_attribute, dimsf_sca, error)
    call check(error)

    ! Closes the attribute.
    call H5Aclose_f(attr_id, error)
    call check(error)

  end function read_text_attribute



  function create_dset(dataset_name, scalar, dimsf_vec, dimsf_sca, group_id, &
                       datatype) result(idx)
    ! Creates a hdf5 dataspace, a namespace and dataset
    ! Returns the index
    ! NB This is a collective procedure
    integer(hsize_t), dimension(3), intent(in), optional :: dimsf_vec
    integer(hsize_t), dimension(2), intent(in), optional :: dimsf_sca
    integer(hid_t), intent(in), optional :: group_id
    integer(hsize_t), dimension(3), target :: chunkdim_vec
    integer(hsize_t), dimension(2), target :: chunkdim_sca
    character(len=*), intent(in) :: dataset_name
    integer(hid_t), intent(in), optional :: datatype
    logical, intent(in), optional :: scalar
    integer(hid_t) :: plist_id, group_id_actual      ! property list identifier
    integer(hid_t) :: dataspace ! dataspace identifier (temporary)
    integer(hid_t) :: datatype_actual
    logical :: scalar_actual
    integer :: idx
    integer :: rank, error
    integer(hsize_t), parameter :: REDSHIFTS_PER_CHUNK = 1


    ! Increments ndsets and sets idx
    ndsets = ndsets+1
    idx = ndsets
    ! Adds name to the list
    dset_names(idx) = dataset_name

    if (present(scalar)) then
      scalar_actual = scalar
    else
      scalar_actual = .false.
    endif

    if (present(datatype)) then
      datatype_actual = datatype
    else
      datatype_actual = H5T_NATIVE_DOUBLE
    endif

    if (present(group_id)) then
      group_id_actual = group_id
    else
      group_id_actual = output_group_id
    endif

    call message('Creating dataset '//dataset_name//', idx = ',val_int= idx, &
                 info=4)
    if (scalar_actual) then
      ! "scalar" here, means no variation with radius
      ! thus, 2 dimensions: galaxy id and time
      rank = 2
      if (lchunking) &
          chunkdim_sca = [ REDSHIFTS_PER_CHUNK, min(dimsf_sca(2), chunksize) ]
      call H5Screate_simple_f(rank, dimsf_sca, dataspace, error)
      call check(error)
    else
      ! 3 dimensions: galaxy id, time and radius
      rank = 3
      if (lchunking) chunkdim_vec = &
              [ REDSHIFTS_PER_CHUNK, dimsf_vec(2), min(dimsf_vec(3), chunksize) ]
      call H5Screate_simple_f(rank, dimsf_vec, dataspace, error)
      call check(error)
    endif

    ! Creates a dataset property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call check(error)

    if (lchunking) then
      ! Sets the size of the chunks
      if (scalar_actual) then
        call h5pset_chunk_f(plist_id, rank, chunkdim_sca, error)
        call check(error)
      else
        call h5pset_chunk_f(plist_id, rank, chunkdim_vec, error)
        call check(error)
      endif
      if (lcompression) then
        ! Sets compression
        call h5pset_deflate_f(plist_id, compression_level, error)
        call check(error)
      endif
    endif
    ! Creates the dataset
    call H5Dcreate_f(group_id_actual, dataset_name, datatype_actual, &
                    dataspace, dset_ids(idx), error,  dcpl_id=plist_id)
    call check(error)
    ! Closes dataspace
    call H5Sclose_f(dataspace, error)
    call check(error)

    return
  end function create_dset


  function open_dset(dataset_name, group_id) result(idx)
    ! Creates a hdf5 dataspace, a namespace and dataset
    ! Returns the index
    ! NB This is a collective procedure
    character(len=*), intent(in) :: dataset_name
    integer(hid_t), optional :: group_id
    integer :: idx, error
    integer(hid_t) :: group_id_actual

    if (present(group_id)) then
      group_id_actual = group_id
    else
      group_id_actual = input_group_id
    endif

    ! Increments ndsets and sets idx
    ndsets = ndsets+1
    idx = ndsets
    ! Adds name to the list
    dset_names(idx) = dataset_name

    ! Opens the dataset
    call H5Dopen_f(group_id_actual, dataset_name, dset_ids(idx), error)
    call check(error)

    return
  end function open_dset


  function find_dset(dset_name) result(idx)
    ! Returns index for a dataset given its name
    character(len=*), intent(in) :: dset_name
    integer :: idx
    integer :: i

    do i=1,ndsets
      if (dset_name==trim(dset_names(i))) then
        idx = i
        return
      endif
    end do
    idx = -1
    return
  end function find_dset

  subroutine check(error)
    integer, intent(in) :: error
    if (error /= 0) then
      call error_message("IO","Error when calling HDF5.", info=0, abort=.true.)
    endif
  end subroutine

end module IO
