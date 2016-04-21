! Contains a module that implements hdf5 IO
module IO
  use hdf5
  implicit none
  
  interface IO_write_dataset ! Overload the IO_write_dataset function
    module procedure IO_write_dataset_scalar
    module procedure IO_write_dataset_vector
  end interface
  
  private
  
  ! Number of galaxies
  integer :: gals_number, grid_points

  ! Parameter setting the maximum possible number of datasets
  integer, parameter :: max_number_of_datasets=55
  
  integer(hid_t) :: file_id          ! hdf5 file identifier
  integer(hid_t) :: file_id_out      ! hdf5 file identifier
  integer(hid_t) :: output_group_id  ! hdf5 group identifier
  integer(hid_t) :: input_group_id   ! hdf5 group identifier
  
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
  public IO_read_dataset_scalar, IO_read_dataset_single

contains

  subroutine IO_start(mpi_comm, mpi_info)
    use grid
    use global_input_parameters
    ! Initializes IO
    character(len=40) :: filepath
    integer, intent(in), optional :: mpi_comm, mpi_info
    integer(hid_t) :: plist_id      ! property list identifier
    integer :: error
    
    ! Reads properties from the global parameters module
    grid_points = nx
    lchunking = p_IO_chunking
    lcompression = p_IO_compression
    compression_level = p_IO_compression_level
    chunksize = p_IO_number_of_galaxies_in_chunks
    lseparate_output = p_IO_separate_output

    if (.not.present(mpi_comm)) then
      stop 'Fatal Error: start_IO, trying to initialize parallel hdf5 IO without a communicator.'
    endif

    if (Initialized) then
      stop 'Fatal Error: start_IO trying to initialize already initialized IO'
    endif

    ! Initializes predefined datatypes
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, error)
    ! Opens the file collectively.
    call h5fopen_f(trim(input_file_name), H5F_ACC_RDWR_F, &
                    file_id, error, access_prp=plist_id)

    ! Reads some properties from the global attributes
    gals_number = read_global_attribute('Number of galaxies')
    ngals = gals_number
    number_of_redshifts = read_global_attribute('Number of snapshots')

    if (lseparate_output) then
      ! Creates the file collectively.
      call h5fcreate_f(trim(output_file_name), H5F_ACC_TRUNC_F, &
                        file_id_out, error, access_prp=plist_id)
    else
      file_id_out = file_id
    endif

    ! Creates output group
    call h5gcreate_f(file_id_out, 'Output', output_group_id, error)
    ! Opens input group
    call h5gopen_f(file_id, 'Input', input_group_id, error)

    ! Closes property list
    call h5pclose_f(plist_id, error)

    Initialized = .true.
    return    
    
  end subroutine IO_start
  

  subroutine IO_start_galaxy(gal_id, info)
    ! Initializes galaxy
    ! At the moment, this is a dummy
    integer, intent(in) :: gal_id, info

    if (.not.Initialized) then
      stop  'Fatal Error: IO_write_meta, IO not initialized'
    endif

  end subroutine IO_start_galaxy


  subroutine IO_write_dataset_scalar(dataset_name, gal_id, info, data, &
                                     units, description)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    double precision, dimension(:), intent(in) :: data
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: description
    integer ::  idx, error
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset
    integer, dimension(1) :: data_shape

    integer(hsize_t), dimension(2) :: dimsf_sca
    integer(hsize_t), dimension(2) :: dimsf_sca_1gal

    data_shape = shape(data)

    ! Sets dataset dimensions.
    dimsf_sca = (/data_shape(1),gals_number/)
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_sca_1gal = (/data_shape(1),1/)

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      idx = create_dset(dataset_name, scalar=.true., dimsf_sca=dimsf_sca)
      ! Also writes the attributes (if needed)
      if (present(units)) &
        call add_text_attribute(dset_ids(idx), 'Units', units)
      if (present(description)) &
        call add_text_attribute(dset_ids(idx), 'Description', description)

      call h5screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      call h5dget_space_f(dset_ids(idx), dataspace_ids(idx), error)

    end if

    ! Selects hyperslab in the file.
    offset = (/0,gal_id-1/) ! Means: in second dimension, start from gal_id-1
    call h5sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_sca_1gal, error)
    ! At last, writes the dataset
    call h5dwrite_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error,&
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))

    return

  end subroutine IO_write_dataset_scalar

subroutine IO_read_dataset_scalar(dataset_name, gal_id, info, data, nrows)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    integer, intent(in), optional :: nrows
    logical :: full
    double precision, dimension(:), intent(out) :: data
    integer ::  idx, error, nrows_actual
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset
    integer, dimension(1) :: data_shape
    integer(hsize_t), dimension(2) :: dimsf_sca
    integer(hsize_t), dimension(2) :: dimsf_sca_1gal
    
    if (present(nrows)) then
        nrows_actual = nrows
        full = .true.
    else
        nrows_actual = gals_number
        full = .false.
    endif
    data_shape = shape(data)
    ! Sets dataset dimensions.
    dimsf_sca = (/data_shape(1),gals_number/)
    ! Sets the dimensions associated with reading a single galaxy
    dimsf_sca_1gal = (/data_shape(1),1/)

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      idx = open_dset(dataset_name)
      if (full) then
        dimsf_sca = (/nrows_actual,1/)
        call h5screate_simple_f(rank, dimsf_sca, memspace_ids(idx), error)
      else
        call h5screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      endif
      call h5dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if

    ! Selects hyperslab in the file.
    offset = (/0,gal_id-1/) ! Means: in second dimension, start from gal_id-1
    if (.not.full) &
      call h5sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                                 offset, dimsf_sca_1gal, error)
    ! At last, reads the dataset
    if (full) then
      call h5dread_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error)
    else
      call h5dread_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error,&
                      file_space_id=dataspace_ids(idx), &
                      mem_space_id=memspace_ids(idx))
    endif
    return

  end subroutine IO_read_dataset_scalar

  subroutine IO_read_dataset_single(dataset_name, gal_id, info, data)
    ! Writes a dataset to disk - scalar version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    double precision, intent(out) :: data
    integer ::  idx, error
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset

    integer(hsize_t), dimension(2) :: dimsf_sca
    integer(hsize_t), dimension(2) :: dimsf_sca_1gal

    ! Sets dataset dimensions.
    dimsf_sca = (/1,gals_number/)
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_sca_1gal = (/1,1/)

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)

    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      idx = open_dset(dataset_name)
      call h5screate_simple_f(rank, dimsf_sca_1gal, memspace_ids(idx), error)
      call h5dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if

    ! Selects hyperslab in the file.
    offset = (/0,gal_id-1/) ! Means: in second dimension, start from gal_id-1
    call h5sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_sca_1gal, error)
    ! At last, reads the dataset
    call h5dread_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_sca, error,&
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))
    return

  end subroutine IO_read_dataset_single


  subroutine IO_write_dataset_vector(dataset_name, gal_id, info, data, &
                                     units, description)
    ! Writes a dataset to disk - vector version
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
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
    dimsf_vec = (/data_shape(1),data_shape(2),gals_number/)
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_vec_1gal = (/data_shape(1),data_shape(2),1/)

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      idx = create_dset(dataset_name, dimsf_vec=dimsf_vec)
      ! Also writes the attributes (if needed)
      if (present(units)) &
        call add_text_attribute(dset_ids(idx), 'Units', units)
      if (present(description)) &
        call add_text_attribute(dset_ids(idx), 'Description', description)

      call h5screate_simple_f(rank, dimsf_vec_1gal, memspace_ids(idx), error) 
      call h5dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
    end if
    
    ! Selects hyperslab in the file.
    offset = (/0,0,gal_id-1/) ! Means: in third dimension, start from gal_id-1
    call h5sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, &
                               offset, dimsf_vec_1gal, error)
    ! At last, writes the dataset
    call h5dwrite_f(dset_ids(idx), H5T_NATIVE_DOUBLE, data, dimsf_vec, error, &
                    file_space_id=dataspace_ids(idx), &
                    mem_space_id=memspace_ids(idx))
    return    
    
  end subroutine IO_write_dataset_vector
    
  
  subroutine IO_finish_galaxy(gal_id, info)
    ! Finishes IO 
    ! At the moment, this is a dummy
    integer, intent(in) :: gal_id, info
    
  end subroutine IO_finish_galaxy
  
  
  subroutine IO_end(info)
    ! Finishes IO (closing everything)
    integer, intent(in) :: info
    integer :: i, error

    ! Closes all dataspaces, namespaces and datasets
    ! (this may be moved to IO_finish_galaxy, if necessary)
    do i=1,ndsets
      if (info>2) print *, 'Closing ', dset_names(i)
!       print *, 'Closing dataspace ', i
      call h5sclose_f(dataspace_ids(i), error)
      
!       print *, 'Closing memspace ', i
      call h5sclose_f(memspace_ids(i), error)
      
!       print *, 'Closing dset ', i
      call h5dclose_f(dset_ids(i), error)
    end do
    ndsets = 0
    ! Closes output group
    if (info>2) print *, "Closing Output group"
    call h5gclose_f(output_group_id, error)
    ! Closes input group
    if (info>2) print *, "Closing input group"
    call h5gclose_f(input_group_id, error)
    ! Closes file
    if (info>2) print *, "Closing input file"
    call h5fclose_f(file_id, error)
    ! Closes file
    if (lseparate_output) then
        if (info>2) print *, "Closing Output file"
      call h5fclose_f(file_id_out, error)
    endif
  end subroutine IO_end  
  
  subroutine add_text_attribute(dset_id, attribute_name, attribute)
    ! Adds new text attributes to an existing dataset
    ! Input: dset_id -> id of an opened dataset
    !        attribute_name -> string with the name
    !        attribute -> strig with the value
    character(len=*), intent(in) ::  attribute_name ! attribute name
    character(len=*), intent(in) ::  attribute ! attribute data
    integer(hid_t) :: dset_id       ! Dataset identifier
    integer(hid_t) :: attr_id       ! attribute identifier
    integer(hid_t) :: aspace_id     ! attribute dataspace identifier
    integer(hid_t) :: atype_id      ! attribute dataspace identifier
    integer(hsize_t), dimension(1) :: adims = (/1/) ! attribute dimension
    integer     ::   arank = 1                      ! attribute rank
    integer(size_t) :: attrlen    ! length of the attribute string
    integer     ::   error ! error flag
    integer(hsize_t), dimension(1) :: data_dims


    attrlen = len(attribute)
    data_dims(1) = 1

    ! Creates scalar data space for the attribute.
    call h5screate_simple_f(arank, adims, aspace_id, error)

    ! Creates datatype for the attribute.
    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
    call h5tset_size_f(atype_id, attrlen, error)

    ! Creates dataset attribute.
    call h5acreate_f(dset_id, attribute_name, atype_id, aspace_id, attr_id, &
                     error)

    ! Writes the attribute data.
    call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)

    ! Closes the attribute.
    call h5aclose_f(attr_id, error)

    ! Terminates access to the data space.
    call h5sclose_f(aspace_id, error)

  end subroutine add_text_attribute


  function read_global_attribute(attribute_name)
    integer(hid_t) :: attr_id,atype_id       ! attribute identifier
    integer     ::   error ! error flag
    integer(hsize_t), dimension(1) :: dimsf_sca
    integer :: read_global_attribute
    character(len=*), intent(in) ::  attribute_name ! attribute name

    dimsf_sca = 0.0

    call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)

    call h5aopen_f(file_id, attribute_name, attr_id, error)

    call h5aread_f(attr_id, atype_id, read_global_attribute, dimsf_sca,error)

    ! Closes the attribute.
    call h5aclose_f(attr_id, error)

  end function read_global_attribute



  function create_dset(dataset_name, scalar, dimsf_vec, dimsf_sca) result(idx)
    ! Creates a hdf5 dataspace, a namespace and dataset
    ! Returns the index
    ! NB This is a collective procedure
    integer(hsize_t), dimension(3), intent(in), optional :: dimsf_vec
    integer(hsize_t), dimension(2), intent(in), optional :: dimsf_sca
    integer(hsize_t), dimension(3), target :: chunkdim_vec
    integer(hsize_t), dimension(2), target :: chunkdim_sca
    integer, pointer :: chunkdim
    character(len=*), intent(in) :: dataset_name
    logical, intent(in), optional :: scalar
    integer(hid_t) :: plist_id      ! property list identifier
    integer(hid_t) :: dataspace ! dataspace identifier (temporary)
    logical :: scalar_actual
    integer :: idx
    integer :: rank, error

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

    if (scalar_actual) then
      ! "scalar" here, means no variation with radius
      ! thus, 2 dimensions: galaxy id and time
      rank = 2
      if (lchunking) &
          chunkdim_sca = (/ dimsf_sca(1), min(dimsf_sca(2), chunksize) /)
      call h5screate_simple_f(rank, dimsf_sca, dataspace, error)
    else
      ! 3 dimensions: galaxy id, time and radius
      rank = 3

      if (lchunking) chunkdim_vec = &
              (/ dimsf_vec(1), dimsf_vec(2)/2, min(dimsf_vec(3), chunksize) /)
      call h5screate_simple_f(rank, dimsf_vec, dataspace, error)
    endif

    ! Creates a dataset property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    if (lchunking) then
      ! Sets the size of the chunks
      if (scalar_actual) then
        call h5pset_chunk_f(plist_id, rank, chunkdim_sca, error)
      else
        call h5pset_chunk_f(plist_id, rank, chunkdim_vec, error)
      endif
      if (lcompression) then
        ! Sets compression
        call h5pset_deflate_f(plist_id, compression_level, error)
      endif
    endif
    ! Creates the dataset
    call h5dcreate_f(output_group_id, dataset_name, H5T_NATIVE_DOUBLE, &
                    dataspace, dset_ids(idx), error,  dcpl_id=plist_id)
    ! Closes dataspace
    call h5sclose_f(dataspace, error)

    return
  end function create_dset


  function open_dset(dataset_name) result(idx)
    ! Creates a hdf5 dataspace, a namespace and dataset
    ! Returns the index
    ! NB This is a collective procedure
    character(len=*), intent(in) :: dataset_name
    integer(hid_t) :: dataspace ! dataspace identifier (temporary)
    integer :: idx
    integer :: rank, error

    ! Increments ndsets and sets idx
    ndsets = ndsets+1
    idx = ndsets
    ! Adds name to the list
    dset_names(idx) = dataset_name

    ! Opens the dataset
    call h5dopen_f(input_group_id, dataset_name, dset_ids(idx), error)

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
end module IO
