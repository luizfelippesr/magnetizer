! Contains a module that implements hdf5 IO
module IO
  use hdf5
  implicit none
  
  interface IO_write_dataset ! Overload the IO_write_dataset function
    module procedure IO_write_dataset_scalar
    module procedure IO_write_dataset_vector
  end interface
  
  private
  
  ! Parameter setting the maximum possible number of datasets
  integer, parameter :: max_number_of_datasets=25
  
  integer(hid_t) :: file_id       ! hdf5 file identifier 
  
  ! The following store dataset names and ids
  character(len=15), dimension(max_number_of_datasets) :: dset_names = ''
  integer(hid_t), dimension(max_number_of_datasets) :: dset_ids 
  integer(hid_t), dimension(max_number_of_datasets) :: dataspace_ids   
  integer(hid_t), dimension(max_number_of_datasets) :: memspace_ids     
  
  integer(hsize_t), dimension(3) :: dimsf_vec
  integer(hsize_t), dimension(2) :: dimsf_sca
  integer(hsize_t), dimension(3) :: dimsf_vec_1gal
  integer(hsize_t), dimension(2) :: dimsf_sca_1gal
  
  integer :: ndsets = 0 ! number of datasets
  logical :: Initialized = .false. ! Initialisation flag
  
  public IO_start, IO_start_galaxy, IO_finish_galaxy, IO_write_dataset, IO_end
        
contains

  subroutine IO_start(path_to_model, output_file, mpi_comm, mpi_info)
    use grid
    use global_input_parameters
    ! Initializes IO
    ! NB at the moment, MPI is required (will be relaxed later)
    character(len=*), intent(in) :: path_to_model, output_file
    character(len=40) :: filename
    integer, intent(in), optional :: mpi_comm, mpi_info
    integer(hid_t) :: plist_id      ! property list identifier 
    integer :: error
    
    if (.not.present(mpi_comm)) then
      error stop 'Fatal Error: start_IO, trying to initialize parallel hdf5 IO without a communicator.'
    endif
    
    if (Initialized) then
      error stop 'Fatal Error: start_IO trying to initialize already initialized IO'
    endif
    
    filename = path_to_model // '/' // output_file
    
    ! Initializes predefined datatypes
    call h5open_f(error) 
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, mpi_comm, mpi_info, error)
    ! Creates the file collectively.
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error, &
                     access_prp=plist_id)
    ! Creates property list
    call h5pclose_f(plist_id, error)
    
    ! Sets dataset dimensions.
    dimsf_vec = (/nx,max_number_of_redshifts,ngals/) 
    dimsf_sca = (/max_number_of_redshifts,ngals/) 
    ! Sets the dimensions associated with writing a single galaxy
    dimsf_vec_1gal = (/nx,max_number_of_redshifts,1/) 
    dimsf_sca_1gal = (/max_number_of_redshifts,1/) 
    
    Initialized = .true.
    
    return    
    
  end subroutine IO_start
  
  subroutine IO_start_galaxy(gal_id, info)
    ! Initializes galaxy 
    ! At the moment, this is a dummy
    integer, intent(in) :: gal_id, info

    if (.not.Initialized) then
      error stop  'Fatal Error: IO_write_meta, IO not initialized'
    endif

  end subroutine IO_start_galaxy

  
  subroutine IO_write_dataset_scalar(dataset_name, gal_id, info, data)
    ! Writes a dataset to disk - scalar version
    
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    double precision, dimension(:), intent(in) :: data
    integer ::  idx, error
    integer, parameter :: rank = 2
    integer(hssize_t), dimension(2) :: offset 

    ! Tries to find a previously opened dataset (-1 signals new)
    idx = find_dset(dataset_name)
    ! If it wasn't previously opened, creates it (collectively)
    if (idx < 0) then
      idx = create_dset(dataset_name, scalar=.true.)
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
    
    
  subroutine IO_write_dataset_vector(dataset_name, gal_id, info, data)
    ! Writes a dataset to disk

    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    double precision, dimension(:,:), intent(in) :: data
    integer ::  idx, error
    integer, parameter :: rank = 3
    integer(hssize_t), dimension(3) :: offset 
! 
!     ! Tries to find a previously opened dataset (-1 signals new)
!     idx = find_dset(dataset_name)
!     ! If it wasn't previously opened, creates it (collectively)
!     if (idx < 0) then
!       idx = create_dset(dataset_name)
!       call h5screate_simple_f(rank, dimsf_vec_1gal, memspace_ids(idx), error) 
!       call h5dget_space_f(dset_ids(idx), dataspace_ids(idx), error)
!     end if
!     
!     ! Selects hyperslab in the file.
!     offset = (/0,0,gal_id-1/) ! Means: in third dimension, start from gal_id-1
!     call h5sselect_hyperslab_f(dataspace_ids(idx), H5S_SELECT_SET_F, & 
!                                offset, dimsf_vec_1gal, error)
!     ! At last, writes the dataset
!     call h5dwrite_f(dset_ids(idx), H5T_NATIVE_integer, data, dimsf_vec, error, &
!                     file_space_id=dataspace_ids(idx), &
!                     mem_space_id=memspace_ids(idx))
                    
    return    
    
  end subroutine IO_write_dataset_vector
    
  
  subroutine IO_finish_galaxy(gal_id, info)
    ! Finishes IO 
    ! At the moment, this is a dummy
    integer, intent(in) :: gal_id, info
    
  end subroutine IO_finish_galaxy
  
  
  subroutine IO_end()  
    ! Finishes IO (closing everything)
    integer :: i, error

    ! Closes all dataspaces, namespaces and datasets
    ! (this may be moved to IO_finish_galaxy, if necessary)
    do i=1,ndsets
      print *, 'Closing dataspace ', i
      call h5sclose_f(dataspace_ids(i), error)
      
      print *, 'Closing memspace ', i
      call h5sclose_f(memspace_ids(i), error)
      
      print *, 'Closing dset ', i
      call h5dclose_f(dset_ids(i), error)
    end do
    ndsets = 0
    ! Closes file
    call h5fclose_f(file_id, error)    

  end subroutine IO_end  
  
  
  function create_dset(dataset_name, scalar) result(idx)
    ! Creates a hdf5 dataspace, a namespace and dataset
    ! Returns the index
    ! NB This is a collective procedure
    character(len=*), intent(in) :: dataset_name
    logical, intent(in), optional :: scalar
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
      call h5screate_simple_f(rank, dimsf_sca, dataspace, error)
    else
      ! 3 dimensions: galaxy id, time and radius
      rank = 3
      call h5screate_simple_f(rank, dimsf_vec, dataspace, error)
    endif
    
    ! Creates the dataset with default properties.
    call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace, &
                     dset_ids(idx), error)
    ! Closes dataspace
    call h5sclose_f(dataspace, error)
  
    return  
  end function create_dset
  
  
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