! Contains a module that implements simple 'write-to-text-file' IO
! May be used as a template for other IO-schemes
module IO
  implicit none
  
  interface IO_write_dataset ! Overload the IO_write_dataset function
    module procedure IO_write_dataset_scalar
    module procedure IO_write_dataset_vector
  end interface
  
  private
  logical :: Initialized = .false.
  
  integer, parameter :: output_unit = 12
  integer, parameter :: ts_unit = 13
  integer, parameter :: diag_unit = 20
  integer, parameter :: columns_unit = 21
  
  character(len=8) :: gal_id_string
  character(len=:), allocatable :: model_path   
    
  public IO_start, IO_start_galaxy,IO_finish_galaxy, IO_write_dataset, IO_end
        
contains

  subroutine IO_start(path_to_model, output_file, mpi_comm, mpi_info)
    ! Initializes IO
    ! The output file is NOT used with this IO settings
    
    character(len=*), intent(in) :: path_to_model, output_file   
    integer, intent(in), optional :: mpi_comm, mpi_info

    
    if (Initialized) then
      error stop 'Fatal Error: start_IO trying to initialize already initialized IO'
    endif
    
    Initialized = .true.
    
    ! Stores model directory in module variable
    allocate(character(len=len_trim(path_to_model)) :: model_path)    
    model_path = path_to_model
    
    return    
    
  end subroutine IO_start
  
  subroutine IO_start_galaxy(gal_id, info)
    ! Initializes galaxy / writes meta-information/diagnostics
    
    integer, intent(in) :: gal_id, info
    character(len=8) :: frmt
    logical, save :: first = .true.

    if (.not.Initialized) then
      error stop  'Fatal Error: IO_write_meta, IO not initialized'
    endif

    if (first) then
      open(unit=diag_unit,file= model_path // '/diagnostic.out', &
           status="replace")
    endif
    
    frmt='(I8.8)'
    write(gal_id_string,frmt) gal_id
    write(diag_unit,*)'Writing output for final timestep to file run',gal_id_string,'.out'
    
    if (info> 0) then
      print*,'Writing time series output to file ts',gal_id_string,'.out'
    endif
    
    ! Opens files
    open(unit=ts_unit,file= model_path // '/ts_' // gal_id_string // &
          '.out', status="replace")
    
    open(unit=columns_unit,file= model_path // '/columns_' // gal_id_string // &
          '.out', status="replace")      
        
  end subroutine IO_start_galaxy
    
    
  subroutine IO_write_dataset_scalar(dataset_name, gal_id, info, data)
    ! Writes a dataset to disk - scalar version
    ! This is a wrapper that calls the vector version
    
    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    double precision, dimension(:), intent(in) :: data
    double precision, dimension(size(data),1) :: data_vec
    
    data_vec(:,1) = data
    call IO_write_dataset_vector(dataset_name, gal_id, info, data_vec)
    
  end subroutine IO_write_dataset_scalar
    
    
  subroutine IO_write_dataset_vector(dataset_name, gal_id, info, data)
    ! Writes a dataset to disk

    character(len=*), intent(in) :: dataset_name
    integer, intent(in) :: gal_id, info
    double precision, dimension(:,:), intent(in) :: data
    
    ! Writes the name of the dataset to a file (to allow for later retrieval)
    write(columns_unit,*) dataset_name

    ! Writes the data to output file
    write(ts_unit,*) data
  
  end subroutine IO_write_dataset_vector
  
  
  subroutine IO_finish_galaxy(gal_id, info)
    ! Finishes IO (closes output files)
    integer, intent(in) :: gal_id, info
    close(ts_unit)
    close(columns_unit)
    
  end subroutine IO_finish_galaxy
  
  
  subroutine IO_end()  
    ! Finishes IO (in this case, closes output file)
    close(diag_unit)
    

  end subroutine IO_end  
  
end module IO