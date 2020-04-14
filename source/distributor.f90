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
module jobs
  use mpi
  use input_params
  use global_input_parameters
  use IO
  use data_transfer
  use messages
  implicit none
  private
  public jobs_prepare, jobs_distribute, jobs_reads_galaxy_list

  ! Module variables
  integer, parameter :: master_rank = 0
  logical :: select_completed = .false.
  logical :: select_incomplete = .true.
  integer :: rank, nproc, info_mpi
  logical :: lresuming_run = .false.
  double precision :: tstart,tfinish
  character(len=100) :: date

  ! Function to be passed as argument
  abstract interface
    function work_routine_template(gal_id, switch, error)
      integer, intent(in) :: gal_id
      logical, intent(in) :: switch
      logical, intent(out) :: error
      double precision :: work_routine_template
    end function work_routine_template
  end interface
  ! Function to be passed as argument
  abstract interface
    subroutine write_routine_template(gal_id, runtime)
      integer, intent(in) :: gal_id
      double precision, intent(in) :: runtime
    end subroutine write_routine_template
  end interface

contains
  function jobs_reads_galaxy_list(list_args_position) result(galaxies)
    integer, dimension(:), allocatable :: galaxies
    logical, dimension(ngals) :: valid
    integer,intent(in) :: list_args_position
    character(len=100) :: command_argument
    integer :: len_arg, i_arg, i

    i_arg = list_args_position
    valid = .false.
    allocate(galaxies(ngals))

    call get_command_argument(i_arg, command_argument)
    len_arg = len_trim(command_argument)

    if (len_arg == 0) then
      ! No galaxies list present!
      do i=1,ngals
        galaxies(i)=i
      enddo
    else

      do i=1,ngals
        galaxies(i) = str2int(command_argument)
        valid(i) = .true.

        i_arg = i_arg+1
        call get_command_argument(i_arg, command_argument)
        if (len_trim(command_argument)==0) exit
      enddo
      galaxies = pack(galaxies,valid)
    end if

  end function jobs_reads_galaxy_list

  subroutine jobs_prepare(completed, incomplete)
    ! Initializes MPI and sets the scene for distributing work between processors
    logical, intent(in), optional :: incomplete, completed
    character(len=100) :: command_argument

    integer :: ierr, rc
    integer,dimension(8) :: time_vals

    if (present(incomplete)) select_incomplete = incomplete
    if (present(completed)) select_completed = completed

    call MPI_INIT(ierr)
    if (ierr/= MPI_SUCCESS) then
      call message('Error starting MPI program. Terminating.')
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
    endif
    !Get the number of processors this job
    call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
    !Get the rank of the processor this thread is running on
    call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
    if (ierr/= MPI_SUCCESS) then
      call message('Error getting processor rank. Terminating.')
      call MPI_Abort(MPI_COMM_WORLD, rc, ierr)
    endif


    if (nproc==1) then
      call welcome_message()
      call message(' ', rank=-1)
      call message('Runnning on a single processor')
    else
      if (rank==master_rank) call welcome_message()
      call message(' ', rank=rank, master_only=.true.)
      call message('Runnning on', val_int=nproc, msg_end='processors', &
                   master_only=.true.)
    endif

    call get_command_argument(1, command_argument)
    if (len_trim(command_argument) == 0) then
      ! Uses example parameter file if nothing was found
      command_argument = 'example/example_global_parameters.in'
      call read_global_parameters(trim(command_argument))
      call message('No parameter file provided. Using standard: '// &
                    trim(command_argument), master_only=.true.,&
                    set_info=info)
    else
      ! Uses specified parameter file
      call read_global_parameters(trim(command_argument))
      call message('Using global parameters file: '// trim(command_argument), &
                  master_only=.true., set_info=info)
    endif

    if (rank == master_rank) then !Only the master (rank 0)
      tstart = MPI_wtime()
    else
      tstart = -1
      tfinish = -1
    endif
    ! Checks whether this is a new run or if one is resuming an older run
    ! Note: resuming runs currently only works for separate output files
    if (p_IO_separate_output) then
      if (rank == master_rank) then
        inquire (file=trim(output_file_name), exist=lresuming_run)
      endif
      ! Broadcast the the result
      call MPI_Bcast(lresuming_run, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    endif

    if (rank==master_rank) then
      call date_and_time(values=time_vals)
      ! Displays date in the format yyyy-mm-dd  hh-mm-ss timezone
      date = str(time_vals(1))//'-'//str(time_vals(2))//'-'//str(time_vals(3)) &
                            //'  '//str(time_vals(5))//':'//str(time_vals(6))&
                            //' (UTC'//str(time_vals(4))//')'
    endif
    ! It is probably a good idea to sync, to avoid different nodes having
    ! slightly different times, for any reason (the attribute creation is
    ! collective, anyway).
    call MPI_Bcast(date, len(date), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Avoids MPI errors
    call MPI_Info_create(info_mpi, IERR)
    call MPI_Info_set(info_mpi, "romio_ds_write", "disable", ierr)
    call MPI_Info_set(info_mpi, "romio_ds_read", "disable", ierr)

    ! Prepares IO (reads ngals, nz's from the hdf5 input file)
    call IO_prepare(MPI_COMM_WORLD, info_mpi, rank, lresuming_run, date)

  end subroutine jobs_prepare

  subroutine jobs_distribute(work_routine, write_routine, func_flag, galaxies_list)
    procedure(work_routine_template) :: work_routine
    procedure(write_routine_template) :: write_routine
    logical, intent(in) :: func_flag
    integer, dimension(:), intent(in) :: galaxies_list
    integer :: igal, nmygals, igal_first, igal_last, igal_finished
    integer, allocatable, dimension(:) :: mygals, allgals
    integer :: i, j, iproc
    logical :: lstop, error
    logical :: start_galaxy = .false.
    integer :: ierr, ncycles, flush_signal, ngals
    integer, parameter :: finished_tag = 0
    integer, parameter :: newjob_tag = 17
    integer, dimension(MPI_STATUS_SIZE) :: status
    double precision :: runtime
    ngals = size(galaxies_list)

    ! Allocates arrays to for logging
    nmygals = 0
    allocate(mygals((ngals+1)*nproc))
    mygals = -17 ! Initializes to bad value
    allocate(allgals((ngals+1)*nproc))
    flush_signal = ngals+42 ! Arbitrary larger-than-ngals value


    call message('Warm up run', master_only=.true.)
    if (rank==master_rank) then
      ! Initializes IO (this also reads ngals from the hdf5 input file)
      call IO_start(MPI_COMM_WORLD, info_mpi, rank, lresuming_run, date)

      ! Before distributing work among workers, this warm-up runs a single galaxy
      ! on the master process.
      ! The main reason is setting the scene for the HDF5 library:
      ! this creates or prepare the file, allowing worker processes to be able
      ! to access it (in read only mode).
      ! Also, this helps catching any obvious problem in the beginning
      do igal=1, ngals
        start_galaxy = IO_start_galaxy(galaxies_list(igal))
        if (decide_run(start_galaxy)) then
          runtime = work_routine(galaxies_list(igal), func_flag, error)
          if (runtime>0) call write_routine(galaxies_list(igal), runtime)
          if (error) cycle
          if (rank == master_rank) then
            nmygals = nmygals + 1
            mygals(nmygals) = galaxies_list(igal)
          endif
          exit
        endif
      enddo
      call IO_end()
      lresuming_run = .true.
      ! Now the run exists!
    endif
    call MPI_Bcast(lresuming_run, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call message('Warm up done', master_only=.true.)

    call IO_start(MPI_COMM_WORLD, info_mpi, rank, lresuming_run, date)

    ! Calculates the number of cycles
    ncycles = max(ngals/p_ncheckpoint,1)
    call message('Total number of cycles',val_int=ncycles, master_only=.true.)

    if (nproc<=1) then
      p_master_works_too=.true.
      p_master_skip=1
      nproc = 1
    endif

    ! ----------
    !   Master
    ! ----------
    if (rank==master_rank .and. ngals/=1) then
      do j=0, ncycles
        ! Finds boundaries of present cycle
        igal_first = 1+j*p_ncheckpoint
        igal_last = min((j+1)*p_ncheckpoint, ngals)
        call message('Cycle',val_int=j+1, master_only=.true.)

        ! Submit initial jobs
        do iproc=1, nproc-1
            igal = iproc-1+igal_first
            if (igal > igal_last) igal=ngals+1
            call MPI_Send(igal, 1, MPI_INTEGER, iproc, newjob_tag, MPI_COMM_WORLD, ierr)
        enddo

        ! Loops sending and receiving
        do igal=nproc+igal_first-1, igal_last
          if ( p_master_works_too .and. &
              modulo(igal, (nproc + int(nproc/p_master_skip))-1) == 0) then
            ! Call dynamo code for each galaxy in the list
            call message('Starting',gal_id=galaxies_list(igal), rank=rank)

            ! Check whether the galaxy was processed previously
            start_galaxy = IO_start_galaxy(galaxies_list(igal))
            if (decide_run(start_galaxy)) then
              ! If it is a new galaxy, runs it!
              runtime = work_routine(galaxies_list(igal), func_flag, error)
              if (runtime>0d0) call write_routine(galaxies_list(igal), runtime)
              nmygals = nmygals + 1
              mygals(nmygals) = galaxies_list(igal)
            else
              call message('Skipping',gal_id=galaxies_list(igal), rank=rank)
            endif
          else
            ! Receives finished work from worker
            call MPI_Recv(igal_finished, 1, MPI_INTEGER, MPI_ANY_SOURCE, finished_tag, &
                          MPI_COMM_WORLD, status, ierr)
            ! Finds out which worker has just finished this work
            iproc = STATUS(MPI_SOURCE)
            ! Receives the data and saves to disk
            if (igal_finished>0) call receive_and_save(galaxies_list(igal_finished), iproc)
            ! Sends new job to that worker
            call MPI_Send(igal, 1, MPI_INTEGER, iproc, newjob_tag, MPI_COMM_WORLD, &
                          ierr)
          endif
        enddo

        ! Checks whether a STOP file exists
        inquire (file='STOP', exist=lstop)
        ! If yes, exits gently
        if (lstop) then
          call message('STOP file found. Will finish this cycle and then exit.', &
                      info=0, master_only=.true.)
        endif
        if (p_max_walltime>0) then
          if ((MPI_wtime()-tstart) > p_max_walltime) then
            lstop = .true.
            call message('Maximum walltime reached! Will finish this cycle and then exit.', &
                      info=0, master_only=.true.)
          endif
        end if

        ! Goes through all workers, requesting them to either stop or flush
        ! the data to the disk
        do iproc=1, nproc-1
            ! Receives whatever message it is sending
            call MPI_Recv(igal_finished, 1, MPI_INTEGER, iproc, finished_tag, &
                          MPI_COMM_WORLD, status, ierr)
            ! Receives the data and saves to disk
            if (igal_finished>0) call receive_and_save(galaxies_list(igal_finished), iproc)

            if (lstop .or. (igal_last>=ngals)) then
              ! If it is time to stop, tells to workers to finish
              call MPI_Send(-1, 1, MPI_INTEGER, iproc, newjob_tag, &
                            MPI_COMM_WORLD, ierr)
            else
              ! Otherwise, sends a signal requesting the worker to flush
              call MPI_Send(flush_signal, 1, MPI_INTEGER, iproc, newjob_tag, &
                            MPI_COMM_WORLD, ierr)
            endif
        enddo

        ! Tells the master to stop
        if (lstop .or. igal_last>=ngals) exit
        call IO_flush() ! Tells the master to checkpoint
      enddo
    ! ----------
    !   Worker
    ! ----------
    else if (ngals/=1) then
        do
          ! Receives new task or a flag
          call MPI_Recv(igal, 1, MPI_INTEGER, 0, newjob_tag, MPI_COMM_WORLD, status, ierr)

          ! If received an exit flag, exits
          if (igal<0) then
            call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
            call message('Exiting', rank=rank)
            exit
          else if (igal<=ngals) then
            ! Call dynamo code for each galaxy in mygals
            call message('Starting',gal_id=galaxies_list(igal), rank=rank, &
                         info=3)
            ! Check whether the galaxy was processed previously
            start_galaxy = IO_start_galaxy(galaxies_list(igal))
            if (decide_run(start_galaxy)) then
              ! If it is a new galaxy, runs it!
              runtime = work_routine(galaxies_list(igal), func_flag, error)
              nmygals = nmygals + 1
              mygals(nmygals) = galaxies_list(igal)
              ! Sends the result (to mark it done)
              if (runtime>0d0) then
                call message('Done. Sending ',gal_id=galaxies_list(igal), &
                             rank=rank, info=3)
                call MPI_Send(igal, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
                call write_routine(galaxies_list(igal), runtime)
              else
                call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
              endif
            else
              call message('Skipping',gal_id=galaxies_list(igal), rank=rank)
              ! Sends the result (to mark it done)
              call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
            endif
          else if (igal==flush_signal) then
            call IO_flush() ! Tells the worker to checkpoint
            call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
          else
            ! Invalid galaxy. Probably, more processors than galaxies!) Nothing is done.
            call MPI_Send(-1, 1, MPI_INTEGER, 0, finished_tag, MPI_COMM_WORLD, ierr)
          endif
      enddo
    endif

    call IO_flush()
    call message('All computing done', info=0)

    if (nmygals>0) then
      if (nmygals<10) then
        print *, trim(str(rank)),':  Finished working on galaxies:', mygals(:nmygals),'.'
      else
        call message('Finished working on', val_int=nmygals, msg_end='galaxies.')
      endif
    else
      call message('Finished without working on any galaxy.')
    endif

    call MPI_Gather(mygals, ngals, MPI_INTEGER, allgals, ngals, &
                    MPI_INTEGER, master_rank, MPI_COMM_WORLD, ierr)

    ! Gets the date
    if (rank == master_rank) then
      call date_and_time(date=date)
    end if
    ! Broadcast the date
    call MPI_Bcast(date, 8, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Prints a small report
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! Finalizes IO
    call IO_end(date)
    call message('', master_only=.true.)
    call message('IO finished', master_only=.true.)

    if (rank == master_rank) then !Only the master (rank 0)
      tfinish= MPI_WTime()
      ! Removes stop file if necessary
      if (lstop) then
        i = 17
        open(unit=i, FILE='STOP')
        close(i, status='delete')
        call message('Removed STOP file.', master_only=.true.)
      endif
    endif

    if (info>1 .and. rank==master_rank) then
      i = 1
      do j=1,ngals*nproc
        igal = allgals(j)
        if (igal>0) then
          mygals(i) = igal
          i = i+1
        endif
      enddo
      ngals = i-1
      if (ngals<200 .and. ngals>1) then
        print *,
        print *, '  Galaxies in this run:'
        print *, mygals(:ngals)
        print *,
      else
        print *,
        print *, '  Number of galaxies in this run:', ngals
        print *,
      endif
    endif

    call message('Total wall time in seconds =',tfinish-tstart, &
                master_only=.true., info=0)
    if ((ngals>1) .and. (ngals /= 0)) then
      call message('Wall time per galaxy =', (tfinish-tstart)/ngals, &
                  master_only=.true., info=0)
      call message('Average CPU per galaxy =', (tfinish-tstart)*nproc/ngals, &
                  master_only=.true., info=0)
    endif
    !Tell the MPI library to release all resources it is using
    call MPI_Finalize(ierr)
    call message('MPI finished', master_only=.true.)
  end subroutine jobs_distribute

  logical pure function decide_run(start_galaxy)
    logical, intent(in) :: start_galaxy
    decide_run = ((start_galaxy .and. select_incomplete) .or. &
                 ((.not.start_galaxy) .and. select_completed))
  end function decide_run

  subroutine welcome_message()
    print *,"   __  __                        _   _               "
    print *,"  |  \/  | __ _  __ _ _ __   ___| |_(_)_______ _ __  "
    print *,"  | |\/| |/ _` |/ _` | '_ \ / _ \ __| |_  / _ \ '__| "
    print *,"  | |  | | (_| | (_| | | | |  __/ |_| |/ /  __/ |    "
    print *,"  |_|  |_|\__,_|\__, |_| |_|\___|\__|_/___\___|_|    "
    print *,"               |___/                                 "
  end subroutine

end module jobs
