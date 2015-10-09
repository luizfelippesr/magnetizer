!*****************************************************
!This is the top level dynamo program.
!It contains the subroutine dynamo_run, which solves the dynamo equations.
!*****************************************************
module dynamo
  use ts_arrays
  use timestep
  use output
  use initial_conditions

  implicit none 
  private 
  
  integer :: it=0, jt=0
  double precision :: cpu_time_start, cpu_time_finish
  double precision, allocatable, dimension(:,:) :: f, dfdt
  
  public dynamo_run
  
  contains
    subroutine dynamo_run(info, gal_id, flag)
      integer, intent(in) :: info, gal_id
      integer, intent(out) :: flag
      logical :: ok
      integer :: fail_count
      character(len=8) :: frmt
      character(len=8) :: gal_id_string
      integer, parameter :: MAX_FAILS=6

      
      frmt='(I8.8)'
      write(gal_id_string,frmt) gal_id

      call cpu_time(cpu_time_start)
      
      ! Sets the number of variables
      call init_var
      ! Allocates f-array (which contains all the data for the calculations)
      if (.not. allocated(f)) allocate(f(nx,nvar)) 
      if (.not. allocated(dfdt)) allocate(dfdt(nx,nvar)) 
      
      print *, 'Galaxy ',gal_id_string
      ! Prepares the grid where the computation will be made
      call construct_grid()
      ! Reads in the model parameters (for the first snapshot)
      call set_input_params(gal_id_string,info)  
      ! Sets other necessary parameters
      call set_calc_params  
      ! Constructs galaxy model for the initial snapshot
      call construct_profiles()
      ! Adds a seed field to the f-array (uses the profile info)
      call init_seed(f)
      ! Calculates |Bz|
      call estimate_Bzmod(f)
      

      ! Loops through the SAM's snapshots
      do it=1,n1  
        print *, 'Galaxy ',gal_id_string, ' it=',it

        ! Initializes the number of steps to the global input value
        nsteps = nsteps_0 
        ! Constructs galaxy model for the present snapshot
        if (it /= 1) call construct_profiles()
          
        ! Will try to solve the equations a few times, with different 
        ! timestep choices, if there is no success, aborts.
        do fail_count=0, MAX_FAILS
          ok = .true.
          
          ! Loops through the timesteps
          do jt=1,nsteps
            print *, 'Galaxy ',gal_id_string, ' jt=',jt
            ! Runs Runge-Kutta time-stepping routine
            call rk(f)  
            
            ! If the magnetic field blows up or something else blows up
            ! flags and exit the loop
            if (isnan(f(4+nxphys/2,1)) .or. maxval(f)>1.d10) then  
              ok = .false.
              exit 
            endif
          end do ! timesteps loop
          
          ! If the time evolution was solved correctly (no blow-ups), exits
          if (ok) exit
          ! Otherwise, the calculation needs to be remade with a smaller
          ! timestep
          if (info>0) then
            print *, 'Galaxy ',gal_id_string, &
            ': NANs or unrealistically large magnetic fields detected, changing time step'
            print *, 'nsteps',nsteps
          endif
          
          ! Doubles the number of timesteps
          nsteps = nsteps * 2
          call set_ts_params()
        
        end do ! try or fail loop
        
        if (ok) then
          ! If the run was sucessful prepares to store and store
          
          ! Impose boundary conditions before writing output
          call impose_bc(f)  
          ! Estimates the value of |B_z| using Div B =0 condition
          call estimate_Bzmod(f)  
          ! Then, stores simulation output in arrays containing the time series
          call make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)  
        else
          ! If the run was unsuccessful, leaves the ts arrays with INVALID 
          ! values, except for ts_t (so that one knows when did things break)
          ts_t(it) = t !lfsr: actually, this very ugly... this variable 
                       !shouldn't be accessible at all!
          dfdt = 0.0  ! Same here 
          
          ! Resets the f array and adds a seed field 
          call init_seed(f)
          ! Calculates |Bz|
          call estimate_Bzmod(f)
        endif

        ! Breaks loop there are no more snapshots in the input
        if (last_output) exit        
        
      end do  ! snapshots loop
      
      call write_output(gal_id)  !Writes final simulation output
      
      call cpu_time(cpu_time_finish)
      
      if (info>0) then
        print *, 'Galaxy ', gal_id_string,': finished after', &
            (cpu_time_finish -cpu_time_start), ' s  CPU time'
      endif

      open(20,file= 'diagnostic.out',status="old",position="append")
      write(20,*) "simulation time in seconds: ", (cpu_time_finish -cpu_time_start)
      close(20)

      flag=-1  !Tell calldynamo that simulation was successful
      iread=0  !Reset iread
    end subroutine dynamo_run
end module dynamo
!*****************************************************
