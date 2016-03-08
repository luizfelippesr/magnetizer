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
  double precision, allocatable, dimension(:,:) :: f_old, dfdt_old
  
  public dynamo_run

  contains
    subroutine dynamo_run(info, gal_id, flag, test_run)
      integer, intent(in) :: info, gal_id
      logical, intent(in) :: test_run
      integer, intent(out) :: flag
      logical :: ok, able_to_construct_profiles
      integer :: fail_count, i
      character(len=8) :: frmt
      character(len=8) :: gal_id_string
      integer, parameter :: MAX_FAILS=30
      double precision, dimension(nx) :: Btmp
      double precision :: this_t

      
      frmt='(I8.8)'
      write(gal_id_string,frmt) gal_id

      call cpu_time(cpu_time_start)
      
      ! Sets the number of variables
      call init_var
      ! Allocates f-array (which contains all the data for the calculations)
      if (.not. allocated(f)) allocate(f(nx,nvar)) 
      if (.not. allocated(dfdt)) allocate(dfdt(nx,nvar)) 
      if (.not. allocated(f_old)) allocate(f_old(nx,nvar))
      if (.not. allocated(dfdt_old)) allocate(dfdt_old(nx,nvar)) 
      
      f=0
      dfdt=0
      
      print *, 'Galaxy ',gal_id_string
      ! Prepares the grid where the computation will be made
      call construct_grid()
      ! Reads in the model parameters (for the first snapshot)
      call set_input_params(gal_id_string,info)  
      ! Sets other necessary parameters
      call set_calc_params  
      ! Constructs galaxy model for the initial snapshot
      able_to_construct_profiles = construct_profiles()
      ! Adds a seed field to the f-array (uses the profile info)
      call init_seed(f)
      ! Calculates |Bz| (uses profile info)
      call estimate_Bzmod(f)
      ! Backs up initial state
      f_old = f
      dfdt_old = dfdt ! this one is probably not necessary...
      
      ! Loops through the SAM's snapshots
      do it=1,n1
        if (info>0) print *, 'Main loop: Galaxy ',gal_id_string, ' it=',it

        ! Initializes the number of steps to the global input value
        nsteps = nsteps_0
        call set_ts_params()
        ! Constructs galaxy model for the present snapshot
        if (it /= 1) then
          if (p_oneSnaphotDebugMode) exit
          Btmp = sqrt(f(:,2)**2 + f(:,1)**2 + Bzmod**2)
          able_to_construct_profiles = construct_profiles(Btmp)
        endif
        ! If unable to construct the profiles, exit and write output
        if (.not. able_to_construct_profiles) then
          last_output = .true.
          call make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
          exit
        endif

        ! Will try to solve the equations a few times, with different
        ! timestep choices, if there is no success, aborts.
        do fail_count=0, MAX_FAILS
          ok = .true.
          ! If a run without magnetic field evolution was requested
          if (test_run) exit

          ! Loops through the timesteps
          do jt=1,nsteps
            this_t = t_Gyr +dt*t0_Gyr*jt
            if (info>2) then
              print *, 'Inner loop: Galaxy ',gal_id_string, ' jt=',jt, &
                      'nsteps=',nsteps, ' fail_count=',fail_count
              print *, '            t (Gyr)',t_Gyr +dt*t0_Gyr*jt
              print *, '            t',t
            endif

            call estimate_Bzmod(f)
            ! CALCULATE MAGNETIC ENERGY (WITHOUT THE FACTOR 1/(8PI))
            Btmp = sqrt(f(:,2)**2 + f(:,1)**2 + Bzmod**2)

            ! Updates all profiles
            able_to_construct_profiles = construct_profiles(sqrt(Btmp))
            ! If not able to construct profiles, flags and exit the loop
            if (.not.able_to_construct_profiles) then
              ok = .false.
              exit
            endif

            ! Runs Runge-Kutta time-stepping routine
            call rk(f, dfdt)
            ! Impose boundary conditions
            call impose_bc(f)
            ! If the magnetic field blows up, flags and exits the loop
            if (isnan(f(nxghost+1+nxphys/2,1)) .or. maxval(f)>1.d10) then
              ok = .false.
              ! Stores the bogus profiles
              ! NB the magnetic field info is out of date
              call make_ts_arrays(it+1,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
              exit
            endif

            if (p_oneSnaphotDebugMode) then
              ! Then, stores simulation output in arrays containing the time series
              call make_ts_arrays(jt,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
            endif
          end do ! timesteps loop
          
          ! If the time evolution was solved correctly (no blow-ups), exits
          ! Also exits if in the one Snaphot Debug Mode (we don't want to
          ! increase nsteps in this case).
          if (ok .or. p_oneSnaphotDebugMode) exit
          
          ! Otherwise, the calculation needs to be remade with a smaller
          ! timestep
          if (info>1) then
            print *, 'Galaxy ',gal_id_string, &
            ': NANs or unrealistically large magnetic fields detected, changing time step'
          endif
          
          ! Doubles the number of timesteps
          nsteps = nsteps * 2
          
          ! Resets the f array
          f = f_old
          dfdt = dfdt_old

          call set_ts_params()
        
        end do ! try or fail loop
        
        if (ok) then
          ! If the run was sucessful prepares to store and store
          f_old = f
          dfdt_old = dfdt
          ! Impose boundary conditions before writing output
          call impose_bc(f)  
          ! Estimates the value of |B_z| using Div B =0 condition
          call estimate_Bzmod(f)  

          ! Then, stores simulation output in arrays containing the time series
          if (p_oneSnaphotDebugMode) then
            i = jt
          else
            i = it
          endif
          call make_ts_arrays(i,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
        else
!           print *, 'FAIL'
!           stop
          ! If the run was unsuccessful, leaves the ts arrays with INVALID 
          ! values, except for ts_t (so that one knows when did things break)
          ts_t_Gyr(it) = t_Gyr !lfsr: actually, this very ugly... this variable 
                           !shouldn't be accessible at all!
          ! Resets the f array and adds a seed field 
          dfdt = 0.0  ! lfsr: Ugly as well... 
          f = 0.0  ! lfsr: Same here...
          call init_seed(f)
          ! Calculates |Bz|
          call estimate_Bzmod(f)
        endif

        ! Breaks loop if there are no more snapshots in the input
        if (last_output .or. p_oneSnaphotDebugMode) exit
        
        ! Reads in the model parameters for the next snapshot
        call set_input_params(gal_id_string,info)  
      end do  ! snapshots loop
      
      call write_output(gal_id)  !Writes final simulation output
      
      call cpu_time(cpu_time_finish)
      
      if (info>0) then
        print *, 'Galaxy ', gal_id_string,': finished after', &
            (cpu_time_finish - cpu_time_start), ' s  CPU time'
      endif
      flag=-1  !Tell calldynamo that simulation was successful
      call reset_input_params()  !Reset iread
      ! Resets the arrays which store the time series
      call reset_ts_arrays()  
    end subroutine dynamo_run
end module dynamo
!*****************************************************
