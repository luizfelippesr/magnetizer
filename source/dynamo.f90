!*****************************************************
!This is the top level dynamo program.
!It contains the subroutine dynamo_run, which solves the dynamo equations.
!*****************************************************
module dynamo
  use ts_arrays
  use timestep
  use output
  use initial_conditions
  use messages
  implicit none 
  private 
  
  integer :: it=0, jt=0
  double precision :: cpu_time_start, cpu_time_finish
  double precision, allocatable, dimension(:,:) :: f, dfdt
  double precision, allocatable, dimension(:,:) :: f_old, dfdt_old

  public dynamo_run

  contains
    subroutine dynamo_run(gal_id, test_run, rank)
      integer, intent(in) :: gal_id
      integer, intent(in), optional :: rank
      logical, intent(in) :: test_run
      logical :: ok, able_to_construct_profiles, elliptical
      integer :: fail_count, rank_actual
      integer, parameter :: MAX_FAILS=3
      double precision, dimension(nx) :: Btmp
      double precision :: this_t

      if (present(rank)) then
          rank_actual = rank
      else
          rank_actual = 0
      endif

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

      ! Prepares the grid where the computation will be made
      call construct_grid()
      ! Reads in the model parameters (for the first snapshot)
      call set_input_params(gal_id)
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
      do it=init_it,n1
        this_t = t_Gyr
        call message('Main loop: it = ', gal_id=gal_id, val_int=it, info=1)

        ! Initializes the number of steps to the global input value
        nsteps = nsteps_0

        ! Traps the case of negligible disks
        ! (a recent major merger can convert the galaxy into an elliptical,
        !  making Mdisk=rdisk=0)
        if (r_disk < r_max_kpc*rmin_over_rmax .or.  &
              Mgas_disk < Mgas_disk_min) then
            elliptical = .true.
            ! Resets the f array and adds a seed field
            dfdt = 0.0
            f = 0.0
            call init_seed(f)
        else
            ! If galaxy was an elliptical during previous snapshot,
            ! reconstruct the profiles
            if (elliptical) then
              able_to_construct_profiles = construct_profiles()
              call init_seed(f)
              call estimate_Bzmod(f)
            endif
            elliptical = .false.
        endif

        call set_ts_params()
        ! Constructs galaxy model for the present snapshot
        if (it /= init_it) then
          if (p_oneSnaphotDebugMode) exit
          if (p_simplified_pressure) then
            able_to_construct_profiles = construct_profiles()
            ! If unable to construct the profiles, exit and write output
            if (.not. able_to_construct_profiles) then
              last_output = .true.
              call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
              exit
            endif
          endif
        endif

        ! Will try to solve the equations a few times, with different
        ! timestep choices, if there is no success, aborts.
        do fail_count=0, MAX_FAILS
          ok = .true.
          ! If a run without magnetic field evolution was requested
          if (test_run .or. elliptical) exit

          ! Loops through the timesteps
          do jt=1,nsteps
            this_t = t_Gyr +dt*t0_Gyr*jt
            call message('Inner loop: jt = ',gal_id=gal_id, val_int=jt, info=2)
            call message('Inner loop: fail_count =',val_int=fail_count,gal_id=gal_id, info=2)
            call message('Inner loop: nsteps = ',gal_id=gal_id, val_int=nsteps, info=3)
            call message('Inner loop: t (Gyr) = ', t_Gyr +dt*t0_Gyr*jt,gal_id=gal_id, info=3)

            call estimate_Bzmod(f)
            ! If not using the simplified pressure, all profiles need to be
            ! recomputed at each timestep, since they depend on the large
            ! scale field.
            if (.not.p_simplified_pressure) then
              ! Computes the total magnetic field at this timestep
              Btmp = sqrt(f(:,2)**2 + f(:,1)**2 + Bzmod**2)
              ! Updates all profiles
              able_to_construct_profiles = construct_profiles(sqrt(Btmp))
              ! If not able to construct profiles, flags and exit the loop
              if (.not.able_to_construct_profiles) then
                call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
                ok = .false.
                exit
              endif
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
              call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
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
          call message('NANs or unrealistically large magnetic fields detected, changing time step', &
                         gal_id=gal_id, info=2)
          
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
          call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
        else
          ! If the run was unsuccessful, leaves the ts arrays with INVALID 
          ! values, except for ts_t (so that one knows when did things break)
          if (.not.p_oneSnaphotDebugMode) &
            call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k, &
                              alp,Uz,Ur,n,Beq,rmax,delta_r, invalid_run=.true.)
          ! Resets the f array and adds a seed field 
          dfdt = 0.0
          f = 0.0
          call init_seed(f)
          ! Calculates |Bz|
          call estimate_Bzmod(f)
        endif

        ! Breaks loop if there are no more snapshots in the input
        if (last_output .or. p_oneSnaphotDebugMode) exit
        
        ! Reads in the model parameters for the next snapshot
        call set_input_params(gal_id)
      end do  ! snapshots loop
      

      call cpu_time(cpu_time_finish)


      call message('Finished after ', (cpu_time_finish - cpu_time_start),  &
                   gal_id= gal_id, msg_end=' s  CPU time', info=1)

      call write_output(gal_id, ok, cpu_time_finish - cpu_time_start)  !Writes final simulation output

      call reset_input_params()  !Reset iread
      ! Resets the arrays which store the time series
      call reset_ts_arrays()  
    end subroutine dynamo_run
end module dynamo
!*****************************************************
