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
  double precision, allocatable, dimension(:,:) :: f
  double precision, allocatable, dimension(:,:) :: f_snapshot_beginning

  public dynamo_run

  contains
    subroutine dynamo_run(gal_id, test_run, rank)
      use interpolation
      double precision :: cpu_time_start
      integer, intent(in) :: gal_id
      integer, intent(in), optional :: rank
      logical, intent(in) :: test_run
      logical :: ok, able_to_construct_profiles, elliptical, timestep_ok, error
      logical :: next_output
      integer :: fail_count, rank_actual
      double precision, dimension(nx) :: Btmp
      double precision :: this_t, t_last_sign_choice
      elliptical = .false. ;  ok = .true.

      t_last_sign_choice = 0d0
      ! Initializes the floor sign randomly
      call random_number(Bfloor_sign)
      Bfloor_sign = sign(1d0,Bfloor_sign-0.5d0)

      if (present(rank)) then
          rank_actual = rank
      else
          rank_actual = 0
      endif

      call cpu_time(cpu_time_start)

      ! Sets the number of variables
      call init_var

      ! Resets the status code
      call reset_status_code(test_run)

      ! Reads in the model parameters (for the first snapshot)
      call set_input_params(gal_id, error)
      if (error) then
        return
      endif
      call construct_grid(r_disk, r_max_kpc_history)
      ! Allocates f-array (which contains all the data for the calculations)
      call check_allocate_f_array(f, nvar)
      call check_allocate(alp); alp = 0.0
      call check_allocate(Bzmod); Bzmod = 0.0

      ! Sets other necessary parameters
      call set_calc_params
      ! Constructs galaxy model for the initial snapshot
      able_to_construct_profiles = construct_profiles()

      if (.not.able_to_construct_profiles) then
        call message('Could not construct profiles for this galaxy.', &
                     gal_id=gal_id, info=1)
        call write_and_finish(cpu_time_start, gal_id, this_t)
        return
      endif

      if (r_disk > p_rdisk_min .and.  Mgas_disk > Mgas_disk_min) then
        ! Adds a seed field to the f-array (uses the profile info)
        call init_seed(f)
        ! Calculates |Bz| (uses profile info)
        call estimate_Bzmod(f)
      endif

      ! Loops through the SAM's snapshots
      do it=init_it,max_it
        this_t = t_Gyr
        call message('Main loop: it = ', gal_id=gal_id, val_int=it, info=2)

        ! Traps the case of negligible disks
        ! (a recent major merger can convert the galaxy into an elliptical,
        !  making Mdisk=rdisk=0)
        if (r_disk < p_rdisk_min .or.  Mgas_disk < Mgas_disk_min) then
            elliptical = .true.
            ! Resets the f array
            f = 0d0
            !! Adds a seed field
            !call init_seed(f)
        else
            ! If galaxy was an elliptical during previous snapshot,
            ! reconstruct the profiles
            if (elliptical) then
              able_to_construct_profiles = construct_profiles()
              call init_seed(f)
              call estimate_Bzmod(f)
            endif
            if (able_to_construct_profiles) then
              elliptical = .false.
            endif
        endif

        next_output = .false.
        ! Constructs galaxy model for the present snapshot
        if (it /= init_it .and. .not.elliptical) then
          if (p_oneSnaphotDebugMode) exit

          if (p_simplified_pressure) then
            call adjust_grid(f, r_disk)
            ! Impose boundary conditions
            call impose_bc(f)
            able_to_construct_profiles = construct_profiles()
            ! If unable to construct the profiles, write output and exit loop
            if (.not. able_to_construct_profiles) then
!                 next_output = .true.
!               call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,&
!                                   alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
              last_output = .true.
              exit
            endif
          endif
        endif

        ! Sets the timestep
        timestep_ok = set_timestep(h, v, etat)
        if (.not. timestep_ok) then
          call message('Timestep became very small! Using maximum number of' &
                      //' timesteps instead.', gal_id=gal_id, info=1)
        endif

        ! Backs up state at the beginning of the snapshot
        call check_allocate_f_array(f_snapshot_beginning, nvar)
        f_snapshot_beginning = f

        ! Will try to solve the equations a few times, with different
        ! timestep choices, if there is no success, aborts.
        do fail_count=0, p_MAX_FAILS
          ok = .true.
          ! If a run without magnetic field evolution was requested or
          ! if it is an elliptical galaxy, exits before the magnetic field
          ! calculations (adjusts unused arrays for output)
          if (test_run .or. elliptical .or. next_output) then
            call check_allocate(alp); alp = 0.0
            call check_allocate(Bzmod); Bzmod = 0.0

            exit
          endif

          ! Loops through the timesteps
          do jt=1,nsteps
            this_t = t_Gyr +dt*t0_Gyr*jt
            call message('Inner loop: jt = ',gal_id=gal_id, val_int=jt, info=3)
            call message('Inner loop: fail_count =', val_int=fail_count, &
                         gal_id=gal_id, info=3)
            call message('Inner loop: nsteps = ', gal_id=gal_id, &
                         val_int=nsteps, info=4)
            call message('Inner loop: t (Gyr) = ', t_Gyr +dt*t0_Gyr*jt, &
                         gal_id=gal_id, info=4)

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
!                 call make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau, &
!                                     alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
                ok = .false.
                exit
              endif
            endif

            ! Alternate floor
            if (p_time_varying_floor) then
              if (this_t-t_last_sign_choice > minval(tau)*t0_Gyr) then
                call random_number(Bfloor_sign)
                Bfloor_sign = sign(1d0,Bfloor_sign-0.5d0)
                t_last_sign_choice = this_t
              endif
            endif

            ! Runs Runge-Kutta time-stepping routine
            call rk(f)

            ! Impose boundary conditions
            call impose_bc(f)
            ! If the magnetic field blows up, flags and exits the loop
            if (maxval(f(:,1:2))>1.d7) then
              ok = .false.
              ! Stores the bogus profiles
              ! NB the magnetic field info is out of date
              call make_ts_arrays(it,this_t,f,Bzmod,alp,rmax)
              exit
            endif

            if (p_oneSnaphotDebugMode) then
              ! Then, stores simulation output in arrays containing the time series
              call make_ts_arrays(jt,this_t,f,Bzmod,alp,rmax)
            endif
          end do ! timesteps loop
          
          ! If the time evolution was solved correctly (no blow-ups), exits
          ! Also exits if in the one Snaphot Debug Mode (we don't want to
          ! increase nsteps in this case).
          if (ok .or. p_oneSnaphotDebugMode) exit
          
          ! Otherwise, the calculation needs to be remade with a smaller
          ! timestep
          call error_message('dynamo_run', 'NANs or unrealistically large '&
                             //'magnetic fields detected, changing time step',&
                             gal_id=gal_id, info=2, code='m')
          
          ! Doubles the number of timesteps
          timestep_ok = set_timestep(reduce_timestep=.true.)
          if (.not. timestep_ok) then
            call error_message('dynamo_run', &
                               'Timestep became too small! Aborting...', &
                               gal_id=gal_id, info=1, code='t')
            ok = .false.
            exit
          endif

          ! Resets the f array to the initial state in the beginning of the
          ! snapshot
          f = f_snapshot_beginning
        end do ! try or fail loop

        if (.not.ok) then
          ! If the run was unsuccessful, leaves the ts arrays with INVALID
          ! values, except for ts_t (so that one knows when did things break)
          call error_message('dynamo_run','Invalid run! Aborting...', &
                             gal_id=gal_id, info=2, code='i')

          if (.not.p_oneSnaphotDebugMode) &
            call make_ts_arrays(it,this_t,f,Bzmod,alp,rmax)
          last_output = .true.
          ! Resets the f array and adds a seed field
          f = 0.0
          call init_seed(f)
        endif
        ! Impose boundary conditions to the final result
        call impose_bc(f)
        ! Estimates the value of |B_z| using Div B = 0 condition
        call estimate_Bzmod(f)

        ! Stores simulation output in arrays containing the time series
        call make_ts_arrays(it,this_t,f,Bzmod,alp,rmax)

        ! Breaks loop if there are no more snapshots in the input
        if (last_output .or. p_oneSnaphotDebugMode) exit

        ! Resets the status code
        call reset_status_code(test_run)

        ! Reads in the model parameters for the next snapshot
        call set_input_params(gal_id, error)
        if (error) exit

      end do  ! snapshots loop

      !Writes final simulation output
      call write_and_finish(cpu_time_start, gal_id, this_t)
    end subroutine dynamo_run

    subroutine write_and_finish(cpu_time_start, gal_id, this_t)
      double precision, intent(in) :: cpu_time_start, this_t
      integer, intent(in) :: gal_id
      double precision :: cpu_time_finish

      !Writes final simulation output
      call cpu_time(cpu_time_finish)

      call estimate_Bzmod(f)

      if (it == 0) it = init_it ! If a problem happened before entering the
                                ! snapshots loop, initialize it properly.
      call make_ts_arrays(it,this_t,f,Bzmod,alp,rmax)

      call write_output(gal_id, cpu_time_finish - cpu_time_start)
      call reset_input_params()  !Resets iread
      ! Resets the arrays which store the time series
      call reset_ts_arrays()
      call message('Finished after ', (cpu_time_finish - cpu_time_start),  &
                   gal_id= gal_id, msg_end='s  CPU time', info=1)

    end subroutine write_and_finish
end module dynamo
!*****************************************************
