! Contains two modules which deal with input parameters
module input_params
  ! Reads and sets input parameters
  use global_input_parameters
  use math_constants

  implicit none

  character (len=200) :: header_time_indep, header_time_dep
  integer :: iread=0

  !
  logical :: last_output = .false.

  ! Time-stepping parameters
  integer :: n1 = -1 !Total number of snapshots (set to -1 to flag it is uninitialized)
  integer :: init_it = -1 ! Index of the initial snapshot (set to -1 to flag it is uninitialized)
  double precision :: tsnap !Time between successive snapshots
  integer :: nsteps=20  !Number of timesteps in between snapshots
  double precision :: dt,t=0,first=0.d0 !Timestep variables (lfsr: a very bad place to define them, indeed)

  double precision, private :: time_between_inputs=0
  
  ! Galaxy parameters
  double precision :: t_Gyr
  double precision :: r_max_kpc, R_kappa
  double precision :: l_sol_kpc, r_l_kpc
  double precision :: v_sol_kms, r_v_kpc
  double precision :: Uz_sol_kms, r_Uz_kpc
  double precision :: h_sol_kpc, r_h_kpc
  double precision :: r_disk, v_disk
  double precision :: r_bulge, v_bulge
  double precision :: r_halo, v_halo, nfw_cs1 
  double precision :: Mstars_disk, Mgas_disk, SFR

  ! All galaxy data (private!)
  integer, private, parameter :: number_of_columns=11 ! Maximum umber of columns in the galaxy input files
  double precision, allocatable, dimension(:,:),private :: galaxy_data
  double precision, allocatable, dimension(:), private :: output_times
  integer, private ::  current_gal_id

  contains
   subroutine set_ts_params()
      use units

      tsnap = time_between_inputs/t0_Gyr
      dt = tsnap/nsteps  !Timestep in units of t0=h0^2/etat0
      print *, 'set_ts_params  ','dt=',dt*t0_Gyr,'Gyr', time_between_inputs/100.
      print *, nsteps
      
    endsubroutine set_ts_params


    subroutine read_input_parameters(gal_id)
      ! Reads the input parameters file to RAM
      use iso_fortran_env
      use IO
      integer, intent(in) :: gal_id
      integer :: i

      if (.not.allocated(galaxy_data)) then
        allocate(galaxy_data(number_of_redshifts,number_of_columns))
        allocate(output_times(number_of_redshifts))
        output_times = 0
        ! Reads the output times
        call IO_read_dataset_scalar('t', gal_id, info, output_times, &
                                    nrows=number_of_redshifts)
      endif

      ! Saves current galaxy identifier
      current_gal_id = gal_id
      ! Resets galaxy data array
      galaxy_data(:,:) = -1

      ! Reads time dependent parameters for this galaxy
      call IO_read_dataset_scalar('r_disk', gal_id, info, galaxy_data(:,1))
      call IO_read_dataset_scalar('v_disk', gal_id, info, galaxy_data(:,2))
      call IO_read_dataset_scalar('r_bulge', gal_id, info, galaxy_data(:,3))
      call IO_read_dataset_scalar('v_bulge', gal_id, info, galaxy_data(:,4))
      call IO_read_dataset_scalar('r_halo', gal_id, info, galaxy_data(:,5))
      call IO_read_dataset_scalar('v_halo', gal_id, info, galaxy_data(:,6))
      call IO_read_dataset_scalar('nfw_cs1', gal_id, info, galaxy_data(:,7))
      call IO_read_dataset_scalar('Mgas_disk', gal_id, info, galaxy_data(:,8))
      call IO_read_dataset_scalar('Mstars_disk', gal_id, info, galaxy_data(:,9))
      call IO_read_dataset_scalar('SFR', gal_id, info, galaxy_data(:,10))
      ! Reads the maximum radius for this galaxy
      call IO_read_dataset_single('r_max', gal_id, info, r_max_kpc)

      ! Sets n1, maximum number of snapshots
      n1 = number_of_redshifts
      ! Sets the initial valid snapshot (uses the disk size as a marker)
      do i=1,n1
        if (galaxy_data(i,1)>=0) then
          init_it = i
          exit
        endif
      enddo
      iread = init_it-1
    endsubroutine read_input_parameters


    subroutine reset_input_params()
      ! Resets the reading of the input parameters file
      iread = 0
      last_output = .false.
    end subroutine reset_input_params


    subroutine set_input_params(gal_id,info)
      ! Reads dimensional input parameters that must be specified and may vary
      ! from galaxy to galaxy and from timestep to timestep

!       character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info, gal_id
      double precision :: next_time_input
      double precision :: current_time_input

      ! Reads the whole file on first access
      if (gal_id /= current_gal_id ) then
        call read_input_parameters(gal_id)
      endif

      iread=iread+1
      
      current_time_input = output_times(iread)

      if ( iread < number_of_redshifts ) then
        next_time_input = output_times(iread+1)
        time_between_inputs = next_time_input-current_time_input
        t = 0 ! At each snapshot, reset the time variable
      else
        last_output = .true.
      endif

      t_Gyr   = next_time_input
      r_disk  = galaxy_data(iread,1)
      v_disk  = galaxy_data(iread,2)
      r_bulge = galaxy_data(iread,3)
      v_bulge = galaxy_data(iread,4)
      r_halo  = galaxy_data(iread,5)
      v_halo  = galaxy_data(iread,6)
      nfw_cs1 = galaxy_data(iread,7)
      Mgas_disk = galaxy_data(iread,8)
      Mstars_disk = galaxy_data(iread,9)
      SFR = galaxy_data(iread,10)
      ! Temporarily setting v_sol_kms to the turbulent speed
      v_sol_kms = p_ISM_sound_speed_km_s * p_ISM_kappa
      l_sol_kpc = p_ISM_turbulent_length

!     DIMENSIONLESS PARAMETERS THAT MUST BE SPECIFIED BUT WILL NOT NORMALLY VARY FROM GALAXY TO GALAXY
!     DIFFUSIVE MAGNETIC HELICITY FLUX (DEFAULTS)
      R_kappa=         1.0d0 !Ratio kappa_t/eta_t of turbulent diffusivities of alpha_m and B

    endsubroutine set_input_params
end module input_params
!*****************************************************
module calc_params
  !Contains parameters that are calculated from the input parameters
  use units
  use math_constants
  use input_params
  use grid
  use input_constants
!
  implicit none
!
  double precision :: etat_sol,etat_sol_kmskpc,etat_sol_cm2s,td_sol,td_sol_kpcskm,td_sol_Gyr,td_sol_s,om0_kmskpc, &
                      Ur_sol_kms,r_sol,lambda,n_sol,r_n,h_sol,r_h,l_sol,r_l,v_sol,r_v,Uz_sol,r_Uz, &
                      Ur_sol,om0,r_om,r1
  contains
    subroutine set_calc_params
!     DIMENSIONAL PARAMETERS THAT CAN BE CALCULATED FROM SPECIFIED DIMENSIONAL PARAMETERS OR THAT ARE NOT NORMALLY VARIED:
!     TURBULENCE
      etat_sol_kmskpc= 1.d0/3*l_sol_kpc*v_sol_kms  !Typical value of etat in units of km kpc/s
      etat_sol_cm2s= 1.d0/3*l_sol_kpc*cm_kpc*v_sol_kms*cm_km  !Typical turbulent diffusivity in units of cm^2/s
      td_sol_kpcskm= h0_kpc**2/etat_sol_kmskpc  !Typical vertical turbulent diffusion timescale in units of kpc s/km
      td_sol_Gyr= h0_kpc**2/etat_sol_cm2s/s_Gyr*cm_kpc*cm_kpc  !Typical vertical turbulent diffusion timescale in units of Gyr
      td_sol_s= td_sol_Gyr*s_Gyr  !Typical vertical turbulent diffusion timescale in units of seconds

!     RADIAL FLOW
      Ur_sol_kms=   0.0d0  !Radial mean velocity at r=r_sol in km/s

!     DIMENSIONLESS PARAMETERS THAT CAN BE CALCULATED OR THAT ARE NOT NORMALLY VARIED:
!     NUMERICAL
      lambda=h0_kpc/r_max_kpc  !Typical aspect ratio of disk

!     TURBULENCE
      h_sol= h_sol_kpc/h0_kpc*h0  !Disk thickness at r=r_sol in units of h0
      r_h= r_h_kpc/r_max_kpc  !Exponential scale radius of disk scale height
      l_sol= l_sol_kpc/h0_kpc*h0  !Size of largest turbulent eddies
      r_l= r_l_kpc/r_max_kpc  !Exponential scale radius of turbulent scale
      v_sol=v_sol_kms/h0_km*t0_s*h0/t0  !Turbulent velocity
      !r_v= r_v_kpc/r_max_kpc  !Exponential scale radius of rms turbulent velocity
      etat_sol=1.d0/3*l_sol*v_sol  !Typical turbulent diffusivity
      td_sol= h_sol**2/etat_sol  !Typical vertical turbulent diffusion timescale

!     RADIAL FLOW
      Ur_sol=Ur_sol_kms/h0_km*t0_s*h0/t0  !Radial mean velocity at r=r_sol

!     SEED FIELD
      r1= r1_kpc/r_max_kpc !Only relevant if Rand_seed=F

!     PHYSICAL GRID
      r_kpc=r*r_max_kpc

  endsubroutine set_calc_params
end module calc_params


