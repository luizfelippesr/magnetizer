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
  integer, parameter :: n1= max_number_of_redshifts!420   !Number of snapshots
  integer :: nread= 1 !Read in new input parameters every nread snapshots
  integer, parameter :: nscreen= 2  !Print output to screen every nscreen*n2 timesteps
  double precision :: tsnap !Time between successive snapshots
  integer :: n2  !Number of timesteps in between snapshots
  double precision :: dt,t=0.d0,first=0.d0, eps_t=0.5!eps_t=0.005d0

  ! Galaxy parameters
  double precision :: r_disk_kpc,R_kappa
  double precision :: l_sol_kpc, r_l_kpc, v_sol_kms, r_v_kpc
  double precision :: n_sol_cm3, r_n_kpc, Uz_sol_kms, r_Uz_kpc
  double precision :: h_sol_kpc, r_h_kpc, Uphi_sol_kms, r_om_kpc

  ! All galaxy data (private!)
  double precision, dimension(max_number_of_redshifts,13), private :: galaxy_data
  character(len=8), private :: current_gal_id_string = 'xxxxxxxx'


  contains
   subroutine set_ts_params(time_between_inputs)
      use units
      double precision, intent(in) :: time_between_inputs

      tsnap = time_between_inputs/t0_Gyr
      n2= max(nint(tsnap/t0_Gyr/eps_t),1)  !Change n2 by changing eps_t
      dt= tsnap/n2  !Timestep in units of t0=h0^2/etat0
      print*,'n2,dt,iread=',n2,dt,iread
      print*,'Uz_sol_kms=',Uz_sol_kms
      print*,''
    endsubroutine set_ts_params


    subroutine read_input_parameters(gal_id_string)
      ! Reads the input parameters file to RAM
      use iso_fortran_env
      character (len=8), intent(in) :: gal_id_string
      integer, parameter :: u_dep = 30
      integer, parameter :: u_indep = 29
      integer i
      integer             :: stat

      current_gal_id_string = gal_id_string

      ! Reads time-dependent parameter values
      open(u_indep, file = trim(path_to_input_directories) // '/input/'  &
                            // trim(model_name) // '/time_indep_params_' &
                            // gal_id_string // '.in', status="old")
      read(u_indep,*) header_time_indep  !Header gives time-independent parameters in order
      read(u_indep,*) r_disk_kpc  !Read time-independent parameter values
      close(u_indep)  !Close file containing time-independent parameter values

      ! Reads time-dependent parameter values
      open(u_dep, file= trim(path_to_input_directories) // '/input/'  &
                      //  trim(model_name) // '/time_dep_params_'     &
                      // gal_id_string // '.in', status="old")
      read(u_dep,*) header_time_dep ! Discards header
      do i=1, max_number_of_redshifts
        read(u_dep, *, iostat=stat) galaxy_data(i,:)
        if (stat == iostat_end) exit
      enddo
      close(u_dep)

      ! Flags abscence of any other values with negatives
      galaxy_data(i:max_number_of_redshifts,:) = -1

    endsubroutine read_input_parameters


    subroutine reset_input_params()
      ! Resets the reading of the input parameters file
      iread = 0
      last_output = .false.
    end subroutine reset_input_params


    subroutine set_input_params(gal_id_string,info)
      ! Reads dimensional input parameters that must be specified and may vary
      ! from galaxy to galaxy and from timestep to timestep

      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
      real :: next_time_input

      if (Read_param) then
        ! Reads the whole file on first access
        if (gal_id_string /= current_gal_id_string) then
          call read_input_parameters(gal_id_string)
        endif

        iread=iread+1

        next_time_input = galaxy_data(iread+1,1)
        if ( next_time_input > 0 ) then
          call set_ts_params(galaxy_data(iread+1,1)-galaxy_data(iread,1))
        else
          last_output = .true.
        endif

        l_sol_kpc    = galaxy_data(iread,2)
        r_l_kpc      = galaxy_data(iread,3)
        v_sol_kms    = galaxy_data(iread,4)
        r_v_kpc      = galaxy_data(iread,5)
        n_sol_cm3    = galaxy_data(iread,6)
        r_n_kpc      = galaxy_data(iread,7)
        Uz_sol_kms   = galaxy_data(iread,8)
        r_Uz_kpc     = galaxy_data(iread,9)
        h_sol_kpc    = galaxy_data(iread,10)
        r_h_kpc      = galaxy_data(iread,11)
        Uphi_sol_kms = galaxy_data(iread,12)
        r_om_kpc     = galaxy_data(iread,13)

        if (info> 0) then
          print *,''
          print *,'Reading new parameters, iread=',iread, ', gal_id_string=',gal_id_string
        endif

      else
        call set_ts_params(0.025d0)
!       NUMERICAL (DEFAULTS)
        r_disk_kpc=     15.0d0  !The maximum radius of the domain in kpc (unit of length along r)
!
!       TURBULENCE (DEFAULTS)
        l_sol_kpc=       0.1d0  !Typical size of largest turbulent eddies, in kpc
        r_l_kpc=       100.0d0!7.5d0  !Exponential scale radius of l in kpc; Relevant only if Var_l=T
        v_sol_kms=      10.0d0  !Typical turbulent velocity in km/s; together with h0_kpc it determ the unit of time
        r_v_kpc=       100.0d0!7.5d0  !Exponential scale radius of l in kpc; Relevant only if Var_v=T
!
!       EQUIPARTITION FIELD (DEFAULTS)
        n_sol_cm3=       1.0d0  !Number density at r=r_sol in cm^-3
        r_n_kpc=         7.5d0  !Exponential scale radius of density in kpc n=n_sol*exp((r-r_sol)/r_n), note that r_Beq=2*r_n if r_v=infty
!
!       VERTICAL WIND (DEFAULTS)
!mark1
        Uz_sol_kms=      20.0d0  !Vertical mean velocity at r=r_sol in km/s (if Var_Uz=F then Uz_kms=Uz_sol_kms)
        r_Uz_kpc=      100.0d0!15.0d0 !Uz_kms=Uz_sol_kms*exp((r_kpc-r_sol_kpc)/r_Uz_kpc) if Var_Uz=T
!
!       SCALE HEIGHT (DEFAULTS)
        h_sol_kpc=       0.5d0 !amplitude of h_kpc=h_sol_kpc*exp((r_kpc-r_sol_kpc)/r_h_kpc)
        r_h_kpc=         9.8d0 !Relevant only if Flaring=1; charac radius of the hyperb/exp varying scale height in kpc
!
!       ROTATION CURVE (DEFAULTS)
        Uphi_sol_kms=  220.0d0  !Circular rotation speed at r=r_sol in km/s
        r_om_kpc=        2.0d0  !Relevant only if Om_Brandt =1; r_om is the charac radius of Brandt profile in kpc
      endif
!
!     DIMENSIONLESS PARAMETERS THAT MUST BE SPECIFIED BUT WILL NOT NORMALLY VARY FROM GALAXY TO GALAXY
!     DIFFUSIVE MAGNETIC HELICITY FLUX (DEFAULTS)
      R_kappa=         1.0d0 !Ratio kappa_t/eta_t of turbulent diffusivities of alpha_m and B
!
!      print*,'r_disk_kpc= ',r_disk_kpc ,'R_kappa=   ',R_kappa
!      print*,'l_sol_kpc=  ',l_sol_kpc  ,'r_l_kpc=   ',r_l_kpc   ,'v_sol_kms=    ',v_sol_kms    ,'r_v_kpc=    ',r_v_kpc
!      print*,'n_sol_cm3=  ',n_sol_cm3  ,'r_n_kpc=   ',r_n_kpc   ,'Uz_sol_kms=   ',Uz_sol_kms   ,'r_Uz_kpc=   ',r_Uz_kpc
!      print*,'h_sol_kpc=  ',h_sol_kpc  ,'r_h_kpc=   ',r_h_kpc   ,'Uphi_sol_kms= ',Uphi_sol_kms ,'r_om_kpc=   ',r_om_kpc
!      print*,''
    endsubroutine set_input_params
end module input_params
