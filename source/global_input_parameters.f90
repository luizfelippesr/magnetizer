! Global (not galaxy specific) input parameters
module global_input_parameters
  ! Contains several switches to control the behaviour of the code

  implicit none

  ! DIRECTORY NAMES
  !Specifies run paths
  character (len=10) :: model_name = 'CE_010'
  character (len=10) :: path_to_input_directories = '../../..'
  character (len=22) :: output_file_name = 'magnetic_galaxies.hdf5'

  integer :: ngals = 10
  integer :: info = 2

  ! PARAMETER INPUTS
  ! Number of timesteps used in the calculation between 2 snapshots
  integer :: nsteps_0 = 30

  ! PARAMETER INPUTS
  ! Set to T to read in parameters at several timesteps; set to F to read in parameters at the start only
  logical :: Time_evol= .true.

  ! ALPHA QUENCHING
  ! Works with alg_quench=F; Set to T for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical :: Dyn_quench= .true.
  ! Works with dyn_quench=F; Set to T for algebraic alpha quenching; F for no quenching
  logical :: Alg_quench= .false.

  ! RANDOM MAGNETIC FIELD
  ! Strength of the rms random magnetic field brms in units of Beq
  logical :: lFloor= .true.

  ! CLOSURE APPROXIMATION
  ! Set to T for FOSA, F for minimal tau approximation
  logical :: Damp= .false.

  ! SEED MAGNETIC FIELD
  ! Set to F for random seed magnetic field, T for some other seed
  logical :: Rand_seed= .false.

  !Seed field amplitude as a fraction of equipartition magnetic field strength
  double precision :: frac_seed= 0.01d0

  ! CEILING ON ALPHA EFFECT
  ! Set to T to put a ceiling for alpha at alpceil*v
  logical :: Alp_ceiling= .true.

  ! BRANDT ROTATION CURVE
  ! Set to T to use a Brandt rotation curve
  logical :: Om_Brandt = .true.

  ! ALPHA^2 EFFECT
  ! Set to T to include alpha^2 effect; set to F to use alpha-omega approximation equations
  logical :: Alp_squared= .true.

  ! KRAUSE'S LAW
  ! Set to T for alpha effect to decrease in proportion to omega (Krause's formula)
  logical :: Krause= .true.

  ! OMEGA EFFECT
  ! Set to T to include Omega effect in dynamo, 0 for alpha^2 dynamo
!   logical :: Shear= .true.

  ! ADVECTION
  ! Set to F to turn off advection of the magnetic field
  logical :: Advect= .true.

  ! TURBULENT DIFFUSION
  logical :: Turb_dif= .true.  !Set to F to turn off turbulent diffusion

  ! Sound speed (in km/s)
  double precision :: p_ISM_sound_speed_km_s = 10d0
  ! Ratio between turbulent velocity and sound speed
  double precision :: p_ISM_kappa = 1d0
  ! Ratio between turbulent pressure and turbulent magnetic field pressure
  double precision :: p_ISM_csi = 1d0
  ! Adiabatic index of the ISM
  double precision :: p_ISM_gamma = 5d0/3d0
  ! Turbulent length (in kpc)
  double precision :: p_ISM_turbulent_length = 0.1
  ! Uses the scale height as an upper limit for the turbulent scale
  logical :: p_limit_turbulent_scale = .true.

  ! Ratio between stellar scale-height and stellar scale-radius
  double precision :: p_stellarHeightToRadiusScale = 1d0/7.3

  ! Ratio between molecular scale-height and molecular scale-radius
  double precision :: p_molecularHeightToRadiusScale = 0.032

  ! Ratio between total gas scale length and stellar scale length
  double precision :: p_gasScaleRadiusToStellarScaleRadius_ratio = 2d0

  ! Molecular fraction calculation
  ! Blitz&Rosolowsky alpha
  double precision :: p_Rmol_alpha = 0.92d0
  ! Blitz&Rosolowsky P0 (in erg/cm^3)
  double precision :: p_Rmol_P0 = 4.787d-12

  ! Outflow calculation ('no_outflow'/'vturb'/'superbubble_simple'/'superbubble')
  character(len=21) :: p_outflow_type = 'superbubble'
  ! Mechanical luminosity associated with the superbubble (in erg/s)
  double precision :: p_outflow_Lsn = 1d38
  ! Ratio between OB associations (superbubbles) rate and supernovae rate
  double precision :: p_outflow_fOB = 0.7
  ! Ratio between supernovae rate and SFR (in \msun^{-1})
  double precision :: p_outflow_etaSN = 9.4d-3
  ! Lifetime of an OB association (in Myr)
  double precision :: p_tOB = 3
  ! Number of supernovae in 1 OB association
  double precision :: p_N_SN1OB = 40
  ! Density of the hot gas (in g/cm^3)
  double precision :: p_outflow_hot_gas_density = 1.7d-27

  ! Check whether the hydrostatic equilibrium solution is correct
  logical :: p_check_hydro_solution = .false.

  ! Runs without solving the dynamo equations
  logical :: p_no_magnetic_fields_test_run = .false.

  ! When set to false, this substitutes any positive shear by 0
  ! (NB positive shears can only arise from the regularisation)
  logical :: p_allow_positive_shears = .false.

  ! Debug mode: all timesteps are included in the output, but
  ! only 1 snapshot is used.
  logical :: p_oneSnaphotDebugMode = .false.

  double precision :: p_turbulent_to_scaleheight_ratio = 0.25
  logical :: p_use_fixed_turbulent_to_scaleheight_ratio = .false.

  namelist /global_pars/ &
    path_to_input_directories, model_name, output_file_name, &
    nsteps_0, &
    Time_evol, &
    info, &
    ngals, &
    Dyn_quench, Alg_quench, &
    lFloor, &
    Damp, &
    Rand_seed, &
    Alp_ceiling, &
    Om_Brandt, &
    Alp_squared, &
    Krause, &
    Advect, &
    Turb_dif, &
    p_ISM_csi, &
    p_ISM_sound_speed_km_s, &
    p_ISM_gamma, &
    p_ISM_turbulent_length, &
    p_limit_turbulent_scale, &
    p_gasScaleRadiusToStellarScaleRadius_ratio, &
    p_stellarHeightToRadiusScale, &
    p_molecularHeightToRadiusScale, &
    p_Rmol_alpha, &
    p_Rmol_P0, &
    p_outflow_type, &
    p_outflow_Lsn, &
    p_outflow_fOB, &
    p_outflow_etaSN, &
    p_tOB, &
    p_N_SN1OB, &
    p_outflow_hot_gas_density, &
    p_check_hydro_solution, &
    p_no_magnetic_fields_test_run, &
    p_allow_positive_shears, &
    p_oneSnaphotDebugMode, &
    p_turbulent_to_scaleheight_ratio, &
    p_use_fixed_turbulent_to_scaleheight_ratio, &
    frac_seed

  contains

  subroutine read_global_parameters(global_pars_filename)
    ! Reads the global parameters file
    implicit none
    character(len=*), intent(in) :: global_pars_filename

    write(*,*) 'Using global parameters file: ',global_pars_filename
    open(8,file=global_pars_filename, status='old')
    read(8,nml=global_pars)
    close(8)
  end subroutine read_global_parameters

end module global_input_parameters
