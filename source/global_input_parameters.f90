! Global (not galaxy specific) input parameters
module global_input_parameters
  ! Contains several switches to control the behaviour of the code

  implicit none

  ! DIRECTORY NAMES
  !Specifies run paths (obsolete)
  character (len=50) :: model_name = 'Magnetized SAM'

  integer :: ngals = -1
  integer :: number_of_redshifts = -1
  integer :: info = 1

  ! IO
  ! Chunking and compression options
  ! NB currently, March/2016, the HDF5 library does not support filters
  ! (including compression) when using parallel IO. Therefore this options
  ! have limited use.
  logical :: p_IO_chunking = .false.
  integer :: p_IO_number_of_galaxies_in_chunks = 10
  logical :: p_IO_compression = .false. ! requires chunking!
  integer :: p_IO_compression_level = 6
  ! Separate files for input and output
  logical :: p_IO_separate_output = .true.
  ! If p_IO_separate_output==True, the use the following for the outputfile
  character (len=80) :: output_file_name = 'magnetized_galaxies_output.hdf5'
  character (len=80) :: input_file_name  = 'magnetized_galaxies_input.hdf5'

  ! PARAMETER INPUTS
  ! Number of timesteps used in the calculation between 2 snapshots
  integer :: nsteps_0 = 180

  ! ALPHA QUENCHING
  ! Works with alg_quench=F; Set to T for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical :: Dyn_quench = .true.
  ! Works with dyn_quench=F; Set to T for algebraic alpha quenching; F for no quenching
  logical :: Alg_quench = .false.

  ! RANDOM MAGNETIC FIELD
  ! Strength of the rms random magnetic field brms in units of Beq
  logical :: lFloor = .true.

  ! CLOSURE APPROXIMATION
  ! Set to T for FOSA, F for minimal tau approximation
  logical :: Damp = .false.

  ! SEED MAGNETIC FIELD
  ! Set to F for random seed magnetic field, T for some other seed
  logical :: Rand_seed = .true.

  !Seed field amplitude as a fraction of equipartition magnetic field strength
  double precision :: frac_seed = 0.01d0

  ! CEILING ON ALPHA EFFECT
  ! Set to T to put a ceiling for alpha at alpceil*v
  logical :: Alp_ceiling = .true.

  ! ALPHA^2 EFFECT
  ! Set to T to include alpha^2 effect; set to F to use alpha-omega approximation equations
  logical :: Alp_squared= .false.

  ! KRAUSE'S LAW
  ! Set to T for alpha effect to decrease in proportion to omega (Krause's formula)
  logical :: Krause= .true.

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
  ! NB if simplified_pressure is on, this correspond to the ratio between
  ! turbulent pressure and _total_ magnetic field pressure
  double precision :: p_ISM_xi = 1d0
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

  ! Outflow calculation ('no_outflow'/'vturb'/'superbubble_simple'/'superbubble/wind')
  character(len=21) :: p_outflow_type = 'no_outflow'
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
  ! Inverse of the depletion timescale of molecular gas (H2+He) in units of Gyr-1
  double precision :: p_outflow_nu0 = 0.5
  ! Exponent in Galform's parametrization of the mass loading
  double precision :: p_outflow_alphahot = -3.2
  ! Velocity scale in Galform's parametrization of the mass loading
  double precision :: p_outflow_Vhot = 425 !

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
  ! If true: uses a fixed l/h ratio, setting it with
  ! l = p_turbulent_to_scaleheight_ratio * h
  logical :: p_use_fixed_turbulent_to_scaleheight_ratio = .false.
  double precision :: p_turbulent_to_scaleheight_ratio = 0.25
  ! Uses a (very) simplified calculation for the mid-plane pressure
  ! where P_B + P_b = \xi P_{turb}
  ! (alternatively, P_B uses the actual B from the dynamo calculation).
  logical :: p_simplified_pressure = .true.

  double precision :: p_rreg_to_rdisk = 0.15

  ! Defines what it means to have a negligible disk
  double precision :: rmin_over_rmax=0.001
  double precision :: Mgas_disk_min=1d4 ! solar masses

  namelist /global_pars/ &
    model_name, &
    output_file_name, input_file_name, &
    nsteps_0, &
    info, &
    Dyn_quench, Alg_quench, &
    lFloor, &
    Damp, &
    Rand_seed, &
    Alp_ceiling, &
    Alp_squared, &
    Krause, &
    Advect, &
    Turb_dif, &
    p_ISM_xi, &
    p_ISM_sound_speed_km_s, &
    p_ISM_gamma, &
    p_ISM_kappa, &
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
    p_outflow_alphahot, &
    p_outflow_Vhot, &
    p_outflow_nu0, &
    p_tOB, &
    p_N_SN1OB, &
    p_outflow_hot_gas_density, &
    p_check_hydro_solution, &
    p_no_magnetic_fields_test_run, &
    p_allow_positive_shears, &
    p_oneSnaphotDebugMode, &
    p_turbulent_to_scaleheight_ratio, &
    p_use_fixed_turbulent_to_scaleheight_ratio, &
    frac_seed, &
    p_simplified_pressure, &
    p_rreg_to_rdisk, &
    p_IO_chunking, &
    p_IO_number_of_galaxies_in_chunks, &
    p_IO_compression, &
    p_IO_compression_level, &
    p_IO_separate_output, &
    Mgas_disk_min, &
    rmin_over_rmax

  contains

  subroutine read_global_parameters(global_pars_filename)
    ! Reads the global parameters file
    implicit none
    character(len=*), intent(in) :: global_pars_filename
    integer, parameter :: u = 17
    open(unit=u,file=global_pars_filename, status='old')
    read(u,nml=global_pars)
    close(u)
  end subroutine read_global_parameters

end module global_input_parameters
