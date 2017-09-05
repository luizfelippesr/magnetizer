! Global (not galaxy specific) input parameters
! Note to developers: this file should be kept strongly commented, as it will
! also serve as documentation for the parameters
module global_input_parameters
  ! Contains several switches to control the behaviour of the code

  implicit none

  !lfsr: The following should be elsewhere
  integer :: ngals = -1
  integer :: number_of_redshifts = -1

  ! -------------------------------------------------------
  ! Input and output file settings
  ! -------------------------------------------------------
  ! Separate files for input and output
  logical :: p_IO_separate_output = .true.
  ! If p_IO_separate_output==True, the use the following for the outputfile
  character (len=80) :: output_file_name = 'magnetized_galaxies_output.hdf5'
  character (len=80) :: input_file_name  = 'magnetized_galaxies_input.hdf5'
  ! Chunking and compression options
  ! NB currently, March/2016, the HDF5 library does not support filters
  ! (including compression) when using parallel IO. Therefore this options
  ! have limited use.
  logical :: p_IO_chunking = .false.
  integer :: p_IO_number_of_galaxies_in_chunks = 10
  logical :: p_IO_compression = .false. ! requires chunking!
  integer :: p_IO_compression_level = 6

  namelist /io_parameters/ &
    output_file_name, input_file_name, p_IO_separate_output, &
    p_IO_chunking, p_IO_number_of_galaxies_in_chunks, p_IO_compression, &
    p_IO_compression_level


  ! -------------------------------------------------------
  ! Run settings and timestepping parameters
  ! -------------------------------------------------------
  character (len=50) :: model_name = 'Magnetized SAM'
  integer :: info = 1
  integer :: nsteps_0 = 1000
  double precision :: p_courant_v = 0.09 !0.4 for Damp=T and nxphys=151
  double precision :: p_courant_eta = 0.09 !0.3 for Damp=T and nxphys=151
  logical :: p_variable_timesteps = .true.
  integer :: p_nsteps_max = 20000
  integer :: p_MAX_FAILS = 6 ! Number of time-step reductions before giving up
  ! Debug mode: all timesteps are included in the output, but
  ! only 1 snapshot is used.
  logical :: p_oneSnaphotDebugMode = .false.
  ! Runs without solving the dynamo equations
  logical :: p_no_magnetic_fields_test_run = .false.
  ! Seed for the random number generator (can be any integer different from 0)
  integer :: p_random_seed = 17
  ! Number of galaxies after which output file is updated
  integer :: p_ncheckpoint = 5000
  ! The master/root process takes part in the calculation? If .true. it will
  ! distribute work and do some work by himself once every
  ! (nproc+int(nproc/p_master_skip)) galaxies.
  logical :: p_master_works_too = .true.
  double precision :: p_master_skip = 2

  namelist /run_parameters/ &
    model_name, info, nsteps_0, p_courant_v, p_courant_eta, &
    p_variable_timesteps, p_nsteps_max, p_oneSnaphotDebugMode, &
    p_no_magnetic_fields_test_run, p_MAX_FAILS, p_random_seed, &
    p_master_works_too, p_master_skip, p_ncheckpoint

  ! -------------------------------------------------------
  ! Grid settings
  ! -------------------------------------------------------
  ! If rdisk(t_i) < rdisk(t_{i-1}), rescale the B and alpha WITH the disk
  ! (this keeps the number of grid points and is equivalent to a change of
  !  units)
  logical :: p_rescale_field_for_shrinking_disks=.true.
  ! If rdisk(t_i) > rdisk(t_{i-1}), rescale the B and alpha WITH the disk
  ! (this keeps the number of grid points and is equivalent to a change of
  !  units)
  logical :: p_rescale_field_for_expanding_disks=.false.
  logical :: p_scale_back_f_array = .true.
  ! Reference grid size (used both initially and for storage)
  integer :: p_nx_ref=101  !151
  ! Maximum possible grid size, in the case of varying number of grid points
  ! (this is mostly for debugging)
  integer :: p_nx_MAX=10000
  ! Uses a fixed dimensional grid, with r_max_kpc set using the the maximum
  ! half mass radius the galaxy reaches over the entire History
  logical :: p_use_fixed_physical_grid=.false.
  ! The maximum radius to use for computations divided by the half mass radius
  ! i.e. rmax = p_rmax_over_rdisk * rdisk
  double precision :: p_rmax_over_rdisk = 2.75d0
  ! Which method to use for the interpolation, valid options are:
  ! 'simple' (linear, without using fgsl), 'linear' (fgsl),
  ! 'cubic_spline' (fgsl) and 'akima' (fgsl)
  character(len=100) :: p_interp_method = 'simple'
  !!  Warning: the fgsl interpolators have been temporarily disabled

  namelist /grid_parameters/ &
    p_rescale_field_for_shrinking_disks, p_rescale_field_for_expanding_disks, &
    p_scale_back_f_array, p_nx_ref, p_nx_MAX, p_use_fixed_physical_grid, &
    p_rmax_over_rdisk, p_interp_method

  ! -------------------------------------------------------
  ! Dynamo equations parameters and switches
  ! -------------------------------------------------------
  ! ALPHA QUENCHING
  ! Works with alg_quench=F; Set to T for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical :: Dyn_quench = .true.
  ! Works with dyn_quench=F; Set to T for algebraic alpha quenching; F for no quenching
  logical :: Alg_quench = .false.

  ! CLOSURE APPROXIMATION
  ! Set to T for FOSA, F for minimal tau approximation
  logical :: Damp = .false.

  ! SEED MAGNETIC FIELD
  !Seed field amplitude as a fraction of equipartition magnetic field strength
  double precision :: frac_seed = 0.01d0
  ! Selects seed type: 'random', 'decaying' or 'fraction'
  ! fraction -> The seed field is a fixed fraction of the equipartition field
  !             Bseed = Bs = frac_seed*Beq
  ! random   -> The value is drew from a Gaussian distribution with variance Bs
  !             Bseed = Bs*random_gaussian
  ! decaying -> The value decays exponentially with radius
  !             Bseed = r/rmax*(1.d0-r/rmax)**p_nn_seed*dexp(-r/p_r_seed_decay)
  character(len=20) :: p_seed_choice = 'random'
  double precision :: p_r_seed_decay = 15.0d0
  integer:: p_nn_seed = 2 !Only relevant if Rand_seed=F
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
  !Set to F to turn off turbulent diffusion
  logical :: Turb_dif= .true.
  ! Imposes Neumann boundary condition at the maximum radius
  ! I.e. at R=Rmax, symmetric boundary
  logical :: p_neumann_boundary_condition_rmax = .false.
  ! ALPHA EFFECT
  !Factor determining strength of alpha effect
  double precision :: C_alp = 1.0d0
  ! FLOOR
  ! Uses a source term to achieve an effective floor for the magnetic field
  logical :: lFloor = .true.
  !Factor determining strength of the floor
  double precision :: C_floor = 1.0d0
  ! Floor localization (Delta r = p_floor_kappa * p_ISM_turbulent_length)
  ! NB: if using p_use_fixed_turbulent_to_scaleheight_ratio,
  !     this gets inconsistent!
  double precision :: p_floor_kappa = 10d0
  ! Stochastically changes the sign of the floor every Delta r
  logical :: p_space_varying_floor = .true.
  ! Stochastically changes the sign of the floor every
  ! Delta t = p_floor_kappa * tau, where
  !   tau = p_ISM_turbulent_length/(p_ISM_sound_speed_km_s*sqrt(p_ISM_kappa))
  logical :: p_time_varying_floor = .true.
  ! RANDOM MAGNETIC FIELD
  !Strength of the rms random magnetic field brms in units of Beq
  double precision :: fmag = 0.5d0

  namelist /dynamo_parameters/ &
    Dyn_quench, Alg_quench, lFloor, Damp, &
    frac_seed, p_seed_choice, p_r_seed_decay, p_nn_seed, &
    Alp_ceiling, Alp_squared, Krause, Advect, Turb_dif, &
    p_neumann_boundary_condition_rmax, C_alp, C_floor, &
    p_floor_kappa, p_space_varying_floor, p_time_varying_floor, &
    fmag


  ! -------------------------------------------------------
  ! Interstellar medium
  ! -------------------------------------------------------
  ! Sound speed (in km/s)
  double precision :: p_ISM_sound_speed_km_s = 10d0
  ! Ratio between turbulent velocity and sound speed
  double precision :: p_ISM_kappa = 1d0
  ! Ratio between turbulent pressure and turbulent magnetic field pressure
  ! NB if simplified_pressure is on, this correspond to the ratio between
  ! turbulent pressure and _total_ magnetic field pressure
  double precision :: p_ISM_xi =  0.25!1d0 
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
  !lfsr: hum... why did we choose 2 as default?

  ! Molecular fraction calculation
  ! Blitz&Rosolowsky alpha
  double precision :: p_Rmol_alpha = 0.92d0
  ! Blitz&Rosolowsky P0 (in erg/cm^3)
  double precision :: p_Rmol_P0 = 4.787d-12

  ! Check whether the hydrostatic equilibrium solution is correct
  logical :: p_check_hydro_solution = .false.

  ! When set to false, this substitutes any positive shear by 0
  ! (NB positive shears can only arise from the regularisation)
  logical :: p_allow_positive_shears = .false.

  ! If true: uses a fixed l/h ratio, setting it with
  ! l = p_turbulent_to_scaleheight_ratio * h
  logical :: p_use_fixed_turbulent_to_scaleheight_ratio = .false.
  double precision :: p_turbulent_to_scaleheight_ratio = 0.25
  ! Uses a (very) simplified calculation for the mid-plane pressure
  ! where P_B + P_b = \xi P_{turb}
  ! (alternatively, P_B uses the actual B from the dynamo calculation).
  logical :: p_simplified_pressure = .true.
  ! Includes the rotation curve correction in the calculation of the midplane pressure
  logical :: p_enable_P2 = .true.
  ! Assumes constant scaleheight for r<rreg
  logical :: p_truncates_within_rreg = .false.
  ! Sets the regularisation radiso for the rotation curves
  double precision :: p_rreg_to_rdisk = 0.1

  ! Ignores the radial correction to the gravitational potential (P2 in
  ! the paper) and solves the related cubic equation for the scaleheight.
  ! Otherwise, GSL's Brent's method root finder will be used.
  logical :: p_use_legacy_cubic_solver = .false.

  ! Minimum density floor (in g/cm^3)
  double precision :: p_minimum_density = 1d-30 ! g cm^-3

  ! Defines what it means to have a negligible disk
  double precision :: p_rdisk_min=0.5 !kpc
  double precision :: Mgas_disk_min=1d4 ! solar masses
  double precision :: rmin_over_rmax=0.001

  logical :: p_use_Pdm = .true.
  logical :: p_use_Pbulge = .true.
  logical :: p_halo_contraction = .true.
  logical :: p_extra_rotation_curve_outputs = .false.
  logical :: p_extra_pressure_outputs = .false.
  ! Keeps P2>0 to avoid possible artifacts which the regularisation may introduce
  ! Usually not necessary.
  logical :: p_P2_workaround = .false.

  ! Legacy, probably unnecessary, options
  ! Adopts a \rho \propto sech^2(z/h) profile for gas and stars
  ! NB Defaul: rho \propto exp(-|z|/h)
  logical :: p_sech2_profile = .false.
  logical :: p_Pdm_numerical = .false.
  logical :: p_Pbulge_numerical = .false.

  namelist /ISM_and_disk_parameters/ &
    p_ISM_sound_speed_km_s, p_ISM_kappa, p_ISM_xi, p_ISM_gamma, &
    p_ISM_turbulent_length, p_limit_turbulent_scale, &
    p_stellarHeightToRadiusScale, p_molecularHeightToRadiusScale, &
    p_gasScaleRadiusToStellarScaleRadius_ratio, &
    p_Rmol_alpha, p_Rmol_P0, p_check_hydro_solution, &
    p_allow_positive_shears, p_use_fixed_turbulent_to_scaleheight_ratio, &
    p_turbulent_to_scaleheight_ratio, p_simplified_pressure, &
    p_rreg_to_rdisk, p_rdisk_min, Mgas_disk_min, rmin_over_rmax, &
    p_use_legacy_cubic_solver, p_enable_P2, p_truncates_within_rreg, &
    p_minimum_density, p_sech2_profile, p_use_Pdm, p_use_Pbulge, &
    p_halo_contraction, p_extra_rotation_curve_outputs, &
    p_extra_pressure_outputs, p_Pdm_numerical, p_Pbulge_numerical, &
    p_P2_workaround


  ! -------------------------------------------------------
  ! Ouflows
  ! -------------------------------------------------------
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
  double precision :: p_outflow_Vhot = 425 ! eventually, this should be read
                                           ! from the input parameters file
  namelist /outflow_parameters/ &
    p_outflow_type, &
    p_outflow_Lsn, p_outflow_fOB, p_outflow_etaSN, p_tOB, p_N_SN1OB, &
    p_outflow_hot_gas_density, p_outflow_nu0, p_outflow_alphahot, &
    p_outflow_Vhot


  contains

  subroutine read_global_parameters(global_pars_filename)
    ! Reads the multiple namelists in global parameters file
    implicit none
    character(len=*), intent(in) :: global_pars_filename
    integer :: u

!     open(newunit=u,file=global_pars_filename) ! Fortran 2008.. requires new gfortran
    u = 17; open(unit=u,file=global_pars_filename) ! Old Fortran
    ! Reads all the namelists
    ! Note: Rewinding makes the order of the namelists in the file unimportant
    read(u,nml=run_parameters); rewind(u)
    read(u,nml=io_parameters); rewind(u)
    read(u,nml=grid_parameters); rewind(u)
    read(u,nml=dynamo_parameters); rewind(u)
    read(u,nml=outflow_parameters); rewind(u)
    read(u,nml=ISM_and_disk_parameters)
    close(u)
  end subroutine read_global_parameters

end module global_input_parameters
