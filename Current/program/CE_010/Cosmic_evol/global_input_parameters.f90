! Global (not galaxy specific) input parameters
module global_input_parameters
  ! Contains several switches to control the behaviour of the code

  implicit none

  ! DIRECTORY NAMES
  !Specifies run paths
  character (len=10) :: model_name = 'CE_010'
  character (len=10) :: path_to_input_directories = '../../..'


  integer :: ngals = 10
  integer :: info = 2
  
  ! PARAMETER INPUTS
  !Set to T to read in parameters from file, set to F to use defaults
  logical :: Read_param= .true.

  ! PARAMETER INPUTS
  ! Set to T to read in parameters at several timesteps; set to F to read in parameters at the start only
  logical :: Time_evol= .true.
  integer, parameter :: max_number_of_redshifts = 10 ! lfsr: this will became a variable later... (TODO)

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

  ! OUTFLOW VELOCITY
  ! Set to T to make Uz decrease with radius; set to F to keep Uz constant with radius
  logical :: Var_Uz= .false.

  ! TURBULENT SCALE
  ! Set to T to make l decrease with radius; set to F to keep l constant with radius
  logical :: Var_l= .false.

  ! TURBULENT VELOCITY
  ! Set to T to make v decrease with radius; set to F to keep v constant with radius
  logical :: Var_v= .false.

  ! DENSITY
  ! Set to T to make n decrease with radius; set to F to keep n constant with radius (Beq^2=4pi*n*mp*v^2)
  logical :: Var_n= .true.

  ! SCALE HEIGHT
  ! Set to T to include variable scale height.
  logical :: Flaring= .true.

  ! SEED MAGNETIC FIELD
  ! Set to F for random seed magnetic field, T for some other seed
  logical :: Rand_seed= .false.

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
  logical :: Shear= .true.

  ! ADVECTION
  ! Set to F to turn off advection of the magnetic field
  logical :: Advect= .true.

  ! TURBULENT DIFFUSION
  logical :: Turb_dif= .true.  !Set to F to turn off turbulent diffusion

  namelist /global_pars/ &
    path_to_input_directories, model_name, &
    Read_param, Time_evol, &
    info, &
    ngals, &
!     max_number_of_redshifts, & ! lfsr: this will be added later  (TODO)
    Dyn_quench, Alg_quench, &
    lFloor, &
    Damp, &
    Var_Uz, &
    Var_l, Var_v, &
    Var_n, &
    Flaring, &
    Rand_seed, &
    Alp_ceiling, &
    Om_Brandt, &
    Alp_squared, &
    Krause, &
    Shear, &
    Advect, &
    Turb_dif

  contains

  subroutine read_global_parameters(global_pars_filename)
    ! Reads the
    implicit none
    character(len=*), intent(in) :: global_pars_filename

    write(*,*) 'Using global parameters file: ',global_pars_filename
    open(8,file=global_pars_filename, status='old')
    read(8,nml=global_pars)
    close(8)
  end subroutine read_global_parameters

end module global_input_parameters
