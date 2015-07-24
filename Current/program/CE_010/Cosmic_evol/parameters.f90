! Input parameters
! At some point we may substitute this by a functions which reads parameter files...)

module direc_names  !Specifies run paths
  implicit none

  character (len=6), parameter :: s0 = 'test'
  character (len=10), parameter :: s1 = '../../..'

end module direc_names


module modules
  ! Contains switches (to be edited before any particular run)
  implicit none


  ! PARAMETER INPUTS
  !Set to T to read in parameters from file, set to F to use defaults
  logical, parameter :: Read_param= .true.

  ! PARAMETER INPUTS
  ! Set to T to read in parameters at several timesteps; set to F to read in parameters at the start only
  logical, parameter :: Time_evol= .true.
  integer, parameter :: max_number_of_redshifts = 10

  ! ALPHA QUENCHING
  ! Works with alg_quench=F; Set to T for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical, parameter :: Dyn_quench= .true.
  ! Works with dyn_quench=F; Set to T for algebraic alpha quenching; F for no quenching
  logical, parameter :: Alg_quench= .false.

  ! RANDOM MAGNETIC FIELD
  ! Strength of the rms random magnetic field brms in units of Beq
  logical, parameter :: lFloor= .true.

  ! CLOSURE APPROXIMATION
  ! Set to T for FOSA, F for minimal tau approximation
  logical, parameter :: Damp= .false.

  ! OUTFLOW VELOCITY
  ! Set to T to make Uz decrease with radius; set to F to keep Uz constant with radius
  logical, parameter :: Var_Uz= .false.

  ! TURBULENT SCALE
  ! Set to T to make l decrease with radius; set to F to keep l constant with radius
  logical, parameter :: Var_l= .false.

  ! TURBULENT VELOCITY
  ! Set to T to make v decrease with radius; set to F to keep v constant with radius
  logical, parameter :: Var_v= .false.

  ! DENSITY
  ! Set to T to make n decrease with radius; set to F to keep n constant with radius (Beq^2=4pi*n*mp*v^2)
  logical, parameter :: Var_n= .true.

  ! SCALE HEIGHT
  ! Set to T to include variable scale height.
  logical, parameter :: Flaring= .true.

  ! SEED MAGNETIC FIELD
  ! Set to F for random seed magnetic field, T for some other seed
  logical, parameter :: Rand_seed= .false.

  ! CEILING ON ALPHA EFFECT
  ! Set to T to put a ceiling for alpha at alpceil*v
  logical, parameter :: Alp_ceiling= .true.

  ! BRANDT ROTATION CURVE
  ! Set to T to use a Brandt rotation curve
  logical, parameter :: Om_Brandt = .true.

  ! ALPHA^2 EFFECT
  ! Set to T to include alpha^2 effect; set to F to use alpha-omega approximation equations
  logical, parameter :: Alp_squared= .true.

  ! KRAUSE'S LAW
  ! Set to T for alpha effect to decrease in proportion to omega (Krause's formula)
  logical, parameter :: Krause= .true.

  ! OMEGA EFFECT
  ! Set to T to include Omega effect in dynamo, 0 for alpha^2 dynamo
  logical, parameter :: Shear= .true.

  ! ADVECTION
  ! Set to F to turn off advection of the magnetic field
  logical, parameter :: Advect= .true.

  ! TURBULENT DIFFUSION
  logical, parameter :: Turb_dif= .true.  !Set to F to turn off turbulent diffusion

end module modules
