!*****************************************************
! This code contains modules, settings and parameters
! used by the galactic dynamo code dynamo.f90.
!*****************************************************
module var  !Contains subroutine for setting the number of variables for which to solve
  use global_input_parameters
!
  implicit none
  integer :: nvar
!
  contains
    subroutine init_var
!     DEFINE ARRAYS FOR PHYSICAL VARIABLES
      if (.not.Dyn_quench) then
        if (.not.Damp) then
          nvar=2  !vars are Br, Bp
        else
          nvar=4  !vars are Br, Bp, Fr, Fp
        endif
      else
        if (.not.Damp) then
          nvar=3  !vars are Br, Bp, alp_m
        else
          nvar=7  !vars are Br, Bp, Fr, Fp, Er, Ep, alp_m
        endif
      endif
    end subroutine init_var
end module var
!*****************************************************
module boundary_conditions  !Specify and implement boundary conditions at r=0, r=R
  use global_input_parameters
  use grid
  use var
!
  implicit none
!
  contains
    subroutine impose_bc(f)
!
      double precision, dimension(nx,nvar), intent(inout) :: f
      integer :: ix
!
      if (nxghost/=0) then
        !Set Br=Bp=0 at r=0 	Dirichlet BC on Br, Bp
        f(nxghost+1 ,1:2)=  0.d0
        if (.not.p_neumann_boundary_condition_rmax) then
          !Set Br=Bp=0 at r=R 	Dirichlet BC on Br, Bp
          f(nx-nxghost,1:2)=  0.d0
        endif
        do ix=1,nxghost
          !Antisymmetric about r=0    Dirichlet BC on Br, Bp: Br=Bp=0 at r=0
          f(ix     ,1:2)= -f(2*(nxghost+1)-ix     ,1:2)
          if (p_neumann_boundary_condition_rmax) then
          !Symmetric     about r=R    Neumann   BC on alp_m: dBrdr=dBpdr=0 at r=R
            f(nx+1-ix,1:2)=  f(nx+1-2*(nxghost+1)+ix,1:2)
          else
            !Antisymmetric about r=R    Dirichlet BC on Br, Bp: Br=Bp=0 at r=R
            f(nx+1-ix,1:2)= -f(nx+1-2*(nxghost+1)+ix,1:2)
          endif

          if (Dyn_quench) then
            !Symmetric     about r=0    Neumann   BC on alp_m: dalp_mdr=0 at r=0
            f(ix     ,nvar  )=  f(2*(nxghost+1)-ix     ,nvar  )
            !Symmetric     about r=R    Neumann   BC on alp_m: dalp_mdr=0 at r=R
            f(nx+1-ix,nvar  )=  f(nx+1-2*(nxghost+1)+ix,nvar  )
          endif
        enddo
      endif
    end subroutine impose_bc
end module boundary_conditions
!*****************************************************
module initial_conditions  !Set initial conditions
  use calc_params
  use grid
  use var
  use boundary_conditions
  use random
  use profiles
!
  implicit none
!
  contains
    subroutine init_seed(f)
      use global_input_parameters
      integer :: iseed,var
      double precision, dimension(nx) :: Bseed
      double precision, dimension(:,:), intent(inout) :: f
      double precision :: r1
      ! Initializes the seed field to a fraction of the equipartition field
      Bseed=frac_seed*Beq
      ! Initializes f
      f(:,:)=0.d0

      select case (trim(p_seed_choice))
        case('fraction')
          ! The field is a fixed fraction of the equipartition field
          f(:,1)=-Bseed
          f(:,2)= Bseed
        case('decaying')
          ! The field
          r1 = p_r_seed_decay/r_max_kpc
          f(:,1)=-Bseed*r*(1.d0-r)**p_nn_seed*dexp(-r/r1)
          f(:,2)= Bseed*r*(1.d0-r)**p_nn_seed*dexp(-r/r1)
        case('random')
          ! The field at each point (and component) is a Gaussian drawing from
          ! a distribution with variance Bseed
          do var=1,2 !Seed r and phi components of B
            do iseed=1,nx
              f(iseed,var)=Bseed(iseed)*random_normal()
            enddo
          enddo
        case default
          print *, 'init_seed: Invalid option ',p_seed_choice
          stop
      end select
      call impose_bc(f)
    end subroutine init_seed
end module initial_conditions
!*****************************************************
module bzcalc  !Calculates |Bz| using Div B=0 in the no-z approximation
  use math_constants
  use calc_params
  use var
  use grid
  use input_params
  use deriv
  use profiles
  implicit none
  double precision, dimension(:), allocatable :: Bzmod

  contains
    subroutine estimate_Bzmod(f)
      double precision, dimension(:,:), intent(in) :: f
      double precision, dimension(nx) :: dBrdr

      ! Checks if Bzmod has the correct shape
      call check_allocate(Bzmod)
      ! Computes dB_r/dr
      dBrdr = xder(f(:,1))
      ! Estimates Bzmod
      Bzmod = lambda*h*abs(f(:,1)/r + dBrdr)
      ! When r->0, use B_r/r \approx dB_r/dr
      Bzmod(nxghost+1) = lambda*h(nxghost+1)*abs(2d0*dBrdr(nxghost+1))

    end subroutine estimate_Bzmod
end module bzcalc
!*****************************************************
module equ  !Contains the partial differential equations to be solved
  use input_params
  use calc_params
  use profiles
  use deriv
  use boundary_conditions
  use bzcalc

  implicit none
  integer, private :: i
  double precision, dimension(:), allocatable :: alp_m, alp
  double precision :: rmax, hmax, lmax!, Ncells

  contains
    subroutine pde(f,dfdt)
      !     LIST OF VARIABLE NAMES FOR f ARRAY
      !
      !     UNDER FOSA          UNDER TAU APPROXIMATION
      !     f(:,1)=Br           f(:,1)=Br
      !     f(:,2)=Bp           f(:,2)=Bp
      !     f(:,3)=alp_m        f(:,3)=Fr
      !                         f(:,4)=Fp
      !                         f(:,5)=Er
      !                         f(:,6)=Ep
      !                         f(:,7)=alp_m
      !
      double precision, dimension(:,:), target, intent(inout) :: f
      double precision, dimension(:,:), target, intent(out) :: dfdt
      ! Convenience pointers to positions in the f-array
      double precision, dimension(:), pointer :: Br, Bp
      double precision, dimension(:), pointer :: Fr, Fp, Er, Ep
      ! Auxiliary arrays
      double precision, dimension(nx) :: dBrdr, d2Brdr2
      double precision, dimension(nx) :: dBpdr,d2Bpdr2
      double precision, dimension(nx) :: dalp_mdr, d2alp_mdr2
      double precision, dimension(nx) :: detatdr
      double precision, dimension(nx) :: Bsqtot, Dyn_gen
      double precision, dimension(nx) :: Dyn_crit, Dtilde, R_U
      double precision, dimension(nx) :: brms, B_floor
      double precision, dimension(nx) :: Afloor, Ncells
      ! Other
      double precision, dimension(nx), target :: zeros_array

      call impose_bc(f)

      zeros_array = 0.0
      call check_allocate(alp_m)
      call check_allocate(alp)

      ! Sets first convenience pointers and variables
      Br => f(:,1)
      Bp => f(:,2)

      if (Damp) then
        Fr => f(:,3)
        Fp => f(:,4)
        if (Dyn_quench) then
          Er => f(:,5)
          Ep => f(:,6)
          alp_m=f(:,7)
        endif
      else
        Fr => zeros_array
        Fp => zeros_array
        if (Dyn_quench) then
          alp_m=f(:,3)
        endif
      endif

      ! Computes derivatives
      dBrdr = xder(Br)
      d2Brdr2 = xder2(Br)
      dBpdr = xder(Bp)
      d2Bpdr2 = xder2(Bp)

      if (Dyn_quench) then
        dalp_mdr = xder(alp_m)
        d2alp_mdr2 = xder2(alp_m)
      endif

      ! CALCULATE MAGNETIC ENERGY (WITHOUT THE FACTOR 1/(8PI))
      Bsqtot= Br**2 + Bp**2 + Bzmod**2

      if (Alg_quench) then
        alp = alp_k/(1.d0 +Bsqtot/Beq**2)  !Formula for simple alpha quenching
      elseif (Dyn_quench) then
        alp = alp_k +alp_m  !Total alpha is equal to the sum of kinetic and magnetic parts
      else
        alp = alp_k
      endif

      detatdr = 1.d0/3*l*dvdr + 1.d0/3*v*dldr

      ! CALCULATE DYNAMO NUMBER
      Dyn_gen = G*alp*h**3/etat**2

      ! IMPOSE MINIMUM (FLOOR) ON B_PHI DUE TO SMALL-SCALE TURBULENT
      ! FLUCTUATING MAGNETIC FIELD
      if (lFloor) then
        !Estimate of critical dynamo number
        Dyn_crit= -(pi/2)**5*(1.d0 +4*C_U*R_U/pi**2)**2  
        !Ratio of Dynamo number to critical dynamo number
        Dtilde= Dyn_gen/Dyn_crit
        !Dimensionless outflow parameter
        R_U= Uz*h/etat
        !Coefficient of Bfloor for source term to get Bp=Bfloor
        Afloor= (pi/2)**2*etat/h**2*(1.d0+4*C_U*R_U/pi**2)*(1.d0-Dtilde)
!         Afloor=1d0
        !Number of turbulent cells in the annular volume
        Ncells= abs(3.d0*r*Delta_r*h/l**3/lambda**2)
        !Small-scale magnetic field strength
        brms= fmag*Beq
        !Floor magnetic field
        B_floor= exp(-Delta_r/2/r)*brms/Ncells**(1d0/2d0)*l/Delta_r*lambda/3
        !!Old version:
        !do i=nxghost+1,nx-nxghost
        !  if (abs(Bp(i))==maxval(abs(Bp))) then
        !    rmax=r(i)  !radius at max of Bp(r)
        !    delta_r= 2*dsqrt(abs(Bp(i)/d2Bpdr2(i)))  !width of Gaussian approx
        !                                             ! to Bp(r)
        !    hmax=h(i)  !scale height at max of Bp(r)
        !    lmax=l(i)  !turbulent scale at max of Bp(r)
        !  endif
        !enddo
        !Ncells= 3.d0*rmax*delta_r*hmax/lmax**3/lambda**2
        !brms= fmag*Beq
        !B_floor= brms/dsqrt(Ncells)*lmax/delta_r*lambda/3 !brms/dsqrt(Ncells)*lmax/delta_r*lambda
        !B_floor=B_floor*abs(r/rmax)**(1.d0/2)*exp(-(r-rmax)**2/2/(delta_r/2)**2)  !multiply by r^(1/2)*(renormalized Gaussian of width delta_r/2)
!mark2
        !alp= alp*(1.d0 +B_floor**2/Bsqtot)  !Formula for simple alpha quenching
      else
        Afloor= 0d0
        B_floor= 0d0
      endif

!     LIST OF VARIABLE NAMES FOR f ARRAY
!
!     UNDER FOSA          UNDER TAU APPROXIMATION
!     f(:,1)=Br           f(:,1)=Br
!     f(:,2)=Bp           f(:,2)=Bp
!     f(:,3)=alp_m        f(:,3)=Fr
!                         f(:,4)=Fp
!                         f(:,5)=Er
!                         f(:,6)=Ep
!                         f(:,7)=alp_m
!
!     METHOD: ALL ARRAYS ARE 2-DIMENSIONAL

!       r(nxghost+1)=0.000001d0  ! LFSR: this seems to be spurious
      dalp_mdr(nxghost+1)=0.d0
!
      if (.not.Damp) then
        ! CASE 1: FOSA (tau-->0 LIMIT)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
        !            Vertical velocity terms
        dfdt(:,1)= -C_U*Uz*Br/h                                                         &
                   ! alpha effect
                   -2d0/pi/h*ctau*alp*Bp                                                &
                   ! Vertical diffusion
                   -pi**2/4/h**2*(ctau+Rm_inv)*etat*Br                                  &
                   ! Radial diffusion indep of dhdr
                   +(ctau+Rm_inv)*etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2)
                   ! Following commented out for simplicity
                   !-lambda*Ur*Br/r -lambda*Ur*dBrdr

                   !Omega effect
        dfdt(:,2)=  G*Br                                                                &
                   !Vertical velocity terms
                   -C_U*Uz*Bp/h                                                         &
                   !alpha^2 effect
                   -2d0/pi/h*ctau*alp*Br                                                &
                   !Vertical diffusion
                   -pi**2/4/h**2*(ctau+Rm_inv)*etat*Bp                                  &
                   !Radial diffusion propto etat
                   +(ctau+Rm_inv)*etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2)           &
                   !Radial diffusion propto detatdr
                   +ctau*detatdr*lambda**2*( Bp/r +dBpdr)                               &
                   !Forcing function that enforces floor at B_floor
                   +Afloor*B_floor
                   ! Following commented out for simplicity
                   !   -lambda*dUrdr*Bp -lambda*Ur*dBpdr

        if (.not.Alp_squared) then
          ! Remove alpha^2 effect if turned off
          dfdt(:,2)= dfdt(:,2) +2.d0/pi/h*ctau*alp*Br
        endif
        if (Dyn_quench) then
          dfdt(:,3)= -2*(h0_kpc/l_kpc)**2*etat*(                                        &
                     !Emf.B term 1 (alpha)
                     ctau*alp*(Br**2+Bp**2+0d0*Bzmod**2)/Beq**2                         &
                     ! Emf.B term 2 (etat)
                     +ctau*3*etat/pi**(3d0/2d0)/h*abs(Dyn_gen)**(1d0/2d0)*Br*Bp/Beq**2  &
                     ! Ohmic dissipation
                     +Rm_inv*alp_m)                                                     &
                     ! Vertical velocity terms
                     -C_a*alp_m*Uz/h                                                    &
                     ! Vertical diffusion
                     +R_kappa*etat*C_d/h**2*alp_m                                       &
                     ! Radial diffusion propto etat, indep of dhdr
                     +R_kappa*etat*lambda**2*(d2alp_mdr2 +dalp_mdr/r)                   &
                     ! Rad diff propto detatdr, indep of dhdr
                     +R_kappa*detatdr*lambda**2*dalp_mdr
        endif
      else
!       CASE 2: MTA (FINITE tau)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
                   !Curl of emf from equation 3
        dfdt(:,1)=  Fr                                                                  &
                   !Vertical velocity
                   -C_U*Uz*Br/h                                                         &
                   !Ohmic vertical diffusion
                   +Rm_inv*etat*(-pi**2/4/h**2*Br                                       &
                   !Ohmic radial diffusion
                   +lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2))                     
                   ! Following commented out for simplicity
                   !-lambda*Ur*Br/r -lambda*Ur*dBrdr
!
                   !Curl of emf from equation 4
        dfdt(:,2)=  Fp                                                                  &
                   !Omega effect
                   +G*Br                                                                &
                   !Vertical velocity
                   -C_U*Uz*Bp/h                                                         &
                   !Ohmic vertical diffusion
                   +Rm_inv*etat*(-pi**2/4/h**2*Bp                                       &
                   !Ohmic radial diffusion
                   +lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2))                   
                   ! Following commented out for simplicity
                   !-lambda*dUrdr*Bp   -lambda*Ur*dBpdr
!
        dfdt(:,3)=  tau**(-1)*(                                                         &
                   !alpha effect
                   -2d0/pi/h*ctau*alp*Bp                                                &
                   !Turbulent vertical diffusion
                   -pi**2/4/h**2*ctau*etat*Br                                           &
                   !Turbulent radial diffusion
                   +ctau*etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2)                    &
                   !Damping term
                   -Fr)                 
                   !Following commented out for now but may be useful later
                   !+ctau*( detatdz*dBrdz -detatdz*lambda*dBzdr)
!
        dfdt(:,4)=  tau**(-1)*(                                                         &
                   !alpha^2 effect
                   -2d0/pi/h*ctau*alp*Br                                                &
                   !Turbulent vertical diffusion
                   -pi**2/4/h**2*ctau*etat*Bp                                           &
                   !Turbulent radial diffusion propto etat
                   +ctau*etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2)                    &
                   !Turbulent radial diffusion propto detatdr
                   +ctau*detatdr*lambda**2*( Bp/r +dBpdr)                               &
                   !Damping term
                   -Fp)
                   !Following commented out for now but may be useful later
                   !+ctau*detatdz*dBpdz
!
        if (.not.Alp_squared) then
          !Remove alpha^2 effect if turned off
          dfdt(:,4)= dfdt(:,4) +tau**(-1)*2d0/pi/h*ctau*alp*Br           
        endif
        if (Dyn_quench) then
          dfdt(:,5)=  tau**(-1)*(                                                       &
                     !alpha part of E.B (r-component)
                      ctau*alp*Br                                                       &
                     !etat part of E.B (r-component)
                     -ctau*etat*pi/2/h*Bp                                               &
                     !Damping term
                     -Er)
!
          dfdt(:,6)=  tau**(-1)*(                                                       &
                     !alpha part of Emf.B (phi-component)
                      ctau*alp*Bp                                                       &
                     !etat part of Emf.B (phi-component)
                     +ctau*etat*pi/2/h*Br*(1d0 +1d0/2/pi**(3d0/2)*abs(Dyn_gen)**(1d0/2)) & 
                     !Damping term
                     -Ep)
!
          dfdt(:,7)= -2*(h0_kpc/l_kpc)**2*etat*(                                         &
                     !Emf.B
                     (Er*Br +Ep*Bp)/Beq**2                                               & 
                     !Ohmic dissipation
                     +Rm_inv*alp_m)                                                      &
                     !Vertical velocity
                     -C_a*alp_m*Uz/h                                                     &
                     !Turbulent vertical diffusion
                     +R_kappa*etat*(C_d/h**2*alp_m                                       &
                     !Turbulent radial diffusion propto etat
                     +lambda**2*(d2alp_mdr2 +dalp_mdr/r))                                &
                     !Turbulent radial diffusion propto detatdr
                     +R_kappa*detatdr*lambda**2*dalp_mdr
                     ! Following commented out for simplicity
                     !-lambda*alp_m*Ur/r -lambda*alp_m*dUrdr -lambda*Ur*dalp_mdr
!
        endif
      endif
    end subroutine pde

end module equ
!*****************************************************
module timestep  !Contains time-stepping routine
  use var
  use grid
  use input_params
  use equ

contains

  subroutine rk(f)

  implicit none

  double precision :: gam1, gam2, gam3, zet1, zet2
  double precision, dimension(nx,nvar), intent(inout) :: f
  double precision, dimension(nx,nvar) :: dfdt
  double precision, dimension(nx,nvar) :: pdef, ftmp
  double precision :: ttmp

    !  Runge Kutta 3rd order time advance
    !  f = f(exact) - 0.0046*dt**3*d3f/dt3

    gam1=8.d0/15 !gam1, gam2 AND gam3 ARE THE COEFFICIENTS OF THE TIMESTEPS AT WHICH dfdt IS CALCULATED
    gam2=5.d0/12
    gam3=3.d0/4
    zet1=-17.d0/60
    zet2=-5.d0/12

    call pde(f,dfdt)
    pdef=dfdt !pde FUNCTION CALCULATES VALUE OF TIME DERIVATIVE ACCORDING TO P.D.E.s
    f=f+dt*gam1*pdef !FIRST GET INTERMEDIATE VALUE OF f (AT t=t+dt*gam1) USING dfdt AT t_0
    t=t+dt*gam1 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet1*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1) USING dfdt AT t_0
    ttmp=t+dt*zet1 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1)

    call pde(f,dfdt)
    pdef=dfdt !NOW CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1
    f=ftmp+dt*gam2*pdef !USE THIS TO GET ANOTHER INTERMEDIATE VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2)
    t=ttmp+dt*gam2 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet2*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2) USING dfdt AT t=t_0+dt*gam1
    ttmp=t+dt*zet2 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2)

    call pde(f,dfdt)
    pdef=dfdt !CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1+dt*zet1+dt*gam2
    f=ftmp+dt*gam3*pdef !USE THIS TO GET THE FINAL VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2+dt*gam3)
    t=ttmp+dt*gam3 !THEN GO TO THAT TIMESTEP

    first=first+1. !COUNTS THE NUMBER OF TIMES RUNGA-KUTTA ROUTINE IS EXCUTED
  end subroutine rk
end module timestep
