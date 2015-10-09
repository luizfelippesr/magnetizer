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
module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by nsteps timesteps)
!
  use global_input_parameters
  use var
  use grid
  use input_params
!
  implicit none
!
  double precision, parameter :: INVALID = -99999d0
  double precision, dimension(n1+1) :: ts_t = INVALID
  double precision, dimension(n1+1) :: ts_rmax = INVALID
  double precision, dimension(n1+1) :: ts_delta_r = INVALID
  double precision, dimension(n1+1,nx) :: ts_Br = INVALID
  double precision, dimension(n1+1,nx) :: ts_Bp = INVALID
  double precision, dimension(n1+1,nx) :: ts_alp_m = INVALID
  double precision, dimension(n1+1,nx) :: ts_Bzmod = INVALID
  double precision, dimension(n1+1,nx) :: ts_h = INVALID
  double precision, dimension(n1+1,nx) :: ts_om = INVALID
  double precision, dimension(n1+1,nx) :: ts_G = INVALID
  double precision, dimension(n1+1,nx) :: ts_l = INVALID
  double precision, dimension(n1+1,nx) :: ts_v = INVALID
  double precision, dimension(n1+1,nx) :: ts_etat = INVALID
  double precision, dimension(n1+1,nx) :: ts_tau = INVALID
  double precision, dimension(n1+1,nx) :: ts_alp_k = INVALID
  double precision, dimension(n1+1,nx) :: ts_alp = INVALID
  double precision, dimension(n1+1,nx) :: ts_Uz = INVALID
  double precision, dimension(n1+1,nx) :: ts_Ur = INVALID
  double precision, dimension(n1+1,nx) :: ts_n = INVALID
  double precision, dimension(n1+1,nx) :: ts_Beq = INVALID
!
  contains
    subroutine make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
!
      integer, intent(in) :: it
      double precision, intent(in) :: t
      double precision, intent(in) :: rmax
      double precision, intent(in) :: delta_r
      double precision, dimension(nx,nvar), intent(in) :: f
      double precision, dimension(nx), intent(in) :: Bzmod
      double precision, dimension(nx), intent(in) :: h
      double precision, dimension(nx), intent(in) :: om
      double precision, dimension(nx), intent(in) :: G
      double precision, dimension(nx), intent(in) :: l
      double precision, dimension(nx), intent(in) :: v
      double precision, dimension(nx), intent(in) :: etat
      double precision, dimension(nx), intent(in) :: tau
      double precision, dimension(nx), intent(in) :: alp_k
      double precision, dimension(nx), intent(in) :: alp
      double precision, dimension(nx), intent(in) :: Uz
      double precision, dimension(nx), intent(in) :: Ur
      double precision, dimension(nx), intent(in) :: n
      double precision, dimension(nx), intent(in) :: Beq

      ! lfsr: previously, here it was it+1 instead of it, why?
      ! lfsr: mabe to store the initial condition?
      ts_t(it)=       t
      ts_rmax(it)=    rmax
      ts_delta_r(it)= delta_r
      ts_Br(it,:)=    f(:,1)
      ts_Bp(it,:)=    f(:,2)
      if (Dyn_quench) then   
        if (.not.Damp) then
          ts_alp_m(it,:)=   f(:,3)
        else
          ts_alp_m(it,:)=   f(:,7)
        endif
      endif
      ts_h(it,:)=     h(:)
      ts_om(it,:)=    om(:)
      ts_G(it,:)=     G(:)
      ts_l(it,:)=     l(:)
      ts_v(it,:)=     v(:)
      ts_etat(it,:)=  etat(:)
      ts_tau(it,:)=   tau(:)
      ts_alp_k(it,:)= alp_k(:)
      ts_alp(it,:)=   alp(:)
      ts_Uz(it,:)=    Uz(:)
      ts_Ur(it,:)=    Ur(:)
      ts_n(it,:)=     n(:)
      ts_Beq(it,:)=   Beq(:)
      ts_Bzmod(it,:)= Bzmod(:)
    end subroutine make_ts_arrays
end module ts_arrays
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
        f(nxghost+1 ,1:2)=  0.d0  !Set Br=Bp=0 at r=0 	Dirichlet BC on Br, Bp
        f(nx-nxghost,1:2)=  0.d0  !Set Br=Bp=0 at r=R 	Dirichlet BC on Br, Bp
        do ix=1,nxghost
          f(ix     ,1:2)= -f(2*(nxghost+1)-ix     ,1:2)  !Antisymmetric about r=0    Dirichlet BC on Br, Bp: Br=Bp=0 at r=0
          f(nx+1-ix,1:2)= -f(nx+1-2*(nxghost+1)+ix,1:2)  !Antisymmetric about r=R    Dirichlet BC on Br, Bp: Br=Bp=0 at r=R
          if (Dyn_quench) then
            f(ix     ,nvar  )=  f(2*(nxghost+1)-ix     ,nvar  )  !Symmetric     about r=0    Neumann   BC on alp_m    : dalp_mdr=0 at r=0
            f(nx+1-ix,nvar  )=  f(nx+1-2*(nxghost+1)+ix,nvar  )  !Symmetric     about r=R    Neumann   BC on alp_m    : dalp_mdr=0 at r=R
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

      integer :: iseed,var
      double precision, dimension(nx) :: Bseed
      double precision, dimension(nx,nvar), intent(inout) :: f

      Bseed=frac_seed*Beq

      f(:,:)=0.d0 !Initialize

      if (Rand_seed) then  
        do var=1,2 !Seed r and phi components of B
          do iseed=1,nx
            f(iseed,var)=Bseed(iseed)*random_normal()
          enddo
        enddo
      else
        ! NB the old 'r_disk' variable is now always 1.0 in code units
!       f(:,1)=-Bseed*r/r_disk*(1.d0-r/r_disk)**nn*dexp(-r/r1)
        f(:,1)=-Bseed*r*(1.d0-r)**nn*dexp(-r/r1)
!       f(:,2)= Bseed*r/r_disk*(1.d0-r/r_disk)**nn*dexp(-r/r1)
        f(:,2)= Bseed*r*(1.d0-r)**nn*dexp(-r/r1)
      endif

      call impose_bc(f)
      
    end subroutine init_seed
end module initial_conditions
!*****************************************************
module deriv  !Finite differencing routine for spatial differentiation (same as PENCIL Code)
  use grid
!
  contains
    function xder(f)
!
      implicit none
!
      integer :: ix
      double precision, dimension(nx), intent(in) :: f
      double precision :: fac
      double precision, dimension(nx) :: xder
!
!     assume uniform mesh
!
      fac=1.d0/(60.*(x(2)-x(1)))
!
      do ix=4,nx-3
        xder(ix)=fac*(+45.*(f(ix+1)-f(ix-1)) & 
                      - 9.*(f(ix+2)-f(ix-2)) &
                      +    (f(ix+3)-f(ix-3)) &
                     )
      enddo
!
      xder(1)=fac*(-147.d0*f(1)+360.d0*f(2)-450.d0*f(3)+400.d0*f(4)-225.d0*f(5)+72.d0*f(6)-10.d0*f(7))
      xder(2)=fac*( -10.d0*f(1)- 77.d0*f(2)+150.d0*f(3)-100.d0*f(4)+ 50.d0*f(5)-15.d0*f(6)+ 2.d0*f(7))
      xder(3)=fac*(   2.d0*f(1)- 24.d0*f(2)- 35.d0*f(3)+ 80.d0*f(4)- 30.d0*f(5)+ 8.d0*f(6)-      f(7))
!
!     outer points
!
      xder(nx  )=fac*(147.d0*f(nx)-360.d0*f(nx-1)+450.d0*f(nx-2)-400.d0*f(nx-3)+225.d0*f(nx-4)-72.d0*f(nx-5) &
                      +10.d0*f(nx-6))
      xder(nx-1)=fac*( 10.d0*f(nx)+ 77.d0*f(nx-1)-150.d0*f(nx-2)+100.d0*f(nx-3)- 50.d0*f(nx-4)+15.d0*f(nx-5) &
                      - 2.d0*f(nx-6))
      xder(nx-2)=fac*( -2.d0*f(nx)+ 24.d0*f(nx-1)+ 35.d0*f(nx-2)- 80.d0*f(nx-3)+ 30.d0*f(nx-4)- 8.d0*f(nx-5) &
                      +      f(nx-6))
    end function xder
!
    function xder2(f)
!
      implicit none
!
      integer :: ix
      double precision, dimension(nx), intent(in) :: f
      double precision :: fac
      double precision, dimension(nx) :: xder2
!
!     assume uniform mesh
!
      fac=1.d0/(180.d0*(x(2)-x(1))**2)
!
      do ix=4,nx-3
        xder2(ix)=fac*(-490.d0*f(ix)             &
                       +270.d0*(f(ix+1)+f(ix-1)) &
                       - 27.d0*(f(ix+2)+f(ix-2)) &
                       +  2.d0*(f(ix+3)+f(ix-3)) &
                      )
      enddo
!
      xder2(1)=fac*(812.d0*f(1)-3132.d0*f(2)+5265.d0*f(3)-5080.d0*f(4)+2970.d0*f(5)-972.d0*f(6)+137.d0*f(7))
      xder2(2)=fac*(137.d0*f(1)- 147.d0*f(2)- 255.d0*f(3)+ 470.d0*f(4)- 285.d0*f(5)+ 93.d0*f(6)- 13.d0*f(7))
      xder2(3)=fac*(-13.d0*f(1)+ 228.d0*f(2)- 420.d0*f(3)+ 200.d0*f(4)+  15.d0*f(5)- 12.d0*f(6)+  2.d0*f(7))
!
!     outer points
!
      xder2(nx  )=fac*( 812.d0*f(nx)-3132.d0*f(nx-1)+5265.d0*f(nx-2)-5080.d0*f(nx-3)+2970.d0*f(nx-4) &
                       -972.d0*f(nx-5) +137.d0*f(nx-6))
      xder2(nx-1)=fac*( 137.d0*f(nx)- 147.d0*f(nx-1)- 255.d0*f(nx-2)+ 470.d0*f(nx-3)- 285.d0*f(nx-4) &
                       + 93.d0*f(nx-5) - 13.d0*f(nx-6))
      xder2(nx-2)=fac*( -13.d0*f(nx)+ 228.d0*f(nx-1)- 420.d0*f(nx-2)+ 200.d0*f(nx-3)+  15.d0*f(nx-4) &
                        - 12.d0*f(nx-5) + 2.d0*f(nx-6))
    end function xder2
end module deriv
!*****************************************************
module galaxy_model
  use input_params
  use calc_params
  use profiles
!
  implicit none
!
  contains
    subroutine construct_galaxy_model(gal_id_string,info)
!
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
!
!     SET UP NUMERICS AND PLOTTING
!
      if (Time_evol) then
        call set_input_params(gal_id_string,info)  !read in model parameters
        call set_calc_params  !set other parameters
      endif
!
      call construct_profiles  !construct profiles of h(r), om(r), Uz(r), Ur(r), l(r), v(r), etat(r), tau(r), n(r), Beq(r), alp_k(r)
    end subroutine construct_galaxy_model
end module galaxy_model
!*****************************************************
module bzcalc  !Calculates |Bz| using Div B=0 in the no-z approximation
  use math_constants
  use calc_params
  use var
  use grid
  use input_params
  use deriv
  use profiles
!
  implicit none
!
  double precision, dimension(nx) :: Bzmod
!
  contains
    subroutine estimate_Bzmod(f)
!
      double precision, dimension(nx,nvar), intent(in) :: f
!
      if (r(nxghost+1)==0.d0) then
        Bzmod= lambda*h*abs(xder(f(:,1)))
      else
        Bzmod= lambda*h*abs(f(:,1)/r +xder(f(:,1)))
      endif
    end subroutine estimate_Bzmod
end module bzcalc
!*****************************************************
module equ  !Contains the partial differential equations to be solved
  use galaxy_model
  use deriv
  use boundary_conditions
  use bzcalc
!  
  implicit none
!
  double precision, dimension(nx) :: Br, Bp, Fr, Fp, Er, Ep, dBrdr, d2Brdr2, dBpdr, d2Bpdr2
  double precision, dimension(nx) :: alp_m, alp, dalp_mdr, d2alp_mdr2, dalpdr, detatdr
  double precision, dimension(nx) :: Bsqtot, DivVishniac, Dyn_gen
  double precision, dimension(nx) :: brms, B_floor
  integer, private :: i
  double precision :: rmax, delta_r, hmax, lmax, Ncells
!
  contains
    subroutine pde(f,dfdt)
!
      double precision, dimension(nx,nvar) :: f, dfdt
      intent(inout) :: f
      intent(out) :: dfdt
!
      call impose_bc(f)
!
!     USE EXPLICIT NAMES
      Br=f(:,1)
      Bp=f(:,2)
!
      if (Damp) then
        Fr=f(:,3)
        Fp=f(:,4)
        if (Dyn_quench) then
          Er=f(:,5)
          Ep=f(:,6)
          alp_m=f(:,7)
        endif
      else
        Fr=0*Br
        Fp=0*Bp
        if (Dyn_quench) then
          alp_m=f(:,3)
        endif
      endif
!
      dBrdr=xder(Br)
      d2Brdr2=xder2(Br)
      dBpdr=xder(Bp)
      d2Bpdr2=xder2(Bp)
!
      if (Dyn_quench) then
        dalp_mdr=xder(alp_m)
        d2alp_mdr2=xder2(alp_m)
      endif
!
!     CALCULATE MAGNETIC ENERGY (WITHOUT THE FACTOR 1/(8PI))
      Bsqtot=Br**2 +Bp**2 +Bzmod**2
!
      if (Alg_quench) then
        alp= alp_k/(1.d0 +Bsqtot/Beq**2)  !Formula for simple alpha quenching
      elseif (Dyn_quench) then
        alp= alp_k +alp_m  !Total alpha is equal to the sum of kinetic and magnetic parts
      else
        alp= alp_k
      endif
!
      detatdr= 1.d0/3*l*dvdr+1.d0/3*v*dldr
!
!     IMPOSE MINIMUM (FLOOR) ON B_PHI DUE TO SMALL-SCALE TURBULENT FLUCTUATING MAGNETIC FIELD
!
!     CALCULATE DYNAMO NUMBER
      Dyn_gen=G*alp*h**3/etat**2
!
      if (lFloor) then
        do i=nxghost+1,nx-nxghost
          if (abs(Bp(i))==maxval(abs(Bp))) then
            rmax=r(i)  !radius at max of Bp(r)
            delta_r= 2*dsqrt(abs(Bp(i)/d2Bpdr2(i)))  !width of Gaussian approx to Bp(r)
            hmax=h(i)  !scale height at max of Bp(r)
            lmax=l(i)  !turbulent scale at max of Bp(r)
          endif
        enddo
        Ncells= 3.d0*rmax*delta_r*hmax/lmax**3/lambda**2
        brms= fmag*Beq
        B_floor= brms/dsqrt(Ncells)*lmax/delta_r*lambda/3 !brms/dsqrt(Ncells)*lmax/delta_r*lambda
        B_floor=B_floor*abs(r/rmax)**(1.d0/2)*exp(-(r-rmax)**2/2/(delta_r/2)**2)  !multiply by r^(1/2)*(renormalized Gaussian of width delta_r/2)
!mark2
        alp= alp*(1.d0 +B_floor**2/Bsqtot)  !Formula for simple alpha quenching
      endif
!
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
!
      r(nxghost+1)=0.000001d0
      dalp_mdr(nxghost+1)=0.d0
!
      if (.not.Damp) then
!       CASE 1: FOSA (tau-->0 LIMIT)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
        dfdt(:,1)=      -C_U*Uz*Br/h -lambda*Ur*Br/r -lambda*Ur*dBrdr         &
                   -2.d0/pi/h*ctau*alp*Bp -pi**2/4/h**2*(ctau+Rm_inv)*etat*Br &
                   +(ctau+Rm_inv)*etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2)                     
!                   +ctau*( detatdz*dBrdz -detatdz*lambda*dBzdr)  !Contains detatdz terms
        dfdt(:,2)= G*Br -C_U*Uz*Bp/h -lambda*dUrdr*Bp -lambda*Ur*dBpdr        &
                   -2.d0/pi/h*ctau*alp*Br -pi**2/4/h**2*(ctau+Rm_inv)*etat*Bp &
                   +(ctau+Rm_inv)*etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2) &
                   +ctau*lambda**2*( detatdr*Bp/r +detatdr*dBpdr)  !Contains detatdr terms 
!                   +ctau*(detatdz*dBpdz)  !Contains detatdz terms
        if (.not.Alp_squared) then
          dfdt(:,2)=dfdt(:,2) +2.d0/pi/h*ctau*alp*Br
        endif
        if (Dyn_quench) then
          dfdt(:,3)= -2*(h0_kpc/l_kpc)**2*etat*(ctau*alp*(Br**2+Bp**2+Bzmod**2)/Beq**2                      &
                     +ctau*3*etat/pi**(3.d0/2)/h*abs(Dyn_gen)**(1.d0/2)*Br*Bp/Beq**2+Rm_inv*alp_m) &
                     -C_a*alp_m*Uz/h -lambda*alp_m*Ur/r -lambda*alp_m*dUrdr -lambda*Ur*dalp_mdr    &
                     +R_kappa*etat*(lambda**2*d2alp_mdr2 +lambda**2/r*dalp_mdr +C_d/h**2*alp_m)
        endif
      else
!       CASE 2: MTA (FINITE tau)--NOTE: dfdt BLOWS UP AT ORIGIN BUT SET IT TO 0 ANYWAY
        dfdt(:,1)=       -C_U*Uz*Br/h -lambda*Ur*Br/r -lambda*Ur*dBrdr +Fr                       &
                   +Rm_inv*(-pi**2/4/h**2*etat*Br +etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2))
        dfdt(:,2)=  G*Br -C_U*Uz*Bp/h                 -lambda*dUrdr*Bp   -lambda*Ur*dBpdr +Fp    &
                   +Rm_inv*(-pi**2/4/h**2*etat*Bp +etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2))
        dfdt(:,3)=  tau**(-1)*(-2.d0/pi/h*ctau*alp*Bp -pi**2/4/h**2*ctau*etat*Br                 &
                   +ctau*etat*lambda**2*(-Br/r**2 +dBrdr/r +d2Brdr2) -Fr)                          
!                   +ctau*( detatdz*dBrdz -detatdz*lambda*dBzdr)  !Contains detatdz terms
        dfdt(:,4)=  tau**(-1)*(-2.d0/pi/h*ctau*alp*Br -pi**2/4/h**2*ctau*etat*Bp                 &
                   +ctau*etat*lambda**2*(-Bp/r**2 +dBpdr/r +d2Bpdr2) -Fp                         &
                   +ctau*lambda**2*( detatdr*Bp/r +detatdr*dBpdr))  !Contains detatdr terms
!                   +ctau*detatdz*dBpdz  !Contains detatdz terms
        if (.not.Alp_squared) then
          dfdt(:,4)=dfdt(:,4) +tau**(-1)*2./pi/h*ctau*alp*Br
        endif
        if (Dyn_quench) then
          dfdt(:,5)=  tau**(-1)*(ctau*alp*Br -ctau*etat*pi/2/h                                                 *Bp -Er)
          dfdt(:,6)=  tau**(-1)*(ctau*alp*Bp +ctau*etat*pi/2/h*(1. +1.d0/2/pi**(3.d0/2)*abs(Dyn_gen)**(1.d0/2))*Br -Ep)
          dfdt(:,7)= -2*(h0_kpc/l_kpc)**2*etat*((Er*Br +Ep*Bp)/Beq**2 +Rm_inv*alp_m)             &
                     -C_a*alp_m*Uz/h -lambda*alp_m*Ur/r -lambda*alp_m*dUrdr -lambda*Ur*dalp_mdr  &
                     +R_kappa*etat*(lambda**2*d2alp_mdr2 +lambda**2/r*dalp_mdr +C_d/h**2*alp_m)
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
!
contains
  subroutine rk(f,dfdt)
!
  implicit none
!
  double precision :: gam1, gam2, gam3, zet1, zet2
  double precision, dimension(nx,nvar), intent(inout) :: f, dfdt
  double precision, dimension(nx,nvar) :: pdef, ftmp
  double precision :: ttmp
!
!  Runge Kutta 3rd order time advance
!  f = f(exact) - 0.0046*dt**3*d3f/dt3
!
    gam1=8.d0/15 !gam1, gam2 AND gam3 ARE THE COEFFICIENTS OF THE TIMESTEPS AT WHICH dfdt IS CALCULATED
    gam2=5.d0/12
    gam3=3.d0/4
    zet1=-17.d0/60
    zet2=-5.d0/12
!
    call pde(f,dfdt)
    pdef=dfdt !pde FUNCTION CALCULATES VALUE OF TIME DERIVATIVE ACCORDING TO P.D.E.s
    f=f+dt*gam1*pdef !FIRST GET INTERMEDIATE VALUE OF f (AT t=t+dt*gam1) USING dfdt AT t_0
    t=t+dt*gam1 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet1*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1) USING dfdt AT t_0
    ttmp=t+dt*zet1 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1)
!
    call pde(f,dfdt)
    pdef=dfdt !NOW CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1
    f=ftmp+dt*gam2*pdef !USE THIS TO GET ANOTHER INTERMEDIATE VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2)
    t=ttmp+dt*gam2 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet2*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2) USING dfdt AT t=t_0+dt*gam1
    ttmp=t+dt*zet2 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2)
!
    call pde(f,dfdt)
    pdef=dfdt !CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1+dt*zet1+dt*gam2
    f=ftmp+dt*gam3*pdef !USE THIS TO GET THE FINAL VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2+dt*gam3)
    t=ttmp+dt*gam3 !THEN GO TO THAT TIMESTEP
!
    first=first+1. !COUNTS THE NUMBER OF TIMES RUNGA-KUTTA ROUTINE IS EXCUTED
  end subroutine rk
end module timestep
!*****************************************************
module diagnostic  !Writes diagnostic information to screen, file
  use math_constants
  use units
  use input_params
  use calc_params
  use grid
  use global_input_parameters
  use var
  use galaxy_model
  use initial_conditions
  use equ
  use bzcalc
!
  implicit none
!
  contains
    subroutine print_info(info)
      
      integer, intent(in) :: info
       
      if (info>1) then
        print*,''
        print*,'UNITS'
        print*,'t0 (Gyr)=',t0_Gyr,'t0 (kpc s/km)=',t0_kpcskm,'h0 (kpc)=',h0_kpc,'etat0 (cm2s)= ',etat0_cm2s
        print*,''
        print*,'NUMERICS:'
        print*,'nvar=     ',nvar     
        print*,'dt=       ',dt       ,'   n1=       ',n1       ,'   nsteps=     ',nsteps
        print*,'dx=       ',dx       ,'   nxphys=   ',nxphys   ,'   nxghost=',nxghost,'   nx=',nx
        print*,'minval(x)=',minval(x),'   maxval(x)=',maxval(x)
        print*,''
!         print*,'PHYSICAL GRID:'
!         print*,'r_kpc=',r_kpc
        print*,''
        print*,'MODULES:' 
        print*,'Dyn_quench=',Dyn_quench,'Alg_quench=   ',Alg_quench   ,'Damp=       ' ,Damp
        print*,'Var_Uz=     ',Var_Uz     ,'Var_l=     ',Var_l     ,'Var_v=        ',Var_v        ,'Var_n=      ',Var_n
        print*,'Flaring=    ',Flaring    ,'Rand_seed= ',Rand_seed ,'Alp_ceiling=  ',Alp_ceiling  ,'Om_Brandt=  ',Om_Brandt
        print*,'Alp_squared=',Alp_squared,'Krause=    ',Krause    ,'Shear=        ',Shear        ,'Advect=     ',Advect
        print*,'Turb_dif=   ',Turb_dif
        print*,''
        print*,'INPUT PARAMETERS:'
        print*,'r_max_kpc= ',r_max_kpc ,'R_kappa=   ',R_kappa
        print*,'l_sol_kpc=  ',l_sol_kpc  ,'r_l_kpc=   ',r_l_kpc   ,'v_sol_kms=    ',v_sol_kms    ,'r_v_kpc=    ',r_v_kpc
        print*,'n_sol_cm3=  ',n_sol_cm3  ,'r_n_kpc=   ',r_n_kpc   ,'Uz_sol_kms=   ',Uz_sol_kms   ,'r_Uz_kpc=   ',r_Uz_kpc
        print*,'h_sol_kpc=  ',h_sol_kpc  ,'r_h_kpc=   ',r_h_kpc   ,'Uphi_halfmass_kms= ',Uphi_halfmass_kms ,'r_disk=   ',r_disk
        print*,''
        print*,'OTHER FREE PARAMETERS:'
        print*,'r_in= ',r_in ,'r_max_kpc=',r_max_kpc,'r_sol_kpc=',r_sol_kpc,'r1_kpc=',r1_kpc
        print*,'ctau= ',ctau ,'nn=        ',nn        ,'lambda=   ',lambda
        print*,'C_alp=',C_alp,'alpceil=   ',alpceil   ,'Rm_inv=   ',Rm_inv
        print*,''
        print*,'IMPORTANT CALCULATED PARAMETERS'
        print*,'etat_sol_cm2s=',etat_sol_cm2s,'td_sol_Gyr=',td_sol_Gyr,'tau_sol_Gyr=',tau_sol_Gyr
        print*,'etat_sol=     ',etat_sol     ,'td_sol=    ',td_sol    ,'tau_sol     ',tau_sol
        print*,''
        print*,'INPUT/OUTPUT INFORMATION'
        print*,'Directory ending   =',trim(model_name)
        print*,'home directory     =',trim(path_to_input_directories)
      endif
!
!     WRITE INFO TO FILE "diagnostic.out"
      open(20,file= 'diagnostic.out',status="replace")
!
      write(20,*),'UNITS'
      write(20,*),'t0 (Gyr)=',t0_Gyr,'t0 (kpc s/km)=',t0_kpcskm,'h0 (kpc)=',h0_kpc,'etat0 (cm2s)= ',etat0_cm2s
      write(20,*),''
      write(20,*),'NUMERICS:'
      write(20,*),'nvar=     ',nvar     
      write(20,*),'dt=       ',dt       ,'   n1=       ',n1       ,'   nsteps=     ',nsteps
      write(20,*),'dx=       ',dx       ,'   nxphys=   ',nxphys   ,'   nxghost=',nxghost,'   nx=',nx
      write(20,*),'minval(x)=',minval(x),'   maxval(x)=',maxval(x)
      write(20,*),''
      write(20,*),'PHYSICAL GRID:'
      write(20,*),'r_kpc=',r_kpc
      write(20,*),''
      write(20,*),'MODULES:' 
      write(20,*),'Dyn_quench=',Dyn_quench,'Alg_quench=   ',Alg_quench   ,'Damp=       ' ,Damp
      write(20,*),'Var_Uz=     ',Var_Uz     ,'Var_l=     ',Var_l     ,'Var_v=        ',Var_v        ,'Var_n=      ',Var_n
      write(20,*),'Flaring=    ',Flaring    ,'Rand_seed= ',Rand_seed ,'Alp_ceiling=  ',Alp_ceiling  ,'Om_Brandt=  ',Om_Brandt
      write(20,*),'Alp_squared=',Alp_squared,'Krause=    ',Krause    ,'Shear=        ',Shear        ,'Advect=     ',Advect
      write(20,*),'Turb_dif=   ',Turb_dif
      write(20,*),''
      write(20,*),'INPUT PARAMETERS:'
      write(20,*),'r_max_kpc= ',r_max_kpc ,'R_kappa=   ',R_kappa
      write(20,*),'l_sol_kpc=  ',l_sol_kpc  ,'r_l_kpc=   ',r_l_kpc   ,'v_sol_kms=    ',v_sol_kms    ,'r_v_kpc=    ',r_v_kpc
      write(20,*),'n_sol_cm3=  ',n_sol_cm3  ,'r_n_kpc=   ',r_n_kpc   ,'Uz_sol_kms=   ',Uz_sol_kms   ,'r_Uz_kpc=   ',r_Uz_kpc
      write(20,*),'h_sol_kpc=  ',h_sol_kpc  ,'r_h_kpc=   ',r_h_kpc   ,'Uphi_halfmass_kms= ',Uphi_halfmass_kms ,'r_disk=   ',r_disk
      write(20,*),''
      write(20,*),'OTHER FREE PARAMETERS:'
      write(20,*),'r_in= ',r_in ,'r_max_kpc=',r_max_kpc,'r_sol_kpc=',r_sol_kpc,'r1_kpc=',r1_kpc
      write(20,*),'ctau= ',ctau ,'nn=        ',nn        ,'lambda=   ',lambda
      write(20,*),'C_alp=',C_alp,'alpceil=   ',alpceil   ,'Rm_inv=   ',Rm_inv
      write(20,*),''
      write(20,*),'IMPORTANT CALCULATED PARAMETERS'
      write(20,*),'etat_sol_cm2s=',etat_sol_cm2s,'td_sol_Gyr=',td_sol_Gyr,'tau_sol_Gyr=',tau_sol_Gyr
      write(20,*),'etat_sol=     ',etat_sol     ,'td_sol=    ',td_sol    ,'tau_sol     ',tau_sol
      write(20,*),''
      write(20,*),'INPUT/OUTPUT INFORMATION'
      write(20,*),'Directory ending   =',trim(model_name)
      write(20,*),'home directory     =',trim(path_to_input_directories)
      write(20,*),''
      close(20)
    end subroutine print_info
end module diagnostic
!*****************************************************
module initial_data_dump
  use math_constants
  use units
  use input_params
  use calc_params
  use grid
  use global_input_parameters
  use var
  use galaxy_model
  use initial_conditions
  use equ
  use bzcalc
!
  contains
    subroutine write_initial_data(f,gal_id_string,info)
!
      double precision, dimension(nx,nvar), intent(in) :: f
      character (len=8), intent(in) :: gal_id_string
      integer, intent(in) :: info
      integer, parameter :: diag_unit = 20
!
!     WRITE DATA TO FILE "param.out"
!       if (info> 0) then
!         print*,''
!         print*,'Writing parameters to file param',gal_id_string,'.out'
!       endif
!       open(diag_unit,file= 'diagnostic.out',status="old",position="append")
!       write(diag_unit,*)''
!       write(diag_unit,*)'Writing parameters to file param',gal_id_string,'.out'
!       close(diag_unit)
!       open(10,file= trim(path_to_input_directories) // '/output/' // trim(model_name) &
!                     // '/param_' // gal_id_string // '.out',status="replace")
!       write(10,*)t0_Gyr,t0_kpcskm,h0_kpc,etat0_cm2s,n0_cm3,B0_mkG
!       write(10,*)nvar,dt,n1,nsteps,dx,nxphys,nxghost,nx
!       write(10,*)l_sol_kpc,r_l_kpc,v_sol_kms,r_v_kpc,n_sol_cm3,r_n_kpc,Uz_sol_kms,r_Uz_kpc, &
!                  h_sol_kpc,r_h_kpc,Uphi_halfmass_kms,r_disk,R_kappa
!       write(10,*)r_in,r_max_kpc,r_sol_kpc,r1_kpc,ctau,nn,lambda,C_alp,alpceil,Rm_inv
!       write(10,*)etat_sol_cm2s,td_sol_Gyr,tau_sol_Gyr,etat_sol,td_sol,tau_sol
!       close(10)
!       if (info> 0) then
!         print*,'Finished writing parameters to file param',gal_id_string,'.out'
!         print*,'Writing initial output to file init',gal_id_string,'.out'
!       endif
!       open(diag_unit,file= 'diagnostic.out',status="old",position="append")
!       write(diag_unit,*)'Finished writing parameters to file param',gal_id_string,'.out'
!       write(diag_unit,*)'Writing initial output to file init',gal_id_string,'.out'
!       write(diag_unit,*)''
!       close(diag_unit)
!
!     WRITE DATA TO FILE "init.out"
!       open(11,file= trim(path_to_input_directories) // '/output/' // trim(model_name) &
!                     // '/init_' // gal_id_string // '.out' ,status="replace")
!       write(11,*)r
!       write(11,*)h
!       write(11,*)om
!       write(11,*)G
!       write(11,*)Uz
!       write(11,*)Ur
!       write(11,*)l
!       write(11,*)v
!       write(11,*)etat
!       write(11,*)tau
!       write(11,*)alp_k
!       write(11,*)n
!       write(11,*)Beq
!       write(11,*)f(:,1)
!       write(11,*)f(:,2)
!       write(11,*)Bzmod
!       write(11,*)alp
!       close(11)
!       if (info> 0) then
!         print*,'Finished writing initial output to file init',gal_id_string,'.out'
!       endif
!       open(diag_unit,file= 'diagnostic.out',status="old",position="append")
!       write(diag_unit,*)'Finished writing initial output to file init',gal_id_string,'.out'
!       close(diag_unit)
    end subroutine write_initial_data
end module initial_data_dump
!*****************************************************
