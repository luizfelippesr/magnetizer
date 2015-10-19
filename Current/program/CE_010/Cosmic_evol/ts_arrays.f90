module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by nsteps timesteps)

  use global_input_parameters
  use var
  use grid
  use input_params

  implicit none

  double precision, parameter :: INVALID = -99999d0
  double precision, dimension(n1) :: ts_t_Gyr = INVALID
  double precision, dimension(n1) :: ts_Dt = INVALID
  double precision, dimension(n1) :: ts_rmax = INVALID
  double precision, dimension(n1) :: ts_delta_r = INVALID
  double precision, dimension(n1,nx) :: ts_Br = INVALID
  double precision, dimension(n1,nx) :: ts_Bp = INVALID
  double precision, dimension(n1,nx) :: ts_alp_m = INVALID
  double precision, dimension(n1,nx) :: ts_Bzmod = INVALID
  double precision, dimension(n1,nx) :: ts_h = INVALID
  double precision, dimension(n1,nx) :: ts_om = INVALID
  double precision, dimension(n1,nx) :: ts_G = INVALID
  double precision, dimension(n1,nx) :: ts_l = INVALID
  double precision, dimension(n1,nx) :: ts_v = INVALID
  double precision, dimension(n1,nx) :: ts_etat = INVALID
  double precision, dimension(n1,nx) :: ts_tau = INVALID
  double precision, dimension(n1,nx) :: ts_alp_k = INVALID
  double precision, dimension(n1,nx) :: ts_alp = INVALID
  double precision, dimension(n1,nx) :: ts_Uz = INVALID
  double precision, dimension(n1,nx) :: ts_Ur = INVALID
  double precision, dimension(n1,nx) :: ts_n = INVALID
  double precision, dimension(n1,nx) :: ts_Beq = INVALID

contains
  subroutine make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)
    ! Saves the results of simulation (obtained for a particular snapshot it)
    ! to arrays, i.e. stores the time evolution (over snapshots) of the run
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

    print *, 'updating ts', t, t_Gyr, iread
    
    ! lfsr: previously, here it was it+1 instead of it, why?
    ! lfsr: mabe to store the initial condition?
    ts_t_Gyr(it) = t_Gyr
    ts_Dt(it) = t*t0_Gyr
    ts_rmax(it) = rmax
    ts_delta_r(it) = delta_r
    ts_Br(it,:) = f(:,1)
    ts_Bp(it,:) = f(:,2)
    if (Dyn_quench) then   
      if (.not.Damp) then
        ts_alp_m(it,:) = f(:,3)
      else
        ts_alp_m(it,:) = f(:,7)
      endif
    endif
    ts_h(it,:) = h(:)
    ts_om(it,:) = om(:)
    ts_G(it,:) = G(:)
    ts_l(it,:) = l(:)
    ts_v(it,:) = v(:)
    ts_etat(it,:) = etat(:)
    ts_tau(it,:) = tau(:)
    ts_alp_k(it,:) = alp_k(:)
    ts_alp(it,:) = alp(:)
    ts_Uz(it,:) = Uz(:)
    ts_Ur(it,:) = Ur(:)
    ts_n(it,:) =  n(:)
    ts_Beq(it,:) = Beq(:)
    ts_Bzmod(it,:) = Bzmod(:)
  end subroutine make_ts_arrays
  
  subroutine reset_ts_arrays()
    ! Resets the time series arrays
    ts_t_Gyr(:) = INVALID
    ts_Dt(:) = INVALID
    ts_rmax(:) = INVALID
    ts_delta_r(:) = INVALID
    ts_Br(:,:) = INVALID
    ts_Bp(:,:) = INVALID
    ts_alp_m(:,:) = INVALID
    ts_h(:,:) = INVALID
    ts_om(:,:) = INVALID
    ts_G(:,:) = INVALID
    ts_l(:,:) = INVALID
    ts_v(:,:) = INVALID
    ts_etat(:,:) = INVALID
    ts_tau(:,:) = INVALID
    ts_alp_k(:,:) = INVALID
    ts_alp(:,:) = INVALID
    ts_Uz(:,:) = INVALID
    ts_Ur(:,:) = INVALID
    ts_n(:,:) = INVALID
    ts_Beq(:,:) = INVALID
    ts_Bzmod(:,:) = INVALID
  end subroutine reset_ts_arrays
end module ts_arrays
