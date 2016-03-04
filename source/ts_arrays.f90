module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by nsteps timesteps)

  use global_input_parameters
  use var
  use grid
  use input_params

  implicit none
  private

  public make_ts_arrays, reset_ts_arrays

  double precision, allocatable, dimension(:),public :: ts_t_Gyr
  double precision, allocatable, dimension(:),public :: ts_Dt
  double precision, allocatable, dimension(:),public :: ts_rmax
  double precision, allocatable, dimension(:),public :: ts_delta_r
  double precision, allocatable, dimension(:,:),public :: ts_Br
  double precision, allocatable, dimension(:,:),public :: ts_Bp
  double precision, allocatable, dimension(:,:),public :: ts_alp_m
  double precision, allocatable, dimension(:,:),public :: ts_Bzmod
  double precision, allocatable, dimension(:,:),public :: ts_h
  double precision, allocatable, dimension(:,:),public :: ts_om
  double precision, allocatable, dimension(:,:),public :: ts_G
  double precision, allocatable, dimension(:,:),public :: ts_l
  double precision, allocatable, dimension(:,:),public :: ts_v
  double precision, allocatable, dimension(:,:),public :: ts_etat
  double precision, allocatable, dimension(:,:),public :: ts_tau
  double precision, allocatable, dimension(:,:),public :: ts_alp_k
  double precision, allocatable, dimension(:,:),public :: ts_alp
  double precision, allocatable, dimension(:,:),public :: ts_Uz
  double precision, allocatable, dimension(:,:),public :: ts_Ur
  double precision, allocatable, dimension(:,:),public :: ts_n
  double precision, allocatable, dimension(:,:),public :: ts_Beq
  double precision, allocatable, dimension(:,:),public :: ts_rkpc

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
    
    if (.not.allocated(ts_t_Gyr)) call allocate_ts_arrays()

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
    ts_rkpc(it,:) = r_kpc(:)
  end subroutine make_ts_arrays
  
  subroutine allocate_ts_arrays()
    ! Initial allocation
    double precision, parameter :: INVALID = -99999d0

      allocate(ts_t_Gyr(n1))
      allocate(ts_Dt(n1))
      allocate(ts_rmax(n1))
      allocate(ts_delta_r(n1))
      allocate(ts_Br(n1,nx))
      allocate(ts_Bp(n1,nx))
      allocate(ts_alp_m(n1,nx))
      allocate(ts_h(n1,nx))
      allocate(ts_om(n1,nx))
      allocate(ts_G(n1,nx))
      allocate(ts_l(n1,nx))
      allocate(ts_v(n1,nx))
      allocate(ts_etat(n1,nx))
      allocate(ts_tau(n1,nx))
      allocate(ts_alp_k(n1,nx))
      allocate(ts_alp(n1,nx))
      allocate(ts_Uz(n1,nx))
      allocate(ts_Ur(n1,nx))
      allocate(ts_n(n1,nx))
      allocate(ts_Beq(n1,nx))
      allocate(ts_Bzmod(n1,nx))
      allocate(ts_rkpc(n1,nx))

      ts_t_Gyr = INVALID
      ts_Dt = INVALID
      ts_rmax = INVALID
      ts_delta_r = INVALID
      ts_Br = INVALID
      ts_Bp = INVALID
      ts_alp_m = INVALID
      ts_Bzmod = INVALID
      ts_h = INVALID
      ts_om = INVALID
      ts_G = INVALID
      ts_l = INVALID
      ts_v = INVALID
      ts_etat = INVALID
      ts_tau = INVALID
      ts_alp_k = INVALID
      ts_alp = INVALID
      ts_Uz = INVALID
      ts_Ur = INVALID
      ts_n = INVALID
      ts_Beq = INVALID
      ts_rkpc = INVALID

    end subroutine allocate_ts_arrays

  subroutine reset_ts_arrays()
    ! Resets the time series arrays
    deallocate(ts_t_Gyr)
    deallocate(ts_Dt)
    deallocate(ts_rmax)
    deallocate(ts_delta_r)
    deallocate(ts_Br)
    deallocate(ts_Bp)
    deallocate(ts_alp_m)
    deallocate(ts_h)
    deallocate(ts_om)
    deallocate(ts_G)
    deallocate(ts_l)
    deallocate(ts_v)
    deallocate(ts_etat)
    deallocate(ts_tau)
    deallocate(ts_alp_k)
    deallocate(ts_alp)
    deallocate(ts_Uz)
    deallocate(ts_Ur)
    deallocate(ts_n)
    deallocate(ts_Beq)
    deallocate(ts_Bzmod)
    deallocate(ts_rkpc)
  end subroutine reset_ts_arrays
end module ts_arrays
