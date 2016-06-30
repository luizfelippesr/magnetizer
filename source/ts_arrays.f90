module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by nsteps timesteps)

  use global_input_parameters
  use var
  use grid
  use input_params
  use messages

  implicit none
  private

  public make_ts_arrays, reset_ts_arrays

  double precision, parameter :: INVALID = -99999d0
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
  subroutine make_ts_arrays(it,this_t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r, invalid_run)
    ! Saves the results of simulation (obtained for a particular snapshot it)
    ! to arrays, i.e. stores the time evolution (over snapshots) of the run
    ! N.B. If the grid was extended (i.e. if the galaxy grew in this snapshot)
    ! result will be interpolated into the standard (coarser) grid
    use interpolation
    implicit none
    integer, intent(in) :: it
    double precision, intent(in) :: this_t
    double precision, intent(in) :: rmax
    double precision, intent(in) :: delta_r
    double precision, dimension(:,:), intent(in) :: f
    double precision, dimension(:), intent(in) :: Bzmod
    double precision, dimension(:), intent(in) :: h
    double precision, dimension(:), intent(in) :: om
    double precision, dimension(:), intent(in) :: G
    double precision, dimension(:), intent(in) :: l
    double precision, dimension(:), intent(in) :: v
    double precision, dimension(:), intent(in) :: etat
    double precision, dimension(:), intent(in) :: tau
    double precision, dimension(:), intent(in) :: alp_k
    double precision, dimension(:), intent(in) :: alp
    double precision, dimension(:), intent(in) :: Uz
    double precision, dimension(:), intent(in) :: Ur
    double precision, dimension(:), intent(in) :: n
    double precision, dimension(:), intent(in) :: Beq
    logical, optional, intent(in) :: invalid_run
    
    if (.not.allocated(ts_t_Gyr)) call allocate_ts_arrays()
    if (size(ts_t_Gyr)<it) call reallocate_ts_arrays()

    ts_t_Gyr(it) = this_t

    if (present(invalid_run)) then
      call message('INVALID RUN', info=2)
      stop
      if (invalid_run) return
    endif

    ts_Dt(it) = t*t0_Gyr
    ts_rmax(it) = rmax
    ts_delta_r(it) = delta_r

    ts_Br(it,:) = rescale_array(f(nxghost+1:nx-nxghost,1), p_nx_ref)
    ts_Bp(it,:) = rescale_array(f(nxghost+1:nx-nxghost,2), p_nx_ref)

    if (Dyn_quench) then
      if (.not.Damp) then
        ts_alp_m(it,:) = rescale_array(f(nxghost+1:nx-nxghost,3), p_nx_ref)
      else
        ts_alp_m(it,:) = rescale_array(f(nxghost+1:nx-nxghost,7), p_nx_ref)
      endif
    endif
    ts_h(it,:) = rescale_array(h(nxghost+1:nx-nxghost), p_nx_ref)
    ts_om(it,:) = rescale_array(om(nxghost+1:nx-nxghost), p_nx_ref)
    ts_G(it,:) = rescale_array(G(nxghost+1:nx-nxghost), p_nx_ref)
    ts_l(it,:) = rescale_array(l(nxghost+1:nx-nxghost), p_nx_ref)
    ts_v(it,:) = rescale_array(v(nxghost+1:nx-nxghost), p_nx_ref)
    ts_etat(it,:) = rescale_array(etat(nxghost+1:nx-nxghost), p_nx_ref)
    ts_tau(it,:) = rescale_array(tau(nxghost+1:nx-nxghost), p_nx_ref)
    ts_alp_k(it,:) = rescale_array(alp_k(nxghost+1:nx-nxghost), p_nx_ref)
    ts_alp(it,:) = rescale_array(alp(nxghost+1:nx-nxghost), p_nx_ref)
    ts_Uz(it,:) = rescale_array(Uz(nxghost+1:nx-nxghost), p_nx_ref)
    ts_Ur(it,:) = rescale_array(Ur(nxghost+1:nx-nxghost), p_nx_ref)
    ts_n(it,:) =  rescale_array(n(nxghost+1:nx-nxghost), p_nx_ref)
    ts_Beq(it,:) = rescale_array(Beq(nxghost+1:nx-nxghost), p_nx_ref)
    ts_Bzmod(it,:) = rescale_array(Bzmod(nxghost+1:nx-nxghost), p_nx_ref)
    ts_rkpc(it,:) = rescale_array(r_kpc(nxghost+1:nx-nxghost), p_nx_ref)
  end subroutine make_ts_arrays
  
  subroutine allocate_ts_arrays()
    ! Initial allocation / initialization
    implicit none
    integer :: max_outputs
    if (.not. p_oneSnaphotDebugMode) then
      ! Default mode: only times corresponding to snapshots are included
      ! in the output.
      max_outputs = n1 ! Gets it from the global variables
    else
      ! In the oneSnaphotDebugMode, all timesteps are included in the output
      ! but only one galform output is used.
      max_outputs = nsteps_0+1
    endif

    allocate(ts_t_Gyr(max_outputs))
    allocate(ts_Dt(max_outputs))
    allocate(ts_rmax(max_outputs))
    allocate(ts_delta_r(max_outputs))
    allocate(ts_Br(max_outputs,nx-2*nxghost))
    allocate(ts_Bp(max_outputs,nx-2*nxghost))
    allocate(ts_alp_m(max_outputs,nx-2*nxghost))
    allocate(ts_h(max_outputs,nx-2*nxghost))
    allocate(ts_om(max_outputs,nx-2*nxghost))
    allocate(ts_G(max_outputs,nx-2*nxghost))
    allocate(ts_l(max_outputs,nx-2*nxghost))
    allocate(ts_v(max_outputs,nx-2*nxghost))
    allocate(ts_etat(max_outputs,nx-2*nxghost))
    allocate(ts_tau(max_outputs,nx-2*nxghost))
    allocate(ts_alp_k(max_outputs,nx-2*nxghost))
    allocate(ts_alp(max_outputs,nx-2*nxghost))
    allocate(ts_Uz(max_outputs,nx-2*nxghost))
    allocate(ts_Ur(max_outputs,nx-2*nxghost))
    allocate(ts_n(max_outputs,nx-2*nxghost))
    allocate(ts_Beq(max_outputs,nx-2*nxghost))
    allocate(ts_Bzmod(max_outputs,nx-2*nxghost))
    allocate(ts_rkpc(max_outputs,nx-2*nxghost))

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
    implicit none

    ! Resets the time series arrays
    if (allocated(ts_t_Gyr)) then
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
    endif
  end subroutine reset_ts_arrays

  subroutine extend_array_vec(ar, length)
    implicit none
    double precision, allocatable, dimension(:,:), intent(inout) :: ar
    double precision, allocatable, dimension(:,:) :: tmp
    integer, dimension(2) :: new_shape, old_shape
    integer, optional, intent(in) :: length
    integer :: how_long

    if (present(length)) then
      how_long = length
    else
      how_long = 10
    endif

    old_shape = shape(ar)
    new_shape = old_shape
    new_shape(1) = new_shape(1)+how_long

    ! Allocates tmp, moves a there
    call move_alloc(ar, tmp)
    allocate(ar(new_shape(1),new_shape(2)))
    ar = INVALID
    ar(:old_shape(1), :) = tmp
    deallocate(tmp)
  end subroutine extend_array_vec

  subroutine extend_array_sca(ar, length)
    implicit none
    double precision, allocatable, dimension(:), intent(inout) :: ar
    double precision, allocatable, dimension(:) :: tmp
    integer :: new_shape, old_shape
    integer, optional, intent(in) :: length
    integer :: how_long

    if (present(length)) then
      how_long = length
    else
      how_long = 1
    endif

    old_shape = size(ar)
    new_shape = old_shape+how_long

    ! Allocates tmp, moves a there
    call move_alloc(ar, tmp)
    allocate(ar(new_shape))
    ar = INVALID
    ar(:old_shape) = tmp
    deallocate(tmp)
  end subroutine extend_array_sca

  subroutine reallocate_ts_arrays()
    implicit none
    ! Moves each ts_array to a temporary place, allocates more space, move the data back.
    ! (There must be a shorter way of writing this!)

    call extend_array_sca(ts_Dt)
    call extend_array_sca(ts_rmax)
    call extend_array_sca(ts_delta_r)
    ! Vectors
    call extend_array_vec(ts_Br)
    call extend_array_vec(ts_Bp)
    call extend_array_vec(ts_alp_m)
    call extend_array_vec(ts_h)
    call extend_array_vec(ts_om)
    call extend_array_vec(ts_G)
    call extend_array_vec(ts_l)
    call extend_array_vec(ts_v)
    call extend_array_vec(ts_etat)
    call extend_array_vec(ts_tau)
    call extend_array_vec(ts_alp_k)
    call extend_array_vec(ts_alp)
    call extend_array_vec(ts_Uz)
    call extend_array_vec(ts_Ur)
    call extend_array_vec(ts_n)
    call extend_array_vec(ts_Beq)
    call extend_array_vec(ts_Bzmod)
    call extend_array_vec(ts_rkpc)
  end subroutine reallocate_ts_arrays
end module ts_arrays
