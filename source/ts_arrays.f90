module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by nsteps timesteps)

  use global_input_parameters
  use var
  use grid
  use input_params

  implicit none
  private

  public make_ts_arrays, reset_ts_arrays

  double precision, parameter :: INVALID = -99999d0
  double precision, allocatable, dimension(:),public :: ts_t_Gyr
  character, allocatable, dimension(:),public :: ts_status_code
  double precision, allocatable, dimension(:),public :: ts_Dt
  double precision, allocatable, dimension(:),public :: ts_rmax
  double precision, allocatable, dimension(:),public :: ts_delta_r
  double precision, allocatable, dimension(:,:),public :: ts_Br
  double precision, allocatable, dimension(:,:),public :: ts_Bp
  double precision, allocatable, dimension(:,:),public :: ts_alp_m
  double precision, allocatable, dimension(:,:),public :: ts_Bzmod
  double precision, allocatable, dimension(:,:),public :: ts_h
  double precision, allocatable, dimension(:,:),public :: ts_om
  double precision, allocatable, dimension(:,:),public :: ts_om_b
  double precision, allocatable, dimension(:,:),public :: ts_om_d
  double precision, allocatable, dimension(:,:),public :: ts_om_h
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
  subroutine make_ts_arrays(it,this_t,f,Bzmod,alp,rmax)
    ! Saves the results of simulation (obtained for a particular snapshot it)
    ! to arrays, i.e. stores the time evolution (over snapshots) of the run
    ! N.B. If the grid was extended (i.e. if the galaxy grew in this snapshot)
    ! result will be interpolated into the standard (coarser) grid
    use interpolation
    use profiles
    use messages, only: status_code
    implicit none
    integer, intent(in) :: it
    double precision, intent(in) :: this_t
    double precision, intent(in) :: rmax
    double precision, dimension(:,:), intent(in) :: f
    double precision, dimension(:), intent(in) :: Bzmod
    double precision, dimension(:), intent(in) :: alp

    if (.not.allocated(ts_t_Gyr)) call allocate_ts_arrays()
    if (size(ts_t_Gyr)<it) call reallocate_ts_arrays()

    ts_t_Gyr(it) = this_t
    ts_status_code(it) = status_code

    ts_Dt(it) = t*t0_Gyr
    ts_rmax(it) = rmax
    ts_delta_r(it) = delta_r

    call rescale_array(f(nxghost+1:nx-nxghost,1), ts_Br(it,:))
    call rescale_array(f(nxghost+1:nx-nxghost,2), ts_Bp(it,:))

    if (Dyn_quench) then
      if (.not.Damp) then
        call rescale_array(f(nxghost+1:nx-nxghost,3), ts_alp_m(it,:))
      else
        call rescale_array(f(nxghost+1:nx-nxghost,7), ts_alp_m(it,:))
      endif
    endif
    call rescale_array(h(nxghost+1:nx-nxghost), ts_h(it,:))
    call rescale_array(Om(nxghost+1:nx-nxghost), ts_om(it,:))

    if (p_extra_rotation_curve_outputs) then
      print *, Om_h
      call rescale_array(Om_h(nxghost+1:nx-nxghost), ts_Om_h(it,:))
      call rescale_array(Om_b(nxghost+1:nx-nxghost), ts_Om_b(it,:))
      call rescale_array(Om_d(nxghost+1:nx-nxghost), ts_Om_d(it,:))
    endif

    call rescale_array(G(nxghost+1:nx-nxghost), ts_G(it,:))
    call rescale_array(l(nxghost+1:nx-nxghost), ts_l(it,:))
    call rescale_array(v(nxghost+1:nx-nxghost), ts_v(it,:))
    call rescale_array(etat(nxghost+1:nx-nxghost), ts_etat(it,:))
    call rescale_array(tau(nxghost+1:nx-nxghost), ts_tau(it,:))
    call rescale_array(alp_k(nxghost+1:nx-nxghost), ts_alp_k(it,:))
    call rescale_array(Uz(nxghost+1:nx-nxghost), ts_Uz(it,:))
    call rescale_array(Ur(nxghost+1:nx-nxghost), ts_Ur(it,:))
    call rescale_array(n(nxghost+1:nx-nxghost), ts_n(it,:))
    call rescale_array(Beq(nxghost+1:nx-nxghost), ts_Beq(it,:))
    call rescale_array(r_kpc(nxghost+1:nx-nxghost), ts_rkpc(it,:))
    call rescale_array(Bzmod(nxghost+1:nx-nxghost), ts_Bzmod(it,:))

    ! alp is computed in the gutsdynamo module (annoyingly differently from
    ! anything else). Therefore, one needs to be careful. This is a good
    ! candidate for some code refactoring.
    if (status_code == 'M' .or. status_code == 'm') &
        call rescale_array(alp(nxghost+1:nx-nxghost), ts_alp(it,:))

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
      max_outputs = nsteps+1
    endif

    allocate(ts_t_Gyr(max_outputs))
    ts_t_Gyr = INVALID
    allocate(ts_status_code(max_outputs))

    ! Initializes the time series
    ts_status_code = '-'
    ! Marks all the possible redshifts with 'not run' code
    ! (but only if init_it had been previously initialized)
    if (init_it>0) ts_status_code(init_it:max_it) = '0'

    allocate(ts_Dt(max_outputs))
    ts_Dt = INVALID
    allocate(ts_rmax(max_outputs))
    ts_rmax = INVALID
    allocate(ts_delta_r(max_outputs))
    ts_delta_r = INVALID
    allocate(ts_Br(max_outputs,p_nx_ref))
    ts_Br = INVALID
    allocate(ts_Bp(max_outputs,p_nx_ref))
    ts_Bp = INVALID
    allocate(ts_alp_m(max_outputs,p_nx_ref))
    ts_alp_m = INVALID
    allocate(ts_h(max_outputs,p_nx_ref))
    ts_h = INVALID
    allocate(ts_om(max_outputs,p_nx_ref))
    ts_om = INVALID

    if (p_extra_rotation_curve_outputs) then
      allocate(ts_om_h(max_outputs,p_nx_ref))
      ts_om_h = INVALID
      allocate(ts_om_d(max_outputs,p_nx_ref))
      ts_om_d = INVALID
      allocate(ts_om_b(max_outputs,p_nx_ref))
      ts_om_b = INVALID
    endif

    allocate(ts_G(max_outputs,p_nx_ref))
    ts_G = INVALID
    allocate(ts_l(max_outputs,p_nx_ref))
    ts_l = INVALID
    allocate(ts_v(max_outputs,p_nx_ref))
    ts_v = INVALID
    allocate(ts_etat(max_outputs,p_nx_ref))
    ts_etat = INVALID
    allocate(ts_tau(max_outputs,p_nx_ref))
    ts_tau = INVALID
    allocate(ts_alp_k(max_outputs,p_nx_ref))
    ts_alp_k = INVALID
    allocate(ts_alp(max_outputs,p_nx_ref))
    ts_alp = INVALID
    allocate(ts_Uz(max_outputs,p_nx_ref))
    ts_Uz = INVALID
    allocate(ts_Ur(max_outputs,p_nx_ref))
    ts_Ur = INVALID
    allocate(ts_n(max_outputs,p_nx_ref))
    ts_n = INVALID
    allocate(ts_Beq(max_outputs,p_nx_ref))
    ts_Beq = INVALID
    allocate(ts_Bzmod(max_outputs,p_nx_ref))
    ts_Bzmod = INVALID
    allocate(ts_rkpc(max_outputs,p_nx_ref))
    ts_rkpc = INVALID

  end subroutine allocate_ts_arrays

  subroutine reset_ts_arrays()
    implicit none

    ! Resets the time series arrays
    if (allocated(ts_t_Gyr)) then
      deallocate(ts_t_Gyr)
      deallocate(ts_status_code)
      deallocate(ts_Dt)
      deallocate(ts_rmax)
      deallocate(ts_delta_r)
      deallocate(ts_Br)
      deallocate(ts_Bp)
      deallocate(ts_alp_m)
      deallocate(ts_h)
      deallocate(ts_om)

      if (p_extra_rotation_curve_outputs) then
        deallocate(ts_om_h)
        deallocate(ts_om_d)
        deallocate(ts_om_b)
      endif

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
