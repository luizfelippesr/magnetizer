!# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
!#
!# This file is part of Magnetizer.
!#
!# Magnetizer is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# Magnetizer is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
!#
module ts_arrays  !Contains subroutine that stores time series data (n1 snapshots, separated by nsteps timesteps)

  use global_input_parameters
  use var
  use grid
  use input_params
  use tsDataObj
  implicit none
  private

  public make_ts_arrays, reset_ts_arrays

  double precision, parameter :: INVALID = -99999d0

  character(len=15), dimension(5),parameter :: scalar_names = [         &
                                                 't_Gyr         ', &
                                                 'dt            ', &
                                                 'rmax          ', &
                                                 'Bmax          ', &
                                                 'Bmax_idx      ']
  character(len=15), dimension(28), parameter :: profile_names = [     &
                                                 'Br           ', &
                                                 'Bp           ', &
                                                 'alp_m        ', &
                                                 'Bzmod        ', &
                                                 'h            ', &
                                                 'Omega        ', &
                                                 'Omega_b      ', &
                                                 'Omega_d      ', &
                                                 'Omega_h      ', &
                                                 'Shear        ', &
                                                 'l            ', &
                                                 'v            ', &
                                                 'etat         ', &
                                                 'tau          ', &
                                                 'alp_k        ', &
                                                 'alp          ', &
                                                 'Uz           ', &
                                                 'Ur           ', &
                                                 'n            ', &
                                                 'Beq          ', &
                                                 'rkpc         ', &
                                                 'P            ', &
                                                 'Pd           ', &
                                                 'Pm           ', &
                                                 'Pstars       ', &
                                                 'Pbulge       ', &
                                                 'Pdm          ', &
                                                 'P2           ']

  character, allocatable, dimension(:),public :: ts_status_code

  type(tsData), public :: ts_data

contains

  subroutine make_ts_arrays(it,this_t,f,Bzmod,alp)
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
    double precision, dimension(:,:), intent(in) :: f
    double precision, dimension(:), intent(in) :: Bzmod
    double precision, dimension(:), intent(in) :: alp
    double precision, dimension(p_nx_ref) :: Btot
    double precision, dimension(p_nx_ref) :: tmp
    integer :: Bmax_idx
    integer :: max_outputs

    ! Initializes the ts_data object if necessary
    if (.not.ts_data%initialized) then
      if (.not. p_oneSnaphotDebugMode) then
        ! Default mode: only times corresponding to snapshots are included
        ! in the output.
        max_outputs = n1 ! Gets it from the global variables
      else
        ! In the oneSnaphotDebugMode, all timesteps are included in the output
        ! but only one galform output is used.
        max_outputs = nsteps+1
      endif

      if (.not.output_everything) then
        ts_data = tsData(scalar_names, profile_names, max_outputs, p_nx_ref, &
                       output_quantities_list)
      else
        ts_data = tsData(scalar_names, profile_names, max_outputs, p_nx_ref)
      endif

      allocate(ts_status_code(max_outputs))
      ! Initializes the time series
      ts_status_code = '-'
      ! Marks all the possible redshifts with 'not run' code
      ! (but only if init_it had been previously initialized)
      if (init_it>0) ts_status_code(init_it:max_it) = '0'
    endif

    ts_status_code(it) = status_code

    call ts_data%set_scalar('dt', it, t*t0_Gyr)

    ! Reads and stores the magnetic field
    call rescale_array(f(nxghost+1:nx-nxghost,1), tmp); call ts_data%set('Br', it, tmp)
    call rescale_array(f(nxghost+1:nx-nxghost,2), tmp); call ts_data%set('Bp', it, tmp)
    call rescale_array(Bzmod(nxghost+1:nx-nxghost), tmp); call ts_data%set('Bzmod', it, tmp)

    if (Dyn_quench) then
      if (.not.Damp) then
        call rescale_array(f(nxghost+1:nx-nxghost,3), tmp); call ts_data%set('alp_m', it, tmp)
      else
        call rescale_array(f(nxghost+1:nx-nxghost,7), tmp); call ts_data%set('alp_m', it, tmp)
      endif
    endif
    call rescale_array(h(nxghost+1:nx-nxghost), tmp); call ts_data%set('h', it, tmp)
    call rescale_array(Om(nxghost+1:nx-nxghost), tmp); call ts_data%set('Omega', it, tmp)
    call rescale_array(G(nxghost+1:nx-nxghost), tmp); call ts_data%set('Shear', it, tmp)
    call rescale_array(l(nxghost+1:nx-nxghost), tmp); call ts_data%set('l', it, tmp)
    call rescale_array(v(nxghost+1:nx-nxghost), tmp); call ts_data%set('v', it, tmp)
    call rescale_array(etat(nxghost+1:nx-nxghost), tmp); call ts_data%set('etat', it, tmp)
    call rescale_array(tau(nxghost+1:nx-nxghost), tmp); call ts_data%set('tau', it, tmp)
    call rescale_array(alp_k(nxghost+1:nx-nxghost), tmp); call ts_data%set('alp_k', it, tmp)
    call rescale_array(Uz(nxghost+1:nx-nxghost), tmp); call ts_data%set('Uz', it, tmp)
    call rescale_array(Ur(nxghost+1:nx-nxghost), tmp); call ts_data%set('Ur', it, tmp)
    call rescale_array(n(nxghost+1:nx-nxghost), tmp); call ts_data%set('n', it, tmp)
    call rescale_array(Beq(nxghost+1:nx-nxghost), tmp); call ts_data%set('Beq', it, tmp)
    call rescale_array(r_kpc(nxghost+1:nx-nxghost), tmp); call ts_data%set('rkpc', it, tmp)

    call rescale_array(Om_h(nxghost+1:nx-nxghost), tmp); call ts_data%set('Omega_h', it, tmp)
    call rescale_array(Om_b(nxghost+1:nx-nxghost), tmp); call ts_data%set('Omega_b', it, tmp)
    call rescale_array(Om_d(nxghost+1:nx-nxghost), tmp); call ts_data%set('Omega_d', it, tmp)

    if (ts_data%output('P')) then
      call rescale_array(P(nxghost+1:nx-nxghost), tmp); call ts_data%set('P', it, tmp)
      call rescale_array(Pd(nxghost+1:nx-nxghost), tmp); call ts_data%set('Pd', it, tmp)
      call rescale_array(Pm(nxghost+1:nx-nxghost), tmp); call ts_data%set('Pm', it, tmp)
      call rescale_array(Pstars(nxghost+1:nx-nxghost), tmp); call ts_data%set('Pstars', it, tmp)
      call rescale_array(Pbulge(nxghost+1:nx-nxghost), tmp); call ts_data%set('Pbulge', it, tmp)
      call rescale_array(Pdm(nxghost+1:nx-nxghost), tmp); call ts_data%set('Pdm', it, tmp)
      call rescale_array(P2(nxghost+1:nx-nxghost), tmp); call ts_data%set('P2', it, tmp)
    endif

    ! For convenience, computes and stores maximum magnetic field value, and the
    ! position of the maximum
    Btot = (ts_data%get_it('Br',it))**2 +  &
           (ts_data%get_it('Bp',it))**2 +  &
           (ts_data%get_it('Bzmod',it))**2
    Btot = sqrt(Btot)

    Bmax_idx = maxloc(Btot, 1)
    call ts_data%set_scalar('Bmax', it, Btot(Bmax_idx))
    tmp = ts_data%get_it('rkpc',it)
    call ts_data%set_scalar('rmax', it, tmp(Bmax_idx))
    ! Later, everything should be updated to accept integer datasets
    call ts_data%set_scalar('Bmax_idx', it, dfloat(Bmax_idx))

    ! alp is computed in the gutsdynamo module (annoyingly differently from
    ! anything else). Therefore, one needs to be careful. This is a good
    ! candidate for some code refactoring.
    if (status_code == 'M' .or. status_code == 'm') &
        call rescale_array(alp(nxghost+1:nx-nxghost), tmp); call ts_data%set('alp', it, tmp)

  end subroutine make_ts_arrays

  subroutine reset_ts_arrays()
    call ts_data%reset()
    ts_status_code = '-'
  end subroutine reset_ts_arrays

end module ts_arrays
