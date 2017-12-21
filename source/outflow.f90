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
! Contains subroutines which compute the outflows
module outflow
  implicit none
  private
  public outflow_speed
contains
  function outflow_speed(r, rho, h, vt, rdisk, vdisk, SFR, Rm, outflow_type)
    ! Computes the outflow speed.
    ! There are currently 5 prescriptions available:
    !   1) no_outflow -> vout=0
    !   2) vturb -> assumes is always vout=vt
    !   3) wind -> compute vout from Galform's mass outflow rate
    !   4) superbubble -> uses the superbubble formalism
    !   5) superbubble_simple -> the same, without the \rho_{hot}/\rho step
    !
    ! Input: r -> radii array, in kpc
    !        rho -> densities array, in g/cm^3
    !        h -> heights array, in kpc
    !        vt -> turbulent velocities array, in km/s
    !        rdisk -> half mass radius of the disk, in kpc
    !        vdisk -> velocity at the half mass radius of the disk, in km/s
    !        SFR -> total star formation rate of the galaxy, in (solar masses)/yr
    !        Rm -> molecular/diffuse gas ratio
    !        Mgas -> gas mass of the disk, in solar masses
    !        Mstars -> stellar mass of the disk, in solar masses
    ! Output: array containing the outflow speed at different radii

    use input_constants
    use global_input_parameters
    use math_constants
    implicit none
    double precision, intent(in) :: rdisk, vdisk, SFR
    double precision, dimension(:), intent(in) :: r, rho, h, vt, Rm
    double precision, dimension(size(r)) :: outflow_speed, n, fm
    double precision :: constant, rs, beta
    character(len=*), optional, intent(in) :: outflow_type
    character(len=35) :: outflow_type_actual
    logical :: average
    double precision, parameter :: VDISK_MIN = 0.01 ! km/s

    if (present(outflow_type)) then
      outflow_type_actual = outflow_type
    else
      outflow_type_actual = p_outflow_type
    endif

    select case (trim(outflow_type_actual))
      case('no_outflow')
        outflow_speed = r*0.0
        return
      case('vturb')
        outflow_speed = vt
        return
      case('wind')
        ! Computes the molecular fraction
        fm = Rm/(1d0+Rm)
        ! Computes the mass loading
        ! (traps possible problems assigning a minimum rotation velocity)
        beta = (max(VDISK_MIN,vdisk)/p_outflow_vhot)**p_outflow_alphahot
        ! This case assumes the outflow rates used in galform are correct
        ! This, together with a Blitz&Rosolowsky star formation law leads
        ! to:
        outflow_speed = 0.5 * p_outflow_nu0 * beta * h * fm
        outflow_speed = outflow_speed * km_kpc / s_Gyr
        return
        ! Note: at the moment, nu0 and alpha_hot and Vhot are input
        ! parameters, but later these will be read from the input hdf5 files
        ! (which in turn contains the parameters used in that particular
        ! galform run).
      case('superbubble_simple')
        average=.false.
      case('superbubble')
        average=.true.
      case default
        print *, outflow_type_actual
        stop 'outflow: Unknown outflow_type'
    end select

    ! The following is for the superbubble case

    ! First, prepares some constants
    rs = rdisk*constDiskScaleToHalfMassRatio
    n = rho/Hmass

    constant = 51.3 * (p_outflow_Lsn/1d38)**(1./3.) * p_outflow_fOB/0.7d0 &
        * (p_outflow_etaSN/9.4e-3) * (p_tOB/3d0) * (40d0/p_N_SN1OB)

    where (n>0d0)
      outflow_speed = constant * (rs/3.0)**(-2) * SFR
      outflow_speed = outflow_speed * n**(-1./3.)
      outflow_speed = outflow_speed * (h/0.2)**(4./3.)
      outflow_speed = outflow_speed * exp(-r/rs)
    elsewhere
      ! Sets outflow to zero in the case of zero density
      ! (n=0 and h<0 are used to flag errors in the pressure calculation)
      outflow_speed = 0d0
    endwhere

    if (average) then
      where (n>0d0)
        outflow_speed = (p_outflow_hot_gas_density/rho) * outflow_speed
      endwhere
    endif

  end function outflow_speed
end module outflow
