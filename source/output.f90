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
! Contains a module that calls the IO routines to output the results
module output
  use IO
  implicit none

contains
  subroutine write_output(gal_id, runtime)
    ! Writes the output
    use global_input_parameters
    use input_params
    use profiles
    use bzcalc
    use ts_arrays
    use units
    use messages
    integer, intent(in) :: gal_id
    double precision, intent(in) :: runtime

    ! Writes log data
    call IO_write_dataset('runtime', gal_id,                        &
                          [runtime],                                &
                          units='s',                                &
                          description='Running time of the galaxy', &
                          is_log=.true.)

    call IO_write_dataset('nsteps', gal_id,                         &
                          [dble(nsteps)],                           &
                          description='Number of timesteps',        &
                          is_log=.true.)

    call IO_write_dataset('status', gal_id,                         &
                          ts_status_code,                           &
                          description='Status codes.',              &
                          is_log=.true.)

    call IO_write_dataset('completed', gal_id,                      &
                          [1d0], is_log=.true.)

    ! Outputs the data!
    call IO_write_dataset('r', gal_id,                            &
                          ts_data%get('rkpc'),                    &
                          units='kpc',                            &
                          description='Radius')


    if (ts_data%output('Br')) &
      call IO_write_dataset('Br', gal_id,                           &
                          ts_data%get('Br') * B0_mkG / B0,          &
                          units='microgauss')

    if (ts_data%output('Bp')) &
      call IO_write_dataset('Bp', gal_id,                           &
                          ts_data%get('Bp') * B0_mkG / B0,          &
                          units='microgauss')

    if (Dyn_quench) then
      if (ts_data%output('alp_m')) &
        call IO_write_dataset('alp_m', gal_id,                      &
                            ts_data%get('alp_m')*h0_km/h0/t0_s*t0 , &
                            units='km/s')
    endif

    if (ts_data%output('Bzmod')) &
      call IO_write_dataset('Bzmod', gal_id,                        &
                          ts_data%get('Bzmod') * B0_mkG / B0,       &
                          units='microgauss')

    if (ts_data%output('h')) &
      call IO_write_dataset('h', gal_id,                            &
                          ts_data%get('h')/h0*h0_kpc*1d3,           &
                          units='pc',                               &
                          description='Disk scaleheight')

    if (ts_data%output('Omega')) &
      call IO_write_dataset('Omega',gal_id,                         &
                          ts_data%get('Omega')*h0_km/h0_kpc/t0_s*t0,&
                          units='km/s/kpc',                         &
                          description='Angular velocity')

    if (ts_data%output('Shear')) &
      call IO_write_dataset('Shear', gal_id,                        &
                          ts_data%get('Shear')*h0_km/h0_kpc/t0_s*t0,&
                          units='km/s/kpc',                         &
                          description='Shear')

    if (ts_data%output('l')) &
      call IO_write_dataset('l', gal_id,                            &
                          ts_data%get('l')/h0*h0_kpc*1d3,           &
                          units='pc',                               &
                          description='Turbulent length')
    if (ts_data%output('v')) &
    call IO_write_dataset('v', gal_id, &
                          ts_data%get('v')*h0_km/h0/t0_s*t0, &
                          units='km/s', description='Turbulent velocity')

    if (ts_data%output('etat')) &
      call IO_write_dataset('etat', gal_id,                              &
                          ts_data%get('etat')*h0_km*h0_kpc/h0**2/t0_s*t0,&
                          units='kpc km/s')

    if (ts_data%output('tau')) &
      call IO_write_dataset('tau', gal_id,ts_data%get('tau'))

    if (ts_data%output('alp_k')) &
      call IO_write_dataset('alp_k', gal_id,                      &
                          ts_data%get('alp_k')*h0_km/h0/t0_s*t0 , &
                          units='km/s')

    if (ts_data%output('Uz')) &
      call IO_write_dataset('Uz', gal_id,                           &
                          ts_data%get('Uz')*h0_km/h0/t0_s*t0,       &
                          units='km/s',                             &
                          description='Vertical velocity (outflow)')

!     if (ts_data%output('Ur')) &
!        call IO_write_dataset('Ur', gal_id,ts_data%get('Ur')*h0_km/h0/t0_s*t0, &
!                           units='km/s', description='Radial velocity')

    if (ts_data%output('n')) &
      call IO_write_dataset('n', gal_id,                              &
                          ts_data%get('n')*n0_cm3/n0,                 &
                          units='cm^-3',                              &
                          description='Number density of diffuse gas')

    if (ts_data%output('Beq')) &
      call IO_write_dataset('Beq', gal_id,                          &
                          ts_data%get('Beq')* B0_mkG / B0,          &
                          units='microgauss',                       &
                          description='Equipartition field')

    if (ts_data%output('rmax')) &
      call IO_write_dataset('rmax', gal_id,                         &
                          ts_data%get_scalar('rmax')/h0*h0_kpc,     &
                          units='kpc',                              &
                          description='Radius of maximum |B|')

    if (ts_data%output('Bmax')) &
      call IO_write_dataset('Bmax', gal_id,                         &
                          ts_data%get_scalar('Bmax')* B0_mkG / B0,  &
                          units='microgauss',                       &
                          description='Maximum B')

    if (ts_data%output('Bmax_idx')) &
      call IO_write_dataset('Bmax_idx', gal_id,                     &
                          ts_data%get_scalar('Bmax_idx'),           &
                          description='Index of maximum |B|')

    if (ts_data%output('alp')) &
      call IO_write_dataset('alp', gal_id,                          &
                          ts_data%get('alp')*h0_km/h0/t0_s*t0 ,     &
                          units='km/s')
    if (ts_data%output('dt')) &
      call IO_write_dataset('dt', gal_id,ts_data%get_scalar('dt'))


    if (ts_data%output('Omega_h')) &
      call IO_write_dataset('Omega_h', gal_id,                       &
                          ts_data%get('Omega_h'),                    &
                          units='km/s/kpc',                          &
                          description='Angular velocity of the halo' &
                          // ' (no regularization).')
    if (ts_data%output('Omega_b')) &
      call IO_write_dataset('Omega_b', gal_id,                       &
                          ts_data%get('Omega_b'),                    &
                          units='km/s/kpc',                          &
                          description='Angular velocity of the bulge'&
                          // ' (no regularization).')
    if (ts_data%output('Omega_d')) &
      call IO_write_dataset('Omega_d', gal_id,                       &
                          ts_data%get('Omega_d'),                    &
                          units='km/s/kpc',                          &
                          description='Angular velocity of the disc' &
                          // ' (no regularization).')


    if (ts_data%output('P')) &
      call IO_write_dataset('P', gal_id, ts_data%get('P'), units='erg cm^-3', &
                            description='Pressure in the midplane.')
    if (ts_data%output('Pd')) &
      call IO_write_dataset('Pd', gal_id, ts_data%get('Pd'), units='erg cm^-3',           &
                            description='Contribution of the gravity of the ' &
                            // 'diffuse gas to the pressure in the midplane')
    if (ts_data%output('Pm')) &
      call IO_write_dataset('Pm', gal_id, ts_data%get('Pm'), units='erg cm^-3',           &
                            description='Contribution of the gravity of the ' &
                            // 'molecular gas to the pressure in the midplane')
    if (ts_data%output('Pstars')) &
      call IO_write_dataset('Pstars', gal_id, ts_data%get('Pstars'), units='erg cm^-3',   &
                            description='Contribution of the gravity of the ' &
                            // 'stars to the pressure in the midplane')
    if (ts_data%output('Pbulge')) &
      call IO_write_dataset('Pbulge', gal_id, ts_data%get('Pbulge'), units='erg cm^-3',   &
                            description='Contribution of the gravity of the ' &
                            // 'stellar bulge to the pressure in the midplane')
    if (ts_data%output('Pdm')) &
      call IO_write_dataset('Pdm', gal_id, ts_data%get('Pdm'), units='erg cm^-3',            &
                            description='Contribution of the gravity of the '    &
                            // 'dark matter halo to the pressure in the midplane')
    if (ts_data%output('P2')) &
      call IO_write_dataset('P2', gal_id, ts_data%get('P2'), units='erg cm^-3',  &
                            description='Correction to the pressure ' &
                            // ' contribution associated with the fact that ' &
                            // 'the disc is not perfectly thin')

    call IO_finish_galaxy(gal_id)

  end subroutine write_output

end module output
