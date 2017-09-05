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
    call IO_write_dataset('Br', gal_id,                             &
                          ts_Br(:,:) * B0_mkG / B0,                 &
                          units='microgauss')

    call IO_write_dataset('Bp', gal_id,                             &
                          ts_Bp(:,:) * B0_mkG / B0,                 &
                          units='microgauss')

    if (Dyn_quench) then
      call IO_write_dataset('alp_m', gal_id,                        &
                            ts_alp_m(:,:)*h0_km/h0/t0_s*t0 ,        &
                            units='km/s')
    endif

    call IO_write_dataset('Bzmod', gal_id,                          &
                          ts_Bzmod(:,:) * B0_mkG / B0,              &
                          units='microgauss')

    call IO_write_dataset('h', gal_id,                              &
                          ts_h(:,:)/h0*h0_kpc*1d3,                  &
                          units='pc',                               &
                          description='Disk scaleheight')

    call IO_write_dataset('Omega',gal_id,                           &
                          ts_om(:,:)*h0_km/h0_kpc/t0_s*t0,          &
                          units='km/s/kpc',                         &
                          description='Angular velocity')

    call IO_write_dataset('Shear', gal_id,                          &
                          ts_G(:,:)*h0_km/h0_kpc/t0_s*t0,           &
                          units='km/s/kpc',                         &
                          description='Shear')

    call IO_write_dataset('l', gal_id,                              &
                          ts_l(:,:)/h0*h0_kpc*1d3,                  &
                          units='pc',                               &
                          description='Turbulent length')
!     call IO_write_dataset('v', gal_id, &
!                           ts_v(:,:)*h0_km/h0/t0_s*t0, &
!                           units='km/s', description='Turbulent velocity')

    call IO_write_dataset('etat', gal_id,                           &
                          ts_etat(:,:)*h0_km*h0_kpc/h0**2/t0_s*t0,  &
                          units='kpc km/s')

    call IO_write_dataset('tau', gal_id,ts_tau(:,:))

    call IO_write_dataset('alp_k', gal_id,                          &
                          ts_alp_k(:,:)*h0_km/h0/t0_s*t0 ,          &
                          units='km/s')

    call IO_write_dataset('Uz', gal_id,                             &
                          ts_Uz(:,:)*h0_km/h0/t0_s*t0,              &
                          units='km/s',                             &
                          description='Vertical velocity (outflow)')

!     call IO_write_dataset('Ur', gal_id,ts_Ur(:,:)*h0_km/h0/t0_s*t0, &
!                           units='km/s', description='Radial velocity')

    call IO_write_dataset('n', gal_id,                              &
                          ts_n(:,:)*n0_cm3/n0,                      &
                          units='cm^-3',                            &
                          description='Number density of diffuse gas')

    call IO_write_dataset('Beq', gal_id,                            &
                          ts_Beq(:,:)* B0_mkG / B0,                 &
                          units='microgauss',                       &
                          description='Equipartition field')

    call IO_write_dataset('rmax', gal_id,                           &
                          ts_rmax/h0*h0_kpc,                        &
                          units='kpc',                              &
                          description='Radius of maximum B_phi')

    call IO_write_dataset('alp', gal_id,                            &
                          ts_alp(:,:)*h0_km/h0/t0_s*t0 ,            &
                          units='km/s')
    call IO_write_dataset('dt', gal_id,ts_dt)

    call IO_write_dataset('r', gal_id,                              &
                          ts_rkpc,                                  &
                          units='kpc',                              &
                          description='Radius')

    if (p_extra_rotation_curve_outputs) then
      call IO_write_dataset('Omega_h', gal_id,                         &
                            ts_Om_h,                                   &
                            units='km/s/kpc',                          &
                            description='Angular velocity of the halo' &
                            // ' (no regularization).')
      call IO_write_dataset('Omega_b', gal_id,                         &
                            ts_Om_b,                                   &
                            units='km/s/kpc',                          &
                            description='Angular velocity of the bulge' &
                            // ' (no regularization).')
      call IO_write_dataset('Omega_d', gal_id,                         &
                            ts_Om_d,                                   &
                            units='km/s/kpc',                          &
                            description='Angular velocity of the disc' &
                            // ' (no regularization).')
    endif

    if (p_extra_pressure_outputs) then
      call IO_write_dataset('P', gal_id, ts_P, units='erg cm^-3', &
                            description='Pressure in the midplane.')
      call IO_write_dataset('Pd', gal_id, ts_Pd, units='erg cm^-3',           &
                            description='Contribution of the gravity of the ' &
                            // 'diffuse gas to the pressure in the midplane')
      call IO_write_dataset('Pm', gal_id, ts_Pm, units='erg cm^-3',           &
                            description='Contribution of the gravity of the ' &
                            // 'molecular gas to the pressure in the midplane')
      call IO_write_dataset('Pstars', gal_id, ts_Pstars, units='erg cm^-3',   &
                            description='Contribution of the gravity of the ' &
                            // 'stars to the pressure in the midplane')
      call IO_write_dataset('Pbulge', gal_id, ts_Pbulge, units='erg cm^-3',   &
                            description='Contribution of the gravity of the ' &
                            // 'stellar bulge to the pressure in the midplane')
      call IO_write_dataset('Pdm', gal_id, ts_Pdm, units='erg cm^-3',            &
                            description='Contribution of the gravity of the '    &
                            // 'dark matter halo to the pressure in the midplane')
      call IO_write_dataset('P2', gal_id, ts_P2, units='erg cm^-3',  &
                            description='Correction to the pressure ' &
                            // ' contribution associated with the fact that ' &
                            // 'the disc is not perfectly thin')
    endif

    call IO_finish_galaxy(gal_id)
    
    ! Writes specific values to screen for a quick check
    call message('ts_Br(1,1)'   ,ts_Br(1,1), info=4)
    call message('ts_Br(1,37)'  ,ts_Br(1,37), info=4)
    call message('ts_Br(n1,1)' ,ts_Br(n1,1), info=4)
    call message('ts_Br(n1,37)',ts_Br(n1,37), info=4)
    
  end subroutine write_output

end module output
