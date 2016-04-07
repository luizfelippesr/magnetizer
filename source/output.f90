! Contains a module that calls the IO routines to output the results
module output
  use IO
  implicit none
  
contains
  subroutine write_output(gal_id)
    ! Writes the output
    use global_input_parameters
    use galaxy_model
    use bzcalc
    use ts_arrays
    use units

    integer, intent(in) :: gal_id
    
    ! Writes any meta information about the run
    call IO_start_galaxy(gal_id, info)
    if (info>2) print *, 'Galaxy',gal_id,' -> IO initialised'

!     call IO_write_dataset('t', gal_id, info,                     &
!                           ts_t_Gyr,                              &
!                           units='Gyr',                           &
!                           description='Time since the Big Bang.')

    call IO_write_dataset('Br', gal_id, info,                    &
                          ts_Br(:,:) * B0_mkG / B0,              &
                          units='microgauss')

    call IO_write_dataset('Bp', gal_id, info,                    &
                          ts_Bp(:,:) * B0_mkG / B0,              &
                          units='microgauss')

    if (Dyn_quench) then
      call IO_write_dataset('alp_m', gal_id, info,               &
                            ts_alp_m(:,:)*h0_km/h0/t0_s*t0 ,     &
                            units='km/s')
    endif

    call IO_write_dataset('Bzmod', gal_id, info,                 &
                          ts_Bzmod(:,:) * B0_mkG / B0,           &
                          units='microgauss')

    call IO_write_dataset('h', gal_id, info,                     &
                          ts_h(:,:)/h0*h0_kpc*1d3,               &
                          units='pc',                            &
                          description='Disk scaleheight')

    call IO_write_dataset('Omega',gal_id,info,                   &
                          ts_om(:,:)*h0_km/h0_kpc/t0_s*t0,       &
                          units='km/s/kpc',                      &
                          description='Angular velocity')

    call IO_write_dataset('Shear', gal_id, info,                 &
                          ts_G(:,:)*h0_km/h0_kpc/t0_s*t0,        &
                          units='km/s/kpc',                      &
                          description='Shear')

    call IO_write_dataset('l', gal_id, info,                     &
                          ts_l(:,:)/h0*h0_kpc*1d3,               &
                          units='pc',                            &
                          description='Turbulent length')
!     call IO_write_dataset('v', gal_id, info, &
!                           ts_v(:,:)*h0_km/h0/t0_s*t0, &
!                           units='km/s', description='Turbulent velocity')

    call IO_write_dataset('etat', gal_id, info,                     &
                          ts_etat(:,:)*h0_km*h0_kpc/h0**2/t0_s*t0,  &
                          units='kpc km/s')

    call IO_write_dataset('tau', gal_id, info,ts_tau(:,:))

    call IO_write_dataset('alp_k', gal_id, info,                    &
                          ts_alp_k(:,:)*h0_km/h0/t0_s*t0 ,          &
                          units='km/s')

    call IO_write_dataset('Uz', gal_id, info,                       &
                          ts_Uz(:,:)*h0_km/h0/t0_s*t0,              &
                          units='km/s',                             &
                          description='Vertical velocity (outflow)')

!     call IO_write_dataset('Ur', gal_id, info,ts_Ur(:,:)*h0_km/h0/t0_s*t0, &
!                           units='km/s', description='Radial velocity')

    call IO_write_dataset('n', gal_id, info,                        &
                          ts_n(:,:)*n0_cm3/n0,                      &
                          units='cm^-3',                            &
                          description='Number density of diffuse gas')

    call IO_write_dataset('Beq', gal_id, info,                      &
                          ts_Beq(:,:)* B0_mkG / B0,                 &
                          units='microgauss',                       &
                          description='Equipartition field')

    call IO_write_dataset('rmax', gal_id, info,                     &
                          ts_rmax/h0*h0_kpc,                        &
                          units='kpc',                              &
                          description='Radius of maximum B_phi')

    call IO_write_dataset('delta_r', gal_id, info,ts_delta_r)

    call IO_write_dataset('alp', gal_id, info,                      &
                          ts_alp(:,:)*h0_km/h0/t0_s*t0 ,            &
                          units='km/s')
    call IO_write_dataset('dt', gal_id, info,ts_dt)

    call IO_write_dataset('r', gal_id, info,                        &
                          ts_rkpc,                                  &
                          units='kpc',                              &
                          description='Radius')

    if (info>2) print *, 'Galaxy  ',gal_id,' -> datasets written'
    print *, 'Galaxy  ',gal_id,' -> datasets written'
    print *, 'Galaxy  ',gal_id,' -> datasets written'
    print *, 'Galaxy  ',gal_id,' -> datasets written'
    print *, 'Galaxy  ',gal_id,' -> datasets written'
    print *, 'Galaxy  ',gal_id,' -> datasets written'
    call IO_finish_galaxy(gal_id, info)
    
    ! Writes specific values to screen for a quick check
    if (info > 2) then
      print *,'ts_Br(1,1)'   ,ts_Br(1,1)
      print *,'ts_Br(1,37)'  ,ts_Br(1,37)
      print *,'ts_Br(n1,1)' ,ts_Br(n1,1)
      print *,'ts_Br(n1,37)',ts_Br(n1,37)
    endif
    
  end subroutine write_output
  
end module output
