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

    integer, intent(in) :: gal_id
    
    ! Writes any meta information about the run
    call IO_start_galaxy(gal_id, info)
    if (info>2) print *, 'Galaxy',gal_id,' -> IO initialised'
    
    call IO_write_dataset('t_Gyr', gal_id, info,ts_t_Gyr)
    call IO_write_dataset('Br', gal_id, info,ts_Br(:,:))
    call IO_write_dataset('Bp', gal_id, info,ts_Bp(:,:))
    if (Dyn_quench) then
      call IO_write_dataset('alp_m', gal_id, info,ts_alp_m(:,:))
    endif
    call IO_write_dataset('Bzmod', gal_id, info,ts_Bzmod(:,:))
    call IO_write_dataset('h', gal_id, info,ts_h(:,:))
    call IO_write_dataset('om', gal_id, info,ts_om(:,:))
    call IO_write_dataset('G', gal_id, info,ts_G(:,:))
    call IO_write_dataset('l', gal_id, info,ts_l(:,:))
    call IO_write_dataset('v', gal_id, info,ts_v(:,:))
    call IO_write_dataset('etat', gal_id, info,ts_etat(:,:))
    call IO_write_dataset('tau', gal_id, info,ts_tau(:,:))
    call IO_write_dataset('alp_k', gal_id, info,ts_alp_k(:,:))
    call IO_write_dataset('Uz', gal_id, info,ts_Uz(:,:))
    call IO_write_dataset('Ur', gal_id, info,ts_Ur(:,:))
    call IO_write_dataset('n', gal_id, info,ts_n(:,:))
    call IO_write_dataset('Beq', gal_id, info, ts_Beq(:,:))
    call IO_write_dataset('rmax', gal_id, info,ts_rmax)
    call IO_write_dataset('delta_r', gal_id, info,ts_delta_r)
    call IO_write_dataset('alp', gal_id, info,ts_alp(:,:))

    call IO_write_dataset('dt', gal_id, info,ts_dt)
    
    call IO_write_dataset('r', gal_id, info, ts_r)
    
    if (info>2) print *, 'Galaxy',gal_id,' -> datasets written'

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
