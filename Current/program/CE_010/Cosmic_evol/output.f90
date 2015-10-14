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
    print *, 'Galaxy',gal_id,' -> IO initialised' ! TO BE REMOVED
    
    call IO_write_dataset('t', gal_id, info,ts_t)
!     call IO_write_dataset('dt', gal_id, info,ts_dt)
    call IO_write_dataset('Br', gal_id, info,ts_Br(:iread-1,:))
    call IO_write_dataset('Bp', gal_id, info,ts_Bp(:iread-1,:))
    if (Dyn_quench) then
      call IO_write_dataset('alp_m', gal_id, info,ts_alp_m(:iread-1,:))
    endif
    call IO_write_dataset('Bzmod', gal_id, info,ts_Bzmod(:iread-1,:))
    call IO_write_dataset('h', gal_id, info,ts_h(:iread-1,:))
    call IO_write_dataset('om', gal_id, info,ts_om(:iread-1,:))
    call IO_write_dataset('G', gal_id, info,ts_G(:iread-1,:))
    call IO_write_dataset('l', gal_id, info,ts_l(:iread-1,:))
    call IO_write_dataset('v', gal_id, info,ts_v(:iread-1,:))
    call IO_write_dataset('etat', gal_id, info,ts_etat(:iread-1,:))
    call IO_write_dataset('tau', gal_id, info,ts_tau(:iread-1,:))
    call IO_write_dataset('alp_k', gal_id, info,ts_alp_k(:iread-1,:))
    call IO_write_dataset('Uz', gal_id, info,ts_Uz(:iread-1,:))
    call IO_write_dataset('Ur', gal_id, info,ts_Ur(:iread-1,:))
    call IO_write_dataset('n', gal_id, info,ts_n(:iread-1,:))
    call IO_write_dataset('Beq', gal_id, info, ts_Beq(:iread-1,:))
    call IO_write_dataset('rmax', gal_id, info,ts_rmax)
    call IO_write_dataset('delta_r', gal_id, info,ts_delta_r)
    call IO_write_dataset('alp', gal_id, info,ts_alp(:iread-1,:))
    print *, 'Galaxy',gal_id,' -> datasets written'  ! TO BE REMOVED
    call IO_finish_galaxy(gal_id, info)
    
    ! Writes specific values to screen for a quick check
    if (info > 2) then
      print *,'ts_Br(1,1)'   ,ts_Br(1,1)
      print *,'ts_Br(1,37)'  ,ts_Br(1,37)
      print *,'ts_Br(n1+1,1)' ,ts_Br(n1+1,1)
      print *,'ts_Br(n1+1,37)',ts_Br(n1+1,37)
      print *,'ts_Br(10,10)',ts_Br(10,10)
    endif
    
  end subroutine write_output
  
end module output
