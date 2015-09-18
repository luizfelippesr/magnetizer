module output_dump
  implicit none
  
contains
  subroutine write_final_output(f,gal_id_string)

    use global_input_parameters
    use galaxy_model
    use bzcalc
    use ts_arrays

    integer, parameter :: output_unit = 12
    integer, parameter :: ts_unit = 13
    integer, parameter :: diag_unit = 20
    double precision, dimension(nx,nvar), intent(in) :: f
    character (len=8), intent(in) :: gal_id_string

    !     WRITE DATA FOR FINAL TIMESTEP TO FILE "run.out"
    if (info> 0) then
      print*,'Writing output for final timestep to file run',gal_id_string,'.out'
    endif
    open(unit=diag_unit,file= 'diagnostic.out',status="old",position="append")
    write(diag_unit,*)'Writing output for final timestep to file run',gal_id_string,'.out'
    close(diag_unit)

    open(unit=output_unit,file= trim(path_to_input_directories) // '/output/' // trim(model_name) &
                  // '/run_' // gal_id_string // '.out',status="replace")
    write(output_unit,*) t
    write(output_unit,*) f(:,1)
    write(output_unit,*) f(:,2)
    if (Dyn_quench) then
      if (.not.Damp) then
        write(output_unit,*) f(:,3)
      else
        write(output_unit,*) f(:,7)
      endif
    endif
    write(output_unit,*) Bzmod
    write(output_unit,*) h
    write(output_unit,*) om
    write(output_unit,*) G
    write(output_unit,*) l
    write(output_unit,*) v
    write(output_unit,*) etat
    write(output_unit,*) tau
    write(output_unit,*) alp_k
    write(output_unit,*) Uz
    write(output_unit,*) Ur
    write(output_unit,*) n
    write(output_unit,*) Beq
    close(output_unit)
!
!     WRITE DATA FOR ALL TIMESTEPS TO FILE "ts.out"
    if (info> 0) then
      print*,'Writing time series output to file ts',gal_id_string,'.out'
    endif
    open(unit=ts_unit,file= trim(path_to_input_directories) // '/output/' // &
                  trim(model_name) // '/ts_' // gal_id_string // '.out',status="replace")
    write(ts_unit,*) ts_t(:iread-1)
    write(ts_unit,*) ts_Br(:iread-1,:)
    write(ts_unit,*) ts_Bp(:iread-1,:)
    if (Dyn_quench) then
      if (.not.Damp) then
        write(ts_unit,*) ts_alp_m(:iread-1,:)
      else
        write(ts_unit,*) ts_alp_m(:iread-1,:)
      endif
    endif
    write(ts_unit,*) ts_Bzmod(:iread-1,:)
    write(ts_unit,*) ts_h(:iread-1,:)
    write(ts_unit,*) ts_om(:iread-1,:)
    write(ts_unit,*) ts_G(:iread-1,:)
    write(ts_unit,*) ts_l(:iread-1,:)
    write(ts_unit,*) ts_v(:iread-1,:)
    write(ts_unit,*) ts_etat(:iread-1,:)
    write(ts_unit,*) ts_tau(:iread-1,:)
    write(ts_unit,*) ts_alp_k(:iread-1,:)
    write(ts_unit,*) ts_Uz(:iread-1,:)
    write(ts_unit,*) ts_Ur(:iread-1,:)
    write(ts_unit,*) ts_n(:iread-1,:)
    write(ts_unit,*) ts_Beq(:iread-1,:)
    write(ts_unit,*) ts_rmax
    write(ts_unit,*) ts_delta_r
    write(ts_unit,*) ts_alp(:iread-1,:)
    close(ts_unit)
    print*,'ts_Br(1,1)'   ,ts_Br(1,1)
    print*,'ts_Br(1,37)'  ,ts_Br(1,37)
    print*,'ts_Br(n1+1,1)' ,ts_Br(n1+1,1)
    print*,'ts_Br(n1+1,37)',ts_Br(n1+1,37)
    print*,'ts_Br(10,10)',ts_Br(10,10)
  end subroutine write_final_output
end module output_dump
