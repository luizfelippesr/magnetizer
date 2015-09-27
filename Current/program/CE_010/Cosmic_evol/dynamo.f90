!*****************************************************
!This is the top level dynamo program.
!It contains the subroutine dynamo_run, which solves the dynamo equations.
!*****************************************************
module dynamo
  use ts_arrays
  use start
  use timestep
  use output

  implicit none 

  integer :: it=0, jt=0
  double precision :: cpu_time_start, cpu_time_finish

  contains
    subroutine dynamo_run(info, gal_id, flag)

      integer, intent(in) :: info, gal_id
      integer, intent(out) :: flag
      integer :: fail_count = 0
      character(len=8) :: frmt
      character(len=8) :: gal_id_string

      frmt='(I8.8)'
      write(gal_id_string,frmt) gal_id

      call cpu_time(cpu_time_start)

100 continue  !LFSR: what is this? A goto is REALLY necessary?
      call init_start(gal_id_string,info)  !Initialize galaxy model and magnetic field
      call reset_input_params() ! In case there was another run before this call

      call make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)  !Store simulation output into arrays

      if (info> 1) then
        print*,'galaxy id=',gal_id
        print*,'it=',it,'t=',t,'t(Gyr)=',t*t0_Gyr,'r(R/2)(kpc)=',x(4+nxphys/2)*r_max_kpc
        print*,'Br(R/2)(mkG)=',f(4+nxphys/2,1),'Bp(R/2)(mkG)=',f(4+nxphys/2,2),'|Bz(R/2)|(mkG)=',Bzmod(4+nxphys/2)
        print*,'p(R/2)(deg)=',180./pi*atan(f(4+nxphys/2,1)/f(4+nxphys/2,2))
        print*, "cpu time in seconds: ", (cpu_time_finish -cpu_time_start)
        print*,'Directory ending   =', model_name
        print*,
      endif
!     TIME-STEPPING
      do it=1,n1  !At every iteration it, the output is written to a file
        do jt=1,n2  !At every iteration jt, a time-step is calculated
!          if (mod((it-1)*n2+jt,nread) == 0) then  !If it*jt is an integer multiple of nread, new input params are read in and a new gal model calcd
!            print*,'timestep (it-1)*n2+jt=',(it-1)*n2+jt,'it=',it
!            call construct_galaxy_model(gal_id_string,info)
!          endif
          call rk(f)  !Run Runge-Kutta time-stepping routine
          if (isnan(f(4+nxphys/2,1)) .or. maxval(f)>1.d10) then  !Change time step if magnetic field blows up
            call cpu_time(cpu_time_finish)
            print*, "cpu time in seconds: ", (cpu_time_finish -cpu_time_start)
            !stop 'NANs detected, change time step'
            print*, 'NANs or unrealistically large magnetic fields detected, changing time step'
            eps_t=0.5*eps_t
            ! Reset reading of input parameters
            call reset_input_params()
            t=0.d0 !reset time
            first=0.d0 !reset timestepping
            deallocate(f) !Reset variables
            deallocate(dfdt) !Reset variable time derivs
            fail_count = fail_count + 1
            if (fail_count > 10) then
              flag=1
              return
            endif
            goto 100 ! LFSR: we would probably be better removing this...
          endif
        enddo

        call impose_bc(f)  !Impose boundary conditions before writing output

        call estimate_Bzmod(f)  !Estimate the value of |B_z| using Div B =0 condition

        call make_ts_arrays(it,t,f,Bzmod,h,om,G,l,v,etat,tau,alp_k,alp,Uz,Ur,n,Beq,rmax,delta_r)  !Store simulation output into arrays

        if (mod(it,nscreen) == 0) then  !Print some of the output to screen/diagnostic file while simulation runs
            print*,''
            print*,'galaxy_id=',gal_id
            print*,'it=',it,'t=',t,'t(Gyr)=',t*t0_Gyr,'r(R/2)(kpc)=',x(4+nxphys/2)*r_max_kpc
          if (info> 1) then
            print*,'Br(R/2)(mkG)=',f(4+nxphys/2,1),'Bp(R/2)(mkG)=',f(4+nxphys/2,2),'|Bz(R/2)|(mkG)=',Bzmod(4+nxphys/2)
          endif
          print*,'p(R/2)(deg)=',180./pi*atan(f(4+nxphys/2,1)/f(4+nxphys/2,2))
          if (info> 1) then
            print*, "cpu time in seconds: ", (cpu_time_finish -cpu_time_start)
            print*,'Directory ending   =', model_name
          endif
          print*,
          open(20,file= 'diagnostic.out',status="old",position="append")
          write(20,*),'it=',it,'t=',t,'t(Gyr)=',t*t0_Gyr,'r(R/2)(kpc)=',x(4+nxphys/2)*r_max_kpc
          write(20,*),'Br(R/2)(mkG)=',f(4+nxphys/2,1),'Bp(R/2)(mkG)=',f(4+nxphys/2,2),'|Bz(R/2)|(mkG)=',Bzmod(4+nxphys/2)
          write(20,*),'p(R/2)(deg)=',180./pi*atan(f(4+nxphys/2,1)/f(4+nxphys/2,2))
          write(20,*), "cpu time in seconds: ", (cpu_time_finish -cpu_time_start)
          write(20,*)'Directory ending   =', model_name
          write(20,*)
          close(20)
        endif
        if (mod(it,nread) == 0) then  !If it*jt is an integer multiple of nread, new input params are read in and a new gal model calcd
          print*,''
          print*,'Uz_sol_kms=',Uz_sol_kms
          call construct_galaxy_model(gal_id_string,info)
          print*,'it=',it
          print*,'Uz_sol_kms=',Uz_sol_kms
          print*,''
        endif

        if (last_output) exit
      enddo

      call write_output(f,gal_id)  !Writes final simulation output

      call cpu_time(cpu_time_finish)
      if (info>0) then
        print*, ''
        print*, 'galaxy_id ', gal_id_string,' simulation cpu time in seconds: ', (cpu_time_finish -cpu_time_start)
        print*, ''
      endif
      open(20,file= 'diagnostic.out',status="old",position="append")
      write(20,*) "simulation time in seconds: ", (cpu_time_finish -cpu_time_start)
      close(20)

      flag=-1  !Tell calldynamo that simulation was successful
      iread=0  !Reset iread
    end subroutine dynamo_run
end module dynamo
!*****************************************************
