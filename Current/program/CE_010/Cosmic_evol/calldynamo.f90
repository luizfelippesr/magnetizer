!This program runs the code in fortran without mpi
!*****************************************************
program calldynamo
  use dynamo
  use input_params
!
  implicit none 
!
  integer :: info=2, igal=1, flag, ngal=1
! intent(in)
!  info:   how much info to print during run
!          options are 0=min, 1=standard, 2=max
!  gal_id: id number of galaxy
!          ranges from 1 to ~300,000
!
! intent(out)
!  flag:   status of run
!          results are -1=run successful, 1=, 2=, 3=
!
! Note: calculated physical variables (ts_t, ts_Br, ts_Bp, ts_Bzmod, ts_alp_m) are written to files, not passed
!
! Format:  call dynamo_run(info, gal_id, flag)
! Standard parameters: (0, 1, flag), where igal is looped over from 1 to ngal
!
do igal=1,ngal
  call dynamo_run(info, igal, flag)
  if (info>1) then
    print*,'flag=',flag
    print*,''
  endif
enddo
!
end program calldynamo
!*****************************************************
