!*****************************************************
program calldynamo
  use dynamo
  use input_params
!
  implicit none 
!
  integer :: igal, flag
! intent(in)
!  info:   how much info to print during run
!          options are 0=none, 1=min, 2=standard, 3=max 
!  gal_id: id number of galaxy
!          ranges from 1 to ~300,000
!
! intent(out)
!  flag:   status of run
!          results are -1=run successful, 1=, 2=, 3=
!
! Note: physical variables (ts_t, ts_Br, ts_Bp, ts_Bzmod, ts_alp_m) are written to files, not passed
!
! Format:  call dynamo_run(info, gal_id, flag)
!  Standard parameters: (2, igal, flag), where igal is looped over within python
!
  call dynamo_run(2, igal, flag)
!
end program calldynamo
!*****************************************************
