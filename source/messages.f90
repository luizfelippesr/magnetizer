! Contains a module that implements a messaging and logging
! str function adapted from the fortran-utils module
! https://github.com/certik/fortran-utils/blob/master/src/utils.f90
!
module messages
  implicit none
  private

  integer :: MPI_rank = -1
  integer :: verbosity_setting = 10
  public message, str, str2int


  interface str
      module procedure str_int, str_dp, str_dp_n
  end interface

  contains

  subroutine message(msg, val, val_int, gal_id, msg_end, rank, master_only, &
                     info, set_info, ndec)
    ! Prints a message
    ! Possible outputs (depending on the presence of optional arguments):
    !    rank: $msg
    !    rank: Galaxy $gal_id - $msg
    !    rank: Galaxy $gal_id - $msg $val
    !    rank: Galaxy $gal_id - $msg $val $msg_end
    !  Or, in the case of a serial run
    !    msg
    !    Galaxy $gal_id - $msg
    !    Galaxy $gal_id - $msg $val
    !    Galaxy $gal_id - $msg $val $msg_end
    !  If the rank variable is present, updates the module variable
    character(len=*), intent(in) :: msg
    double precision, optional, intent(in) :: val
    integer, optional, intent(in) :: val_int, set_info, info
    integer, optional, intent(in) :: gal_id, rank, ndec
    character(len=*), optional, intent(in) :: msg_end
    logical, optional, intent(in) :: master_only
    character(len=5000) :: p_rank, p_msg_end, p_val, p_gal
    integer :: verbosity
    verbosity = 1
    p_rank=''; p_val=''; p_gal=''; p_msg_end=''

    ! Stores module variables if they are present
    if (present(rank))  MPI_rank = rank
    if (present(set_info)) verbosity_setting = set_info

    ! Exits if the verbosity level (info) is greater the verbosity setting
    if (present(info)) verbosity = info
    if (verbosity >= verbosity_setting) return

    ! Writes $rank: unless single precessor or master_only=True
    if (MPI_rank>=0)  p_rank = str(MPI_rank)//':  '
    if (present(master_only)) then
      if (master_only .and. (MPI_rank > 0)) then
        return
      else
        p_rank =''
      endif
    endif

    ! Adds numerical informtion to the output
    if (present(val_int)) then
      p_val = str(val_int)
    else
      if (present(ndec)) then
        if (present(val))  p_val = str(val,ndec)
      else
        if (present(val))  p_val = str(val)
      endif
    endif

    ! Adds a possible suffix
    if (present(msg_end))  p_msg_end = msg_end

    ! Prints, accounting for the presence of gal_id
    if (present(gal_id))  then
      p_gal = 'Galaxy '//str(gal_id)//' -'
      print *, trim(p_rank),'  ',trim(p_gal),' ',msg,' ',trim(p_val),' ',trim(p_msg_end)
    else
      print *, trim(p_rank),'  ',msg,' ',trim(p_val),' ',trim(p_msg_end)
    endif
  end subroutine message

  pure integer function str_int_len(i) result(sz)
    ! Returns the length of the string representation of 'i'
    integer, intent(in) :: i
    integer, parameter :: MAX_STR = 100
    character(MAX_STR) :: s
    ! If 's' is too short (MAX_STR too small), Fortan will abort with:
    ! "Fortran runtime error: End of record"
    write(s, '(i0)') i
    sz = len_trim(s)
  end function str_int_len

  pure function str_int(i) result(s)
    ! Converts integer "i" to string
    integer, intent(in) :: i
    character(len=str_int_len(i)) :: s
    write(s, '(i0)') i
  end function str_int

  pure integer function str_dp_len(r, fmt) result(sz)
    ! Returns the length of the string representation of 'i'
    double precision, intent(in) :: r
    character(len=*), intent(in) :: fmt
    integer, parameter :: MAX_STR = 100
    character(MAX_STR) :: s
    ! If 's' is too short (MAX_STR too small), Fortan will abort with:
    ! "Fortran runtime error: End of record"
    write(s, fmt) r
    sz = len_trim(s)
  end function str_dp_len

  pure function str_dp(r) result(s)
    ! Converts the real number "r" to string with 7 decimal digits.
    double precision, intent(in) :: r
    character(len=*), parameter :: fmt="(f0.6)"
    character(len=str_dp_len(r, fmt)) :: s
    write(s, fmt) r
  end function str_dp

  pure function str_dp_n(r, n) result(s)
  ! Converts the real number "r" to string with 'n' decimal digits.
  double precision, intent(in) :: r
  integer, intent(in) :: n
  character(len=str_dp_len(r, "(f0." // str_int(n) // ")")) :: s
  write(s, "(f0." // str_int(n) // ")") r
  end function str_dp_n

  function str2int(str) result(i)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer :: i
    integer :: stat

    read(str,*,iostat=stat)  i
    if (stat /= 0) then
      print *, 'Error: conversion of string ',str,' failed!'
      stop
    endif
  end function str2int

end module messages
