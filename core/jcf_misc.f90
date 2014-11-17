module jcf_misc
  implicit none 

  public :: jcf_set_log_fid
  public :: jcf_set_cnf_fid
  public :: jcf_avail_fid

  integer,public :: jcf_log_fid = 731
  integer,public :: jcf_cnf_fid = 741


  integer :: LUN_MIN = 751

contains
  subroutine jcf_set_log_fid( fid )
    integer, intent(in) :: fid
    jcf_log_fid = fid
    return
  end subroutine jcf_set_log_fid


  subroutine jcf_set_cnf_fid( fid )
    integer, intent(in) :: fid
    jcf_cnf_fid = fid
    return
  end subroutine jcf_set_cnf_fid


  function jcf_avail_fid() result(res)
    integer :: res

    integer :: i
    logical :: opened

    i = LUN_MIN

    loop:do
      inquire(unit=i, opened = opened )
      if ( .not. opened ) exit loop
      i = i+1
    end do loop

    res = i
  end function jcf_avail_fid


end module jcf_misc
