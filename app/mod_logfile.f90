module mod_logfile
  implicit none 
  
  public :: open_logfile

  integer, public, parameter :: LOG_FID = 22

contains
  subroutine open_logfile(logfile)
    character (len=*), intent(in) :: logfile

    logical :: flag = .false.

    inquire(unit=log_fid,opened=flag)
    if ( flag ) return

    open(UNIT = LOG_FID, FILE = logfile)

  end subroutine open_logfile


end module mod_logfile


