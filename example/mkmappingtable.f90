program prg_mkmappingtable
  use nicamio_grid, only : init_nicamio, &
       nicamio_read_grid_info => read_grid_info, &
       nicamio_search_target_polygons => search_target_polygons, &
       cal_coefficient_nicam_to_io_area, &
       cal_coefficient_nicam_to_io_monte_carlo, &
       interpolate_nicam_to_io
  use mod_logfile, only : LOG_FID, open_logfile

  implicit none
  integer, parameter :: STR_LEN = 128
  character(len=STR_LEN) :: file_name = "./mappingtable.cnf"  
  character(len=32)      :: mapping_target
  integer :: coef_check_mode = 1
  character(len=STR_LEN) :: coef_cal_type = "AREA_RATIO"
  integer :: div_level = 5
  character(len=STR_LEN) :: mapping_table_nicam_to_io   = "mapping_table_nicam_to_io.txt"

  namelist / mappingtable / mapping_target, &
       coef_check_mode, &
       coef_cal_type, &
       div_level, &
       mapping_table_nicam_to_io


  integer :: ierror

  integer,parameter :: ctl_fid = 11


  open(CTL_FID,             &
       file=trim(file_name),   &
       form='formatted',    &
       status='old',        &
       iostat=ierror)
  if(ierror/=0) then
    write(*,*) 'Cannot open PARAMETER file!'
    stop
  end if
  rewind(ctl_fid)
  read(ctl_fid,nml=mappingtable,end=999)
999 close(ctl_fid)
  write(*,nml=mappingtable)

  call open_logfile("mk_mappingtable.log")

  write(LOG_FID, *) "Msg : Sub[prg_makemappingtable]/Mod[main]"
  write(LOG_FID, *) "----- read configure file"
  write(LOG_FID, *) "----- file_name : ", trim(file_name)
  write(LOG_FID, *) "------- configuration -------"
  write(LOG_FID, *) "mapping_target : ", trim(mapping_target)
  write(LOG_FID, '(A,I5)') " coef_check_mode : ", coef_check_mode
  write(LOG_FID, *) "coef_cal_type : ", trim(coef_cal_type)
  write(LOG_FID, '(A,I5)') " div_level : ", div_level
  write(LOG_FID, *) "mapping_table_nicam_to_io   : ", trim(mapping_table_nicam_to_io)
  write(LOG_FID, *)


  select case(trim(mapping_target))
  case ("NICAMIO")
    call init_nicamio(file_name)
    call nicamio_read_grid_info()
    call nicamio_search_target_polygons()

    select case ( trim(coef_cal_type) )
    case ( "AREA_RATIO" )
      call cal_coefficient_nicam_to_io_area(mapping_table_nicam_to_io, coef_check_mode) 
    case ( "MONTE_CARLO" )
      call cal_coefficient_nicam_to_io_monte_carlo(mapping_table_nicam_to_io, coef_check_mode, div_level)
    case default 
      write(0,*) "coef_cal_type error:"//trim(mapping_target)
      stop
    end select
  case default
    write(0,*) "mapping_target error:"//trim(mapping_target)
    stop
  end select

  write(LOG_FID, *) "Msg : Sub[prg_makemappingtable]/Mod[main]"
  write(LOG_FID, *) "----- mk_mappingtable finish"

  close(LOG_FID)

  stop

end program prg_mkmappingtable
