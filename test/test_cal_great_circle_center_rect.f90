program test_cal_great_circle_center_rect
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8

  real(DP_K) :: lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4
  real(DP_K) :: lon,lat

  call init_sphere_lib()

  do
    write(*,'(A)')'input lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4:'
    read(*,*,err=999) lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4
    write(*,'(1x,4(F21.15,5x))')lon1,lat1,lon2,lat2
    write(*,'(1x,4(F21.15,5x))')lon3,lat3,lon4,lat4
    call cal_great_circle_center_rect( &
         & lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4, &
         & lat, lon )
    write(*,'(A,2(F21.15,5x))')'center point:',lon,lat

    write(*,*)
  end do

999 write(*,*)'done.'
  call exit(0)

end program test_cal_great_circle_center_rect

