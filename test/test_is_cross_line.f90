program test_is_cross_line
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8

  real(DP_K) :: lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4
  real(DP_K) :: lon,lat
  logical :: result

  call init_sphere_lib()

  do
    write(*,*)'input lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4'
    read(*,*,err=999) lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4
    write(*,*)lon1,lat1,lon2,lat2
    write(*,*)lon3,lat3,lon4,lat4
    result=is_cross_line(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4,lat,lon)
    write(*,*)'is_cross_line,lon,lat:',result,lon,lat
  end do

999 call exit(0)

end program test_is_cross_line

