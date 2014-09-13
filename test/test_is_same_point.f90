program test_is_same_point
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8

  real(DP_K) :: lon1,lat1,lon2,lat2
  real(DP_K) :: lon,lat

  call init_sphere_lib()

  do
    write(*,'(A)',advance='NO')'input lon1,lat1,lon2,lat2:'
    read(*,*,err=999) lon1,lat1,lon2,lat2
    write(*,'(1x,4(F21.15,5x))')lon1,lat1,lon2,lat2
    write(*,*)'is_same_point:',is_same_point(lat1,lon1,lat2,lon2)
  end do

999 write(*,*)'done.'
  call exit(0)

end program test_is_same_point

