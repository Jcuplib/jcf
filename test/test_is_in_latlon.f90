program test_is_in_latlon
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8

  real(DP_K) :: lon1,lat1,lon2,lat2
  real(DP_K) :: lon,lat

  integer :: ios=0

  call init_sphere_lib()

  do
    write(*,'(A)')'input lon1,lat1,lon2,lon,lat:'
    read(*,*,iostat=ios) lon1,lat1,lon2,lat2,lon,lat
    if ( ios .ne. 0 ) exit
    write(*,'(1x,4(F21.15,5x))')lon1,lat1,lon2,lat2
    write(*,'(1x,2(F21.15,5x))')lon,lat
    write(*,'(A,L2)')'is_in_latlon:',is_in_latlon(lon1,lat1,lon2,lat2,lon,lat)

    write(*,*)
  end do

  write(*,*) 'done.'
  call exit(0)

end program test_is_in_latlon

