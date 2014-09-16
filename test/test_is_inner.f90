program test_is_inner
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8

  real(DP_K) :: lon1,lat1,lon2,lat2,lon3,lat3
  real(DP_K) :: lon,lat

  integer :: ios=0

  call init_sphere_lib()

  do
    write(*,*)'input lon1,lat1,lon2,lat2,lon3,lat3,lon,lat:'
    read(*,*,iostat=ios) lon1,lat1,lon2,lat2,lon3,lat3,lon,lat
    if ( ios .ne. 0 ) exit
    write(*,'(1x,6(F21.15,5x))')lon1,lat1,lon2,lat2,lon3,lat3
    write(*,'(1x,2(F21.15,5x))')lon,lat
    write(*,'(A,L2)')'is_inner:',is_inner(lon1,lat1,lon2,lat2,lon3,lat3,lon,lat)

    write(*,*)
  end do

  write(*,*)'done.'
  call exit(0)

end program test_is_inner

