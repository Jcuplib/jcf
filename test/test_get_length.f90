program test_get_length
  use jcf_sphere_lib
  implicit none 

  integer,parameter::DP_K=8

  real(DP_K) :: lon1, lat1, lon2, lat2

  call init_sphere_lib()

  do
    write(*,'(A)')'input lon1,lat1,lon2,lat2 (Q for quit):'
    read(*,*,err=999) lon1,lat1,lon2,lat2
    write(*,'(1x,2(F21.15,5x))')lon1,lat1
    write(*,'(1x,2(F21.15,5x))')lon2,lat2
    write(*,'(A,F21.15)')'length:',get_length( lat1,lon1,lat2,lon2)

    write(*,*)
  end do

999 write(*,*)'done.'
  call exit(0)

end program test_get_length

