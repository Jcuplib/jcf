program test_xyz2latlon
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8
  integer,parameter :: qp_k=16

  real(dp_k) :: dlon, dlat
  real(qp_k) :: qlon, qlat

  real(dp_k) :: xd, yd, zd
  real(qp_k) :: xq, yq, zq


  call init_sphere_lib()

  do 
    write(*,*)'input xq,yq,zq (Q for quit):'
    read(*,*,err=999) xq,yq,zq
    xd = real( xq,dp_k )
    yd = real( yq,dp_k )
    zd = real( zq,dp_k )

    write(*,'(1x,3(F21.15,5x))') xd,yd,zd
    write(*,'(1x,3(F40.34,6x))') xq,yq,zq

    call xyz2latlon( xd,yd,zd,dlat,dlon)
    write(*,*)'xyz2latlon_double:'
    write(*,'(2(F21.15,5x))') dlon,dlat

    call xyz2latlon( xq,yq,zq,dlat,dlon)
    write(*,*)'xyz2latlon_quad_double:'
    write(*,'(2(F21.15,5x))') dlon,dlat

    call xyz2latlon( xq,yq,zq,qlat,qlon)
    write(*,*)'xyz2latlon_quad:'
    write(*,'(2(F40.34,6x))') qlon,qlat

    write(*,*)
  end do

999 write(*,*)'done.'
  call exit(0)


end program test_xyz2latlon

