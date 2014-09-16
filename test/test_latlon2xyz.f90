program test_latlon2xyz
  use jcf_sphere_lib
  implicit none 

  integer,parameter :: dp_k=8
  integer,parameter :: qp_k=16

  real(dp_k) :: dlon, dlat
  real(qp_k) :: qlon, qlat

  real(dp_k) :: xd, yd, zd
  real(qp_k) :: xq, yq, zq

  integer :: ios=0

  call init_sphere_lib()

  do 
    write(*,*)'input dlon,dlat (Q for quit):'
    read(*,*,iostat=ios) dlon, dlat
    if ( ios .ne. 0 ) exit
    qlon = real( dlon,qp_k )
    qlat = real( dlat,qp_k )

    write(*,'(1x,2(F21.15,5x))') dlon, dlat
    write(*,'(1x,2(F40.34,6x))') qlon, qlat

    call latlon2xyz( dlat,dlon,xd,yd,zd )
    write(*,*)'latlon2xyz_double:'
    write(*,'(3(F21.15,5x))') xd,yd,zd

    call latlon2xyz( dlat,dlon,xq,yq,zq )
    write(*,*)'latlon2xyz_double_quad:'
    write(*,'(3(F40.34,6x))') xq,yq,zq

    call latlon2xyz( qlat,qlon,xq,yq,zq )
    write(*,*)'latlon2xyz_quad:'
    write(*,'(3(F40.34,6x))') xq,yq,zq

    write(*,*)
  end do

  write(*,*) 'done.'
  call exit(0)

end program test_latlon2xyz
