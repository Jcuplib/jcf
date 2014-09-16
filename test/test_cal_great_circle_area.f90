program test_cal_great_circle_area
  use jcf_sphere_lib
  implicit none 
  integer,parameter :: DP_K=8

  integer :: NP
  real(DP_K),allocatable  :: lon(:),lat(:)
  real(DP_K) :: area
  integer :: i

  integer :: ios=0

  call init_sphere_lib()

  do
    write(*,'(A)',advance='NO')'input NP (Q for quit):'
    read(*,*,iostat=ios)NP
    if ( ios .ne. 0 ) exit
    write(*,'(I2)')NP
    if ( NP<3 ) exit
    allocate( lon(NP),lat(NP) )
    write(*,'(A,I2,A)')'input (lon,lat) line by line(Q for quit):'
    do i=1,NP
      read(*,*,iostat=ios) lon(i),lat(i)
    end do
      write(*,'(2A25)')'lon(deg)','lat(deg)'
    do i=1,NP
      write(*,'(2F25.15)')lon(i),lat(i)
    end do
    area = cal_great_circle_area(np,lat,lon)
    write(*,'(A,F25.15)')'area:',area
    deallocate( lon,lat )
    write(*,*)
  end do

  write(*,*)'done.'
  call exit(0)

end program test_cal_great_circle_area

