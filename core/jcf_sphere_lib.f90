!------------------------------------------------------------------------------------
!> This module difines API for interpolation coefficient calculation,
!! especially for area/length calcuration on sperical cordinate.
!!
!! Calcuration of great circle related are based on
!! the Appendix A of Lauritzen et al.(2008)
!! > LAURITZEN, Peter H.; NAIR, Ramachandran D. Monotone and
!! > conservative cascade mapping between spherical grids (CaRS):
!! > regular latitude-longitude and cubed-sphere grids. Monthly Weather
!! > Review, 2008, 136.4: 1416-1432.
!!
!! \note
!! All latitude/longitude in argument list are in degree.
!!
!! \note
!! All subroutines except for init_sphere_lib are defined as public,
!! but should be used in library only.
!!
!! \author Takashi ARAKAWA <arakawa@rist.jp>
!------------------------------------------------------------------------------------

module jcf_sphere_lib
implicit none

  public :: cal_inverse_matrix
  public :: latlon2xyz
  public :: xyz2latlon

  public :: init_sphere_lib     !< initialize sphere_lib module
  public :: cal_great_circle_area ! real(kind=8) function (num_of_point, lat, lon)
  public :: get_length ! real(kind=8) function (lat1, lon1, lat2, lon2)
  public :: is_same_point  ! logical function (lat1, lon1, lat2, lon2)
  public :: is_on_line ! logical function(lat1, lon1, lat2, lon2, latc, lonc)
  public :: is_same_line ! logical function (lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, cross_lat, cross_lon)
  public :: is_cross_line ! logical function (lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, cross_lat, cross_lon)
  public :: cal_great_circle_center_rect !< calculate center of great circle
  public :: is_in_latlon                 !
  public :: is_inner                     !

  private :: cal_great_circle_triangle_area
  private :: cal_great_circle_dist
  private :: harversin


  !> genelic subroutine to  calculate inverse matrix
  interface cal_inverse_matrix
    module procedure cal_inverse_matrix_double, cal_inverse_matrix_quad
  end interface

  !> genelic subroutine to convert (lat,lon) to cartesian coordinate (X,Y,Z)
  interface latlon2xyz
    module procedure latlon2xyz_double, latlon2xyz_double_quad, latlon2xyz_quad
  end interface

  !> genelic subroutine to convert cartesian coordinate (X,Y,Z) to (lat,lon)
  interface xyz2latlon
    module procedure xyz2latlon_double, xyz2latlon_quad_double, xyz2latlon_quad
  end interface

  real(kind=8),  public, save :: PI !< the circular constant \f$\pi\f$
  real(kind=16), public, save :: PIQ !< the circular constant \f$\pi\f$ in quadruple precision.
  real(kind=8),  public, save :: D2R !< coefficient for degree to radian.

  logical,private :: is_initialized = .false.

contains

!===================================================================================
!> initialize this library
subroutine init_sphere_lib()
  implicit none

  PI = atan(1.d0)*4.d0
  PIQ= atan(1.q0)*4.q0
  D2R =PI/180.d0

  is_initialized = .true.

end subroutine init_sphere_lib

!===================================================================================
!> cal great circle triangle area.
!!
!! cite Lauritzen 2008 appendix A
real(kind=8) function cal_great_circle_area(num_of_point, lat, lon)
  implicit none
  integer, intent(IN) :: num_of_point !< number of point
  real(kind=8), intent(IN) :: lat(:)  !< array of latitude
  real(kind=8), intent(IN) :: lon(:)  !< array of longitude

  real(kind=8) :: area
  integer :: i

  area = 0.d0

  do i = 1, num_of_point-2
    area = area + cal_great_circle_triangle_area(lat(1), lon(1), lat(i+1), lon(i+1), lat(i+2), lon(i+2))
  end do

  cal_great_circle_area = area

end function cal_great_circle_area

!===================================================================================
!> cal great circle trinangle area
!!
!! cite Lauritzen 2008 appendix A
real(kind=8) function cal_great_circle_triangle_area(lat1, lon1, lat2, lon2, lat3, lon3)
  real(kind=8), intent(IN) :: lat1, lon1, lat2, lon2, lat3, lon3
  real(kind=8) :: a, b, c, s
  real(kind=8) :: temp1

  a = cal_great_circle_dist(lat1*D2R, lon1*D2R, lat2*D2R, lon2*D2R)
  b = cal_great_circle_dist(lat1*D2R, lon1*D2R, lat3*D2R, lon3*D2R)
  c = cal_great_circle_dist(lat3*D2R, lon3*D2R, lat2*D2R, lon2*D2R)
  
  s = 0.5d0*(a+b+c)

  temp1 =  max(tan(0.5d0*s)*tan((s-a)*0.5d0)*tan((s-b)*0.5d0)*tan((s-c)*0.5d0), 0.d0)
  cal_great_circle_triangle_area = 4.d0*atan(sqrt(temp1)) ! tan(0.5d0*s)*tan((s-a)*0.5d0)*tan((s-b)*0.5d0)*tan((s-c)*0.5d0)))

end function cal_great_circle_triangle_area

!===================================================================================
real(kind=8) function cal_great_circle_dist(lat1, lon1, lat2, lon2)
  implicit none
  real(kind=8), intent(IN) :: lat1, lon1, lat2, lon2

  cal_great_circle_dist = 2.d0*asin(sqrt(harversin(lat2-lat1)+cos(lat1)*cos(lat2)*harversin(lon2-lon1)))

end function cal_great_circle_dist

!===================================================================================

real(kind=8) function harversin(theta)
  implicit none
  real(kind=8), intent(IN) :: theta

  harversin = (1.d0-cos(theta))*0.5d0

end function harversin

!===================================================================================
!> cal distance of two points on the shere
real(kind=8) function get_length(lat1, lon1, lat2, lon2)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1
  real(kind=8), intent(IN) :: lon1 !< longitude of point1
  real(kind=8), intent(IN) :: lat2 !< latitude of point2
  real(kind=8), intent(IN) :: lon2 !< longitude of point2

  real(kind=8) :: x1, y1, z1, x2, y2, z2
  real(kind=8) :: inner_product
  call latlon2xyz(lat1, lon1, x1, y1, z1)
  call latlon2xyz(lat2, lon2, x2, y2, z2)

  inner_product = x1*x2+y1*y2+z1*z2
  if (inner_product >= 1.0d0) then
    get_length = 0.d0
  else
    get_length = acos(x1*x2+y1*y2+z1*z2)*0.5d0/PI
  end if

end function get_length

!===================================================================================
!> return .true. if the point is on the line defined by two points
logical function is_on_line(lat1, lon1, lat2, lon2, latc, lonc)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1       
  real(kind=8), intent(IN) :: lon1 !< longitude of point1      
  real(kind=8), intent(IN) :: lat2 !< latitude of point2       
  real(kind=8), intent(IN) :: lon2 !< longitude of point2      
  real(kind=8), intent(IN) :: latc !< latitude of target point 
  real(kind=8), intent(IN) :: lonc !< longitude of target point

  real(kind=16) :: x1, y1, z1, x2, y2, z2, xc, yc, zc
  real(kind=16) :: theta1, theta2, theta3

  call latlon2xyz(lat1, lon1, x1, y1, z1)
  call latlon2xyz(lat2, lon2, x2, y2, z2)
  call latlon2xyz(latc, lonc, xc, yc, zc)

  theta1 = acos(x1*x2+y1*y2+z1*z2)
  theta2 = acos(x1*xc+y1*yc+z1*zc)
  theta3 = acos(x2*xc+y2*yc+z2*zc)
  !write(0,*) "is_on_line ", theta1, theta2, theta3, theta1-theta2-theta3
  is_on_line = (theta1 >= theta2+theta3-1.d-12)

end function is_on_line

!===================================================================================
!> return .true. if two points are on the same spot.
logical function is_same_point(lat1, lon1, lat2, lon2)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1 
  real(kind=8), intent(IN) :: lon1 !< longitude of point1
  real(kind=8), intent(IN) :: lat2 !< latitude of point2 
  real(kind=8), intent(IN) :: lon2 !< longitude of point2

  real(kind=8) :: mlon1, mlon2
  real(kind=8) :: distance

  !write(0,*) "is_same_point1 ", lat1, lon1, lat2, lon2

  mlon1 = mod(lon1, 360.d0)
  mlon2 = mod(lon2, 360.d0)

  if ((lat2==lat1).and.(mlon2==mlon1)) then
    is_same_point = .true.
    return
  end if

  !write(0,*) "is_same_point2 ", lat2-lat1, lon2-lon1

  if ((abs(lat2-lat1) <= 1.0D-8).and.(abs(mlon2-mlon1) <= 1.0D-8)) then
    is_same_point = .true.
    return
  end if
  
  distance = get_length(lat1, mlon1, lat2, mlon2)

  !write(0,*) "is_same_point3 ", distance

  is_same_point = (distance < 1.0D-15)

end function is_same_point

!===================================================================================
!> return .true. if two lines (lon1, lat1)-(lon2, lat2) and
!! (lon3, lat3)-(lon4, lat4) are the same.
logical function is_same_line(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1 
  real(kind=8), intent(IN) :: lon1 !< longitude of point1
  real(kind=8), intent(IN) :: lat2 !< latitude of point2 
  real(kind=8), intent(IN) :: lon2 !< longitude of point2
  real(kind=8), intent(IN) :: lat3 !< latitude of point3 
  real(kind=8), intent(IN) :: lon3 !< longitude of point3
  real(kind=8), intent(IN) :: lat4 !< latitude of point4 
  real(kind=8), intent(IN) :: lon4 !< longitude of point4

  real(kind=16) :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd
  real(kind=16) :: matA(3,3), matI(3,3), datA

  call latlon2xyz(lat1, lon1, xa, ya, za)
  call latlon2xyz(lat2, lon2, xb, yb, zb)
  call latlon2xyz(lat3, lon3, xc, yc, zc)
  call latlon2xyz(lat4, lon4, xd, yd, zd)
  matA(1,1) = xb-xa ; matA(1,2) = xc-xd ; matA(1,3) = xd
  matA(2,1) = yb-ya ; matA(2,2) = yc-yd ; matA(2,3) = yd
  matA(3,1) = zb-za ; matA(3,2) = zc-zd ; matA(3,3) = zd

  call cal_inverse_matrix(matA, matI, datA)

  if (abs(datA) < 1.0q-32) then
    is_same_line = .true.
  else
    is_same_line = .false.
  end if

end function is_same_line

!===================================================================================
!> return .true. if two lines (lon1, lat1)-(lon2, lat2) and
!! (lon3,lat3)-(lon4, lat4) are crossing.
logical function is_cross_line(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, cross_lat, cross_lon)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1 
  real(kind=8), intent(IN) :: lon1 !< longitude of point1
  real(kind=8), intent(IN) :: lat2 !< latitude of point2 
  real(kind=8), intent(IN) :: lon2 !< longitude of point2
  real(kind=8), intent(IN) :: lat3 !< latitude of point3 
  real(kind=8), intent(IN) :: lon3 !< longitude of point3
  real(kind=8), intent(IN) :: lat4 !< latitude of point4 
  real(kind=8), intent(IN) :: lon4 !< longitude of point4
  real(kind=8), intent(OUT) :: cross_lat !< latitude of cross point
  real(kind=8), intent(OUT) :: cross_lon !< longitude of cross point 

  real(kind=8) :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd
  real(kind=8) :: cross_x, cross_y, cross_z
  real(kind=16) :: matA(3,3), matI(3,3), datA
  real(kind=16) :: alpha, beta, t, length

  call latlon2xyz(lat1, lon1, xa, ya, za)
  call latlon2xyz(lat2, lon2, xb, yb, zb)
  call latlon2xyz(lat3, lon3, xc, yc, zc)
  call latlon2xyz(lat4, lon4, xd, yd, zd)
  matA(1,1) = xb-xa ; matA(1,2) = xc-xd ; matA(1,3) = xd
  matA(2,1) = yb-ya ; matA(2,2) = yc-yd ; matA(2,3) = yd
  matA(3,1) = zb-za ; matA(3,2) = zc-zd ; matA(3,3) = zd

  call cal_inverse_matrix(matA, matI, datA)

  alpha = (xb*matI(1,1)+yb*matI(1,2)+zb*matI(1,3))/datA
  t = xb*matI(3,1)+yb*matI(3,2)+zb*matI(3,3)
  beta = (xb*matI(2,1)+yb*matI(2,2)+zb*matI(2,3))/t

  is_cross_line = .false.
  cross_lat = 999.d0
  cross_lon = 999.d0

  if ((1.0d0+1.d-12 >= dble(alpha)).and.(dble(alpha) >= -1.d-45)) then
    if ((1.0d0+1.d-12 >= dble(beta)).and.(dble(beta) >= -1.d-45)) then
      is_cross_line = .true.
      cross_x = alpha*xa+(1.d0-alpha)*xb
      cross_y = alpha*ya+(1.d0-alpha)*yb
      cross_z = alpha*za+(1.d0-alpha)*zb
      length = sqrt(cross_x*cross_x+cross_y*cross_y+cross_z*cross_z)
      cross_x = cross_x/length
      cross_y = cross_y/length
      cross_z = cross_z/length
      call xyz2latlon(cross_x, cross_y, cross_z, cross_lat, cross_lon)
    end if
  end if

  !write(0,*) "is_cross_line ", lon1, lat1
  !write(0,*) "is_cross_line ", lon2, lat2
  !write(0,*) "is_cross_line ", lon3, lat3
  !write(0,*) "is_cross_line ", lon4, lat4
  !write(0,*) "is_cross_line ", is_cross_line
  !write(0,*) "is_cross_line ", datA, alpha, beta
  !write(0,*) "is_cross_line ", cross_lon, cross_lat
  !write(0,*)

end function is_cross_line

!===================================================================================
!===================================================================================
!> @brief 
!> calculate center of rectangle on great circle
subroutine cal_great_circle_center_rect(&
     & lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, clat, clon)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1 
  real(kind=8), intent(IN) :: lon1 !< longitude of point1
  real(kind=8), intent(IN) :: lat2 !< latitude of point2 
  real(kind=8), intent(IN) :: lon2 !< longitude of point2
  real(kind=8), intent(IN) :: lat3 !< latitude of point3 
  real(kind=8), intent(IN) :: lon3 !< longitude of point3
  real(kind=8), intent(IN) :: lat4 !< latitude of point4 
  real(kind=8), intent(IN) :: lon4 !< longitude of point4
  real(kind=8), intent(OUT) :: clon !< latitude of the center point
  real(kind=8), intent(OUT) :: clat !< longitude of the center point
  real(kind=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
  real(kind=8) :: cx, cy, cz

  if (.not.is_initialized) then
    write(0,*) "init_sphere_lib is not called"
    stop
  end if

  call latlon2xyz(lat1, lon1, x1, y1, z1)
  call latlon2xyz(lat2, lon2, x2, y2, z2)
  call latlon2xyz(lat3, lon3, x3, y3, z3)
  call latlon2xyz(lat4, lon4, x4, y4, z4)

  cx = (x1+x2+x3+x4)*0.25d0
  cy = (y1+y2+y3+y4)*0.25d0
  cz = (z1+z2+z3+z4)*0.25d0

  call xyz2latlon(cx, cy, cz, clat, clon)

end subroutine cal_great_circle_center_rect

!!$ THIS SUBROUTINE IS OBSOLETE !!
!!$
!===================================================================================
!!$!> calculate maximum latitude and longitude of the line on a great circle.
!!$
!!$subroutine cal_great_circle_max_lat(lat1, lon1, lat2, lon2, mlat, mlon)
!!$  implicit none
!!$  real(kind=8), intent(IN) :: lat1 !< latitude of point1 
!!$  real(kind=8), intent(IN) :: lon1 !< longitude of point1
!!$  real(kind=8), intent(IN) :: lat2 !< latitude of point2 
!!$  real(kind=8), intent(IN) :: lon2 !< longitude of point2
!!$  real(kind=8), intent(OUT) :: mlat !< maximum latitude
!!$  real(kind=8), intent(OUT) :: mlon !< center longitude
!!$  real(kind=8) :: x1, y1, z1, x2, y2, z2, mx, my, mz
!!$  real(kind=8) :: a, b
!!$
!!$  call latlon2xyz(lat1, lon1, x1, y1, z1)
!!$  call latlon2xyz(lat2, lon2, x2, y2, z2)
!!$
!!$  a = (z2*y1-z1*y2)/(x1*y2-x2*y1)
!!$  b = (z2*x1-z1*x2)/(x2*y1-x1*y2)
!!$
!!$  mlon = atan2(a, b)
!!$  mlat = atan2(1.d0, sqrt(a*a*b*b))
!!$
!!$end subroutine cal_great_circle_max_lat

!===================================================================================
!> return .true. if given point(in lat-lon) is in the rectangle given pair of lat,lon.
logical function is_in_latlon(lon1, lat1, lon2, lat2, tlon, tlat)
  implicit none
  real(kind=8), intent(IN) :: lat1 !< latitude of point1 
  real(kind=8), intent(IN) :: lon1 !< longitude of point1
  real(kind=8), intent(IN) :: lat2 !< latitude of point2 
  real(kind=8), intent(IN) :: lon2 !< longitude of point2
  real(kind=8), intent(IN) :: tlon !< latitude of tested point
  real(kind=8), intent(IN) :: tlat !< longitude of tested point

  real(kind=8) :: lon1m, lon2m, tlonm

  if (tlon==360.d0) then
    tlonm = tlon
  else
     tlonm = mod(tlon+360.d0, 360.d0)
  end if

  if (lon2==360.d0) then
    lon2m = 360.d0
  else
    lon2m = mod(lon2+360.d0, 360.d0)
  end if

  lon1m = mod(lon1+360.d0, 360.d0)
 
  is_in_latlon = .false.

  if (tlat < min(lat1, lat2)) return
  if (tlat > max(lat1, lat2)) return
  if (tlonm < min(lon1m, lon2m)) return
  if (tlonm > max(lon1m, lon2m)) return

  is_in_latlon = .true.
  
end function is_in_latlon

!===================================================================================
!> return .true. if given point is in the triangle formed on by three points.
logical function is_inner(lon1, lat1, lon2, lat2, lon3, lat3, tlon, tlat)
  implicit none
  real(kind=8), intent(IN) :: lon1 !< longitude of triangle point 1
  real(kind=8), intent(IN) :: lat1 !< latitude of triangle point 1 
  real(kind=8), intent(IN) :: lon2 !< longitude of triangle point 2
  real(kind=8), intent(IN) :: lat2 !< latitude of triangle point 2
  real(kind=8), intent(IN) :: lon3 !< longitude of triangle point 3
  real(kind=8), intent(IN) :: lat3 !< latitude of triangle point 3
  real(kind=8), intent(IN) :: tlon !< longitude of tested point
  real(kind=8), intent(IN) :: tlat !< latitude of tested point

  real(kind=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
  real(kind=8) :: tx, ty, tz
  real(kind=8) :: detA
  real(kind=8) :: alpha, beta, gamma

  is_inner = .false.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This min, max check is incorrect.
  ! Because, triangle side is great circle.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !if (tlat < min(lat1, lat2, lat3)) return
  !if (tlat > max(lat1, lat2, lat3)) return

  !write(0,*) "is_inner 1 ", PI

  call latlon2xyz(lat1, lon1, x1, y1, z1)

  !write(0,*) "is_inner 2 ", PI

  call latlon2xyz(lat2, lon2, x2, y2, z2)

  call latlon2xyz(lat3, lon3, x3, y3, z3)

  call latlon2xyz(tlat, tlon, tx, ty, tz)

  detA = x1*y2*z3+x3*y1*z2+x2*y3*z1-x3*y2*z1-x1*y3*z2-x2*y1*z3
  detA = sign(1.d0,detA)

!!$  alpha = ((y2*z3-y3*z2)*tx+(z2*x3-z3*x2)*ty+(x2*y3-x3*y2)*tz)/detA
!!$  beta = ((y3*z1-y1*z3)*tx+(z3*x1-z1*x3)*ty+(x3*y1-x1*y3)*tz)/detA
!!$  gamma = ((y1*z2-y2*z1)*tx+(z1*x2-z2*x1)*ty+(x1*y2-x2*y1)*tz)/detA

  alpha = sign(1.d0,((y2*z3-y3*z2)*tx+(z2*x3-z3*x2)*ty+(x2*y3-x3*y2)*tz))
  if (alpha*detA<0.d0) return

  beta  = sign(1.d0,((y3*z1-y1*z3)*tx+(z3*x1-z1*x3)*ty+(x3*y1-x1*y3)*tz))
  if (beta*detA<0.d0) return

  gamma = sign(1.d0,((y1*z2-y2*z1)*tx+(z1*x2-z2*x1)*ty+(x1*y2-x2*y1)*tz))
  if (gamma*detA<0.d0) return

  is_inner = .true.
  
end function is_inner

!===================================================================================
!> calc inverse matrix (for double precision)
!!
!! Matrix A,I are defined as follows;
!!
!!      a11 a12 a13
!!      a21 a22 a23
!!      a31 a32 a33
subroutine cal_inverse_matrix_double(A, I, datA)
  implicit none
  real(kind=8), intent(IN) :: A(3, 3)  !< 3x3 matrix
  real(kind=8), intent(OUT) :: I(3, 3) !< inverse of matrix A
  real(kind=8), intent(OUT) :: datA    !< a determinant of A

  datA = A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(3,1)*A(1,2)*A(2,3) &
        -A(1,1)*A(3,2)*A(2,3)-A(3,1)*A(2,2)*A(1,3)-A(2,1)*A(1,2)*A(3,3)

  I(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
  I(1,2) = A(1,3)*A(3,2)-A(1,2)*A(3,3)
  I(1,3) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
  I(2,1) = A(2,3)*A(3,1)-A(2,1)*A(3,3)
  I(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
  I(2,3) = A(1,3)*A(2,1)-A(1,1)*A(2,3)
  I(3,1) = A(2,1)*A(3,2)-A(2,2)*A(3,1)
  I(3,2) = A(1,2)*A(3,1)-A(1,1)*A(3,2)
  I(3,3) = A(1,1)*A(2,2)-A(1,2)*A(2,1)

end subroutine cal_inverse_matrix_double

!===================================================================================
!> calc inverse matrix (for quadruple precision)
!!
!! Matrix A,I are defined as follows;
!!
!!      a11 a12 a13
!!      a21 a22 a23
!!      a31 a32 a33
subroutine cal_inverse_matrix_quad(A, I, datA)
  implicit none
  real(kind=16), intent(IN) :: A(3, 3)  !< 3x3 matrix
  real(kind=16), intent(OUT) :: I(3, 3) !< inverse of matrix A
  real(kind=16), intent(OUT) :: datA    !< a determinant of A

  datA = A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(3,1)*A(1,2)*A(2,3) &
        -A(1,1)*A(3,2)*A(2,3)-A(3,1)*A(2,2)*A(1,3)-A(2,1)*A(1,2)*A(3,3)

  I(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
  I(1,2) = A(1,3)*A(3,2)-A(1,2)*A(3,3)
  I(1,3) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
  I(2,1) = A(2,3)*A(3,1)-A(2,1)*A(3,3)
  I(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
  I(2,3) = A(1,3)*A(2,1)-A(1,1)*A(2,3)
  I(3,1) = A(2,1)*A(3,2)-A(2,2)*A(3,1)
  I(3,2) = A(1,2)*A(3,1)-A(1,1)*A(3,2)
  I(3,3) = A(1,1)*A(2,2)-A(1,2)*A(2,1)

end subroutine cal_inverse_matrix_quad

!===================================================================================
!> converte (lat,lon) to (x, y, z) in the cartesian coordinate.
subroutine latlon2xyz_double(lat, lon, x, y, z)
  implicit none
  real(kind=8), intent(IN) :: lat !< latitude [deg]
  real(kind=8), intent(IN) :: lon !< longitude [deg]
  real(kind=8), intent(OUT) :: x !< x in cartesian
  real(kind=8), intent(OUT) :: y !< y in cartesian
  real(kind=8), intent(OUT) :: z !< z in cartesian

!!$d  x = cos(lat/180.d0*PI)*cos(lon/180.d0*PI)
!!$d  y = cos(lat/180.d0*PI)*sin(lon/180.d0*PI)
!!$d  z = sin(lat/180.d0*PI)
  x = cos(lat*D2R)*cos(lon*D2R)
  y = cos(lat*D2R)*sin(lon*D2R)
  z = sin(lat*D2R)

end subroutine latlon2xyz_double

!===================================================================================
!> converte (lat,lon) to (x, y, z) in the cartesian coordinate.
!!
!! \note x,y,z are quadruple(kind=16).
subroutine latlon2xyz_double_quad(lat, lon, x, y, z)
  implicit none
  real(kind=8), intent(IN) :: lat !< latitude [deg]
  real(kind=8), intent(IN) :: lon !< longitude [deg]
  real(kind=16), intent(OUT) :: x !< x in cartesian
  real(kind=16), intent(OUT) :: y !< y in cartesian
  real(kind=16), intent(OUT) :: z !< z in cartesian

  x = cos(lat/180.q0*piq)*cos(lon/180.q0*piq)
  y = cos(lat/180.q0*piq)*sin(lon/180.q0*piq)
  z = sin(lat/180.q0*piq)

end subroutine latlon2xyz_double_quad

!===================================================================================
!> converte (lat,lon) to (x, y, z) in the cartesian coordinate.
!!
!! \note all arguments are quadruple(kind=16).
subroutine latlon2xyz_quad(lat, lon, x, y, z)
  implicit none
  real(kind=16), intent(IN) :: lat !< latitude [deg]
  real(kind=16), intent(IN) :: lon !< longitude [deg]
  real(kind=16), intent(OUT) :: x  !< x in cartesian
  real(kind=16), intent(OUT) :: y  !< y in cartesian
  real(kind=16), intent(OUT) :: z  !< z in cartesian

  x = cos(lat/180.q0*piq)*cos(lon/180.q0*piq)
  y = cos(lat/180.q0*piq)*sin(lon/180.q0*piq)
  z = sin(lat/180.q0*piq)

end subroutine latlon2xyz_quad

!===================================================================================
!> converte (x, y, z) in the cartesian coordinate to (lat,lon).
subroutine xyz2latlon_double(x, y, z, lat, lon)
  implicit none
  real(kind=8), intent(IN) :: x    !< x in cartesian
  real(kind=8), intent(IN) :: y    !< y in cartesian
  real(kind=8), intent(IN) :: z    !< z in cartesian
  real(kind=8), intent(OUT) :: lat !< latitude [deg]
  real(kind=8), intent(OUT) :: lon !< longitude [deg]

  real(kind=8),parameter :: EPSLN=1.d-12
!
  if ( abs(x) < EPSLN .and. abs(y) < EPSLN ) then
    lon = 0.d0
    if ( z > 0.d0 ) then
      lat = 90.d0
    else
      lat = -90.d0
    end if
    return
  end if


  lat = atan2(z, sqrt(x*x+y*y))
  !lon = acos(x/cos(lat))
  lon = acos(x/sqrt(x*x+y*y))

  lat = lat*180.d0/pi
  lon = lon*180.d0/pi

  if (y<0) lon = 360.d0-lon

  if ( abs(x) < EPSLN ) then
    write(33,*)'xyz2latlon_double:',x,y,z,lon,lat
  end if

end subroutine xyz2latlon_double

!===================================================================================
!> converte (x, y, z) in the cartesian coordinate to (lat,lon).
!!
!> \note x,y,z are quadruple (kind=16)
subroutine xyz2latlon_quad_double(x, y, z, lat, lon)
  implicit none
  real(kind=16), intent(IN) :: x   !< x in cartesian
  real(kind=16), intent(IN) :: y   !< y in cartesian
  real(kind=16), intent(IN) :: z   !< z in cartesian
  real(kind=8), intent(OUT) :: lat !< latitude [deg]
  real(kind=8), intent(OUT) :: lon !< longitude [deg]

  lat = atan2(z, sqrt(x*x+y*y))
  !lon = acos(x/cos(lat))
  lon = acos(x/sqrt(x*x+y*y))

  lat = lat*180.d0/pi
  lon = lon*180.d0/pi

  if (y<0) lon = 360.d0-lon

end subroutine xyz2latlon_quad_double

!===================================================================================
!> converte (x, y, z) in the cartesian coordinate to (lat,lon).
!!
!> \note all arguments are quadruple (kind=16)
subroutine xyz2latlon_quad(x, y, z, lat, lon)
  implicit none
  real(kind=16), intent(IN) :: x   !< x in cartesian
  real(kind=16), intent(IN) :: y   !< y in cartesian
  real(kind=16), intent(IN) :: z   !< z in cartesian
  real(kind=16), intent(OUT) :: lat !< latitude [deg]
  real(kind=16), intent(OUT) :: lon !< longitude [deg]

  lat = atan2(z, sqrt(x*x+y*y))
  !lon = acos(x/cos(lat))
  lon = acos(x/sqrt(x*x+y*y))

  lat = lat*180.q0/piq
  lon = lon*180.q0/piq

  if (y<0) lon = 360.q0-lon

end subroutine xyz2latlon_quad

!===================================================================================

end module jcf_sphere_lib
