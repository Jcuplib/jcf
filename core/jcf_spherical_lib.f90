!> Module for spherical geometry, especially for spherical triangle and polygons.
!!
!! This module implements subroutines/functions used for interpolation
!! coefficient calculation, especially for area/length calcuration on
!! a spehere.
!!
!! Calcuration of great circle related are based on
!! the Appendix A of [Lauritzen et al.(2007)](http://journals.ametsoc.org/doi/abs/10.1175/2007MWR2181.1).
!! > LAURITZEN, Peter H.; NAIR, Ramachandran D. ''Monotone and
!! > conservative cascade mapping between spherical grids (CaRS):
!! > regular latitude-longitude and cubed-sphere grids''. Monthly Weather
!! > Review, 2007, 136.4: 1416-1432.
!!
!!
!! In this library, it is assumed that a point is on the shpere with
!! it radius is 1, its coordinates is represented in (lon,lat) and/or
!! cartesian coord (x,y,z). Subroutine xyz2lonlat() and lonlat2xyz()
!! convert them mutually.
!!
!! \note
!! All latitude/longitude in argument list are in degree.
!!
!! \note
!! Most of subroutines are defined as public,
!! but should be used in this library only, except for init_sphere_lib().
!!
!! \copyright
!!   The MIT License (MIT),  (c) 2013-2014 RIST.
!!   See the [LICENSE](@ref license) file.
!!   
!! 
!! \author Takashi ARAKAWA <arakawa@rist.jp>
!!
!! \version
!!   - ver.0.2.0 (2013/XX/XX T.Arakawa) Initial public version
!!   - ver.0.2.1 (2014/08/29 T.Inoue)   Refactored.
!!   .
!!
module jcf_spherical_lib
  implicit none

  private !! is default
  public :: dp_k, qp_k

  public :: init_spherical_lib
  public :: lonlat2xyz
  public :: xyz2lonlat

  public :: get_length_on_great_circle
  public :: get_area_of_spherical_polygon
  public :: get_area_of_spherical_triangle
  public :: get_intersection_point
  public :: get_center_of_spherical_rectangle

  public :: is_on_the_line
  public :: is_the_same_point
  public :: is_the_same_line
  public :: is_intersecting_lines
  public :: is_in_lonlat
  public :: is_inner

  !! Are below should be public ??
  public :: cal_inverse_matrix
  public :: haversine


  !> Genelic subroutine to convert (lat,lon) to cartesian coordinate (X,Y,Z)
  interface lonlat2xyz
    module procedure lonlat2xyz_double, lonlat2xyz_double_quad, lonlat2xyz_quad
  end interface lonlat2xyz


!!$  interface latlon2xyz
!!$    module procedure latlon2xyz_double, latlon2xyz_double_quad, latlon2xyz_quad
!!$  end interface

  !> Genelic subroutine to convert cartesian coordinate (X,Y,Z) to (lat,lon)
  interface xyz2lonlat
    module procedure xyz2lonlat_double, xyz2lonlat_quad_double, xyz2lonlat_quad
  end interface xyz2lonlat

  !!\internal
  !> Genelic subroutine to  calculate inverse matrix.
  interface cal_inverse_matrix
    module procedure cal_inverse_matrix_double, cal_inverse_matrix_quad
  end interface cal_inverse_matrix
  !!\endinternal

  integer,parameter :: SP_K = 4 !< kind parameter for single precision
  integer,parameter :: DP_K = 8 !< kind parameter for double precision
  integer,parameter :: QP_K = 16 !< kind parameter for quadruple precision


!!$#ifdef F2003
!!$  real(kind=DP_K), public,parameter :: PI  =atan(1.d0)*4.d0
!!$  real(kind=DP_K), public,parameter :: D2R =PI/180.d0
!!$  real(kind=DP_K), public,parameter :: R2D =180.d0/PI
!!$  real(kind=QP_K), public,parameter :: PIQ =atan(1.0_QP_K)*4.0_QP_K
!!$  real(kind=QP_K), public,parameter :: D2RQ =PIQ/180.0_QP_K
!!$  real(kind=QP_K), public,parameter :: R2DQ =180.0_QP_K/PIQ
!!$#else
  real(kind=DP_K), public,save :: PI   !! the circular constant \f$\pi\f$
  real(kind=DP_K), public,save :: D2R  !! coefficient for degree to radian.
  real(kind=DP_K), public,save :: R2D  !! coefficient for radian to degree.
  real(kind=QP_K), public,save :: PIQ !! the circular constant \f$\pi\f$ in quadruple precision.
  real(kind=QP_K), public,save :: D2RQ  !! coefficient for degree to radian in quadruple precision.
  real(kind=QP_K), public,save :: R2DQ  !! coefficient for radian to degree in quadruple precision.
!!$#endif



  logical,private :: is_initialized = .false.

contains

  !===================================================================================
  !> Initialize this module.
  !!
  !! Call this routine first.
  subroutine init_spherical_lib()
    implicit none

!!$#ifdef F2003
!!$#else
    PI = atan(1.d0)*4.d0
    D2R =PI/180.d0
    R2D =180.d0/PI
    PIQ= atan(1.0_QP_K)*4.0_QP_K
    D2RQ=PIQ/180.0_QP_K
    R2DQ=180.0_QP_K/PIQ
!!$#endif

!!$#ifdef DEBUG
!!$    write(0,*)'PI :',PI
!!$    write(0,*)'PIQ:',PIQ
!!$#endif

    is_initialized = .true.

  end subroutine init_spherical_lib


  !===================================================================================
  !> Convert (lon,lat) to (x, y, z) in the cartesian coordinate.
  !!
  subroutine lonlat2xyz_double(lon, lat, x, y, z)
    implicit none
    real(kind=DP_K), intent(IN) :: lon !< longitude [deg]
    real(kind=DP_K), intent(IN) :: lat !< latitude [deg]
    real(kind=DP_K), intent(OUT) :: x !< x coord in cartesian
    real(kind=DP_K), intent(OUT) :: y !< y coord in cartesian
    real(kind=DP_K), intent(OUT) :: z !< z coord in cartesian

    x = cos(lat*D2R)*cos(lon*D2R)
    y = cos(lat*D2R)*sin(lon*D2R)
    z = sin(lat*D2R)

  end subroutine lonlat2xyz_double

  !===================================================================================
  !> Convert (lon,lat) to (x, y, z) in the cartesian coordinate.
  !!
  !! \note x,y,z are quadruple(kind=QP_K).
  subroutine lonlat2xyz_double_quad(lon, lat, x, y, z)
    implicit none
    real(kind=DP_K), intent(IN) :: lon !< longitude [deg]
    real(kind=DP_K), intent(IN) :: lat !< latitude [deg]
    real(kind=QP_K), intent(OUT) :: x !< x coord in cartesian
    real(kind=QP_K), intent(OUT) :: y !< y coord in cartesian
    real(kind=QP_K), intent(OUT) :: z !< z coord in cartesian

    x = cos(lat*D2RQ)*cos(lon*D2RQ)
    y = cos(lat*D2RQ)*sin(lon*D2RQ)
    z = sin(lat*D2RQ)

  end subroutine lonlat2xyz_double_quad

  !===================================================================================
  !> Converte (lat,lon) to (x, y, z) in the cartesian coordinate.
  !!
  !! \note all arguments are quadruple(kind=QP_K).
  subroutine lonlat2xyz_quad(lon, lat, x, y, z)
    implicit none
    real(kind=QP_K), intent(IN) :: lon !< longitude [deg]
    real(kind=QP_K), intent(IN) :: lat !< latitude [deg]
    real(kind=QP_K), intent(OUT) :: x  !< x coord in cartesian
    real(kind=QP_K), intent(OUT) :: y  !< y coord in cartesian
    real(kind=QP_K), intent(OUT) :: z  !< z coord in cartesian

    x = cos(lat*D2RQ)*cos(lon*D2RQ)
    y = cos(lat*D2RQ)*sin(lon*D2RQ)
    z = sin(lat*D2RQ)

  end subroutine lonlat2xyz_quad



  !===================================================================================
  !> Convert (x, y, z) in the cartesian coordinate to (lon,lat).
  !!
  !OCL NOFP_RELAXED
  subroutine xyz2lonlat_double(x, y, z, lon, lat)
    implicit none
    real(kind=DP_K), intent(IN) :: x    !< x coord in cartesian
    real(kind=DP_K), intent(IN) :: y    !< y coord in cartesian
    real(kind=DP_K), intent(IN) :: z    !< z coord in cartesian
    real(kind=DP_K), intent(OUT) :: lon !< longitude [deg]
    real(kind=DP_K), intent(OUT) :: lat !< latitude [deg]

    real(kind=DP_K),parameter :: EPSLN=1.d-12

    if ( abs(x) < EPSLN .and. abs(y) < EPSLN ) then
      lon = 0.d0
      lat = sign( 90.d0, z )
      return
    end if

    lat = atan2(z, sqrt(x*x+y*y))
    lon = acos(x/sqrt(x*x+y*y))

    lat = lat*R2D
    lon = lon*R2D

    if (y<0) lon = 360.d0-lon

!!$#ifdef DEBUG
!!$    if ( abs(x) < EPSLN ) then
!!$      write(33,*)'xyz2lonlat_double:',x,y,z,lon,lat
!!$    end if
!!$#endif

  end subroutine xyz2lonlat_double



  !===================================================================================
  !> Convert (x, y, z) in the cartesian coordinate to (lat,lon).
  !!
  !> \note x,y,z are quadruple (kind=QP_K)
  subroutine xyz2lonlat_quad_double(x, y, z, lon, lat)
    implicit none
    real(kind=QP_K), intent(IN) :: x   !< x in cartesian
    real(kind=QP_K), intent(IN) :: y   !< y in cartesian
    real(kind=QP_K), intent(IN) :: z   !< z in cartesian
    real(kind=DP_K), intent(OUT) :: lon !< longitude [deg]
    real(kind=DP_K), intent(OUT) :: lat !< latitude [deg]

    real(kind=QP_K),parameter :: EPSLN=1.q-16

    if ( abs(x) < EPSLN .and. abs(y) < EPSLN ) then
      lon = 0.d0
      lat = sign( 90.0_QP_K, z )
      return
    end if

    lat = atan2(z, sqrt(x*x+y*y))
    lon = acos(x/sqrt(x*x+y*y))

    lat = lat*R2D
    lon = lon*R2D

    if (y<0) lon = 360.d0-lon

  end subroutine xyz2lonlat_quad_double


  !===================================================================================
  !> Convert (x, y, z) in the cartesian coordinate to (lat,lon).
  !!
  !> \note all arguments are quadruple (kind=QP_K)
  subroutine xyz2lonlat_quad(x, y, z, lon, lat)
    implicit none
    real(kind=QP_K), intent(IN) :: x   !< x in cartesian
    real(kind=QP_K), intent(IN) :: y   !< y in cartesian
    real(kind=QP_K), intent(IN) :: z   !< z in cartesian
    real(kind=QP_K), intent(OUT) :: lon !< longitude [deg]
    real(kind=QP_K), intent(OUT) :: lat !< latitude [deg]

    real(kind=QP_K),parameter :: EPSLN=1.q-16

    if ( abs(x) < EPSLN .and. abs(y) < EPSLN ) then
      lon = 0.d0
      lat = sign( 90.0_QP_K, z )
      return
    end if

    lat = atan2(z, sqrt(x*x+y*y))
    lon = acos(x/sqrt(x*x+y*y))

    lat = lat*R2DQ
    lon = lon*R2DQ

    if (y<0) lon = 360.0_QP_K-lon

  end subroutine xyz2lonlat_quad


  !===================================================================================
  !> Calc distance of two points on a sphere on a great circle.
  !!
  !! If optional argument `normalize` is present and is .true., result
  !! distance is normalized by the circumference of a circle,
  !! ie. res=[0,1], not res=[0,2*PI]
  !!
  !! \note renamed from get_length() and cal_great_circle_dist()
  function get_length_on_great_circle( lon1, lat1, lon2, lat2, normalize ) result(res)
    implicit none
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    logical,intent(in),optional :: normalize !< if .true. result is normalized.
    real(kind=DP_K) :: res

    !! \note
    !! Method below is not acceptable, caz calculation of
    !! inner_product is not precisely 1.0 even when both points are
    !! same, which cause non-zero distance between same points.
!!$  real(kind=DP_K) :: x1, y1, z1, x2, y2, z2
!!$  real(kind=DP_K) :: inner_product
!!$  call latlon2xyz(lat1, lon1, x1, y1, z1)
!!$  call latlon2xyz(lat2, lon2, x2, y2, z2)
!!$  inner_product = x1*x2+y1*y2+z1*z2
!!$  res = acos(inner_product)

    real(kind=DP_K) :: lo1, la1, lo2, la2 !! lon/lat in radian

    !!
    !! lon must be in 0...360.
    lo1 = mod(lon1,360.d0) ; if ( lo1 < 0.d0 ) lo1=lo1+360.d0
    lo2 = mod(lon2,360.d0) ; if ( lo2 < 0.d0 ) lo2=lo2+360.d0

    lo1 = lo1*D2R
    lo2 = lo2*D2R

    la1 = lat1*D2R
    la2 = lat2*D2R


    res = 2.d0*asin(sqrt(haversine(la2-la1)+cos(la1)*cos(la2)*haversine(lo2-lo1)))

    if ( present(normalize) ) then
      if ( normalize ) res = res*0.5d0/PI
    end if

  end function get_length_on_great_circle


  !===================================================================================
  !> Calc area of a spherical polygon formed by the arcs of great circles.
  !!
  !! Polygon is devided in several spherical triangles, each area of them is
  !! calc'ed by get_area_of_spherical_triangle().
  !!
  !! See Lauritzen 2008 appendix A.
  !!
  !! \note renamed from cal_great_circle_area().
  function get_area_of_spherical_polygon(np, lon, lat) result(res)
    implicit none
    integer, intent(IN) :: np !< number of point
    real(kind=DP_K), intent(IN) :: lat(:)  !< array of latitude
    real(kind=DP_K), intent(IN) :: lon(:)  !< array of longitude
    real(kind=DP_K) :: res
    real(kind=DP_K) :: area
    integer :: i

    area = 0.d0

    do i = 1, np-2
      area = area + get_area_of_spherical_triangle( &
           & lon(1), lat(1), &
           & lon(i+1), lat(i+1), &
           & lon(i+2), lat(i+2))
    end do

    res = area

  end function get_area_of_spherical_polygon


  !===================================================================================
  !> Calc area of a spherical triangle.
  !!
  !! See Lauritzen 2008 appendix A.
  !!
  !! \note renamed from cal_great_circle_triangle_area().
  function get_area_of_spherical_triangle( &
       & lon1, lat1,&
       & lon2, lat2,&
       & lon3, lat3 )&
       & result (res)
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]  
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg] 
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]  
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg] 
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]  
    real(kind=DP_K) :: res

    real(kind=DP_K) :: a, b, c, s
    real(kind=DP_K) :: temp1

!!$    a = cal_great_circle_dist(lat1*D2R, lon1*D2R, lat2*D2R, lon2*D2R)
!!$    b = cal_great_circle_dist(lat1*D2R, lon1*D2R, lat3*D2R, lon3*D2R)
!!$    c = cal_great_circle_dist(lat3*D2R, lon3*D2R, lat2*D2R, lon2*D2R)
    a = get_length_on_great_circle(lon1, lat1, lon2, lat2)
    b = get_length_on_great_circle(lon1, lat1, lon3, lat3)
    c = get_length_on_great_circle(lon3, lat3, lon2, lat2)

    s = 0.5d0*(a+b+c)
!!$#ifdef DEBUG
!!$    write(0,*)'dbg:get_area_of_spherical_triangle:a,b,c,s'
!!$    write(0,*)a,b,c,s
!!$#endif
    ! res = tan(0.5d0*s)*tan((s-a)*0.5d0)*tan((s-b)*0.5d0)*tan((s-c)*0.5d0)))
    temp1 =  max(tan(0.5d0*s)*tan((s-a)*0.5d0)*tan((s-b)*0.5d0)*tan((s-c)*0.5d0), 0.d0)
    res = 4.d0*atan(sqrt(temp1)) 

  end function get_area_of_spherical_triangle






  !!===================================================================================
  !> Calc (lon,lat) of an intersection point of 2 line segments 
  !! (lon1,lat1)--(lon2,lat2) and (lon3,lat3)--(lon4,lat4).
  !!
  !! If given two lines does NOT intersect, return (999,999).
  !!
  !! \note renamed from is_cross_line.
  !!
  !! In comparison of alpha,beta with 1 and 0, EPS1,EPS2 are now quadruple precision,
  !! So result may differ from former is_cross_line()
  !!
  subroutine get_intersection_point(&
       & lon, lat, do_intersect, &
       & lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) 
    implicit none
    real(kind=DP_K), intent(OUT) :: lat !< latitude of cross point [deg]
    real(kind=DP_K), intent(OUT) :: lon !< longitude of cross point [deg]
    logical,         intent(OUT) :: do_intersect !< .true. if intersecting.
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4 [deg]
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 [deg]

    real(kind=DP_K) :: x, y, z
    real(kind=QP_K) :: alpha, beta, t, length

    real(kind=QP_K) :: x1, y1, z1, x2, y2, z2

    real(kind=QP_K) :: EPS1=1.q-12
    real(kind=QP_K) :: EPS2=1.q-45


    call check_intersect( alpha, beta,&
         & lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4, &
         & x1, y1, z1, x2, y2, z2 )

    lon=999.d0
    lat=999.d0
    do_intersect = .false.
!!$    if ((1.0d0+EPS1 >= dble(alpha)).and.(dble(alpha) >= -EPS2)) then
!!$      if ((1.0d0+EP1 >= dble(beta)).and.(dble(beta) >= -EPS2)) then
    if ((EPS1 >= (alpha-1.0_QP_K)).and.(alpha >= -EPS2)) then
      if ((EPS1 >= (beta-1.0_QP_K)).and.(beta >= -EPS2)) then
          x = alpha*x1+(1.0_QP_K-alpha)*x2
          y = alpha*y1+(1.0_QP_K-alpha)*y2
          z = alpha*z1+(1.0_QP_K-alpha)*z2
          length = sqrt(x*x+y*y+z*z)
          x = x/length
          y = y/length
          z = z/length
          call xyz2lonlat(x,y,z,lon,lat)
          do_intersect=.true.
      end if
    end if

    return
  end subroutine get_intersection_point



  

  !!========================================================================
  !> Calculate center of spherical rectangle.
  !!
  !! \note renamed from cal_great_circle_center_rect()
  !!
  subroutine get_center_of_spherical_rectangle( &
       & lon, lat, &
       & lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4)

    real(kind=DP_K), intent(OUT) :: lon !< latitude of the center point [deg]
    real(kind=DP_K), intent(OUT) :: lat !< longitude of the center point [deg]
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4 [deg]
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 [deg]


    real(kind=DP_K) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    real(kind=DP_K) :: cx, cy, cz

!!$    if (.not.is_initialized) then
!!$      write(0,*) "init_sphere_lib is not called"
!!$      stop
!!$    end if

    call lonlat2xyz(lon1, lat1, x1, y1, z1)
    call lonlat2xyz(lon2, lat2, x2, y2, z2)
    call lonlat2xyz(lon3, lat3, x3, y3, z3)
    call lonlat2xyz(lon4, lat4, x4, y4, z4)

    cx = (x1+x2+x3+x4)*0.25d0
    cy = (y1+y2+y3+y4)*0.25d0
    cz = (z1+z2+z3+z4)*0.25d0

    call xyz2lonlat(cx, cy, cz, lon, lat)

    return
  end subroutine get_center_of_spherical_rectangle



  !===================================================================================
  !> Return .true. if the point is on the line defined by two points
  function is_on_the_line(lon1, lat1, lon2, lat2, lonc, latc) result(res)
    implicit none
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]      
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lonc !< longitude of target point [deg]
    real(kind=DP_K), intent(IN) :: latc !< latitude of target point [deg]
    logical :: res

    real(kind=QP_K) :: x1, y1, z1, x2, y2, z2, xc, yc, zc
    real(kind=QP_K) :: theta1, theta2, theta3
    real(kind=QP_K) :: EPSLN=1.e-12_QP_K

    call lonlat2xyz(lon1, lat1, x1, y1, z1)
    call lonlat2xyz(lon2, lat2, x2, y2, z2)
    call lonlat2xyz(lonc, latc, xc, yc, zc)

    theta1 = acos(x1*x2+y1*y2+z1*z2)
    theta2 = acos(x1*xc+y1*yc+z1*zc)
    theta3 = acos(x2*xc+y2*yc+z2*zc)
!!$#ifdef DEBUG
!!$    write(0,*) "dbg:is_on_the_line ", theta1, theta2, theta3, (theta2+theta3)-theta1
!!$#endif
!!$    is_on_line = (theta1 >= theta2+theta3-1.d-12)
    res = ( (theta2+theta3)-theta1 .le. EPSLN)
    
  end function is_on_the_line

  !===================================================================================
  !> Return .true. if two points are on the same spot.
  function is_the_same_point(lon1, lat1, lon2, lat2) result (res)
    implicit none
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    logical :: res 
    real(kind=DP_K) :: mlon1, mlon2
    real(kind=DP_K) :: distance

    real(kind=DP_K),parameter :: EPS1=1.0d-8
    real(kind=DP_K),parameter :: EPS2=1.0d-15

    !write(0,*) "is_same_point1 ", lat1, lon1, lat2, lon2

    mlon1 = mod(lon1, 360.d0); if ( mlon1 < 0.d0 ) mlon1 = mlon1 + 360.d0
    mlon2 = mod(lon2, 360.d0); if ( mlon2 < 0.d0 ) mlon2 = mlon2 + 360.d0

    if ((lat2==lat1).and.(mlon2==mlon1)) then
      res = .true.
      return
    end if

    !write(0,*) "is_same_point2 ", lat2-lat1, lon2-lon1

    if ((abs(lat2-lat1) <= EPS1 ).and.(abs(mlon2-mlon1) <= EPS1)) then
      res = .true.
      return
    end if

    distance = get_length_on_great_circle(mlon1, lat1, mlon2, lat2)

    !write(0,*) "is_same_point3 ", distance

    res = ( distance < EPS2 )

  end function is_the_same_point

  !===================================================================================
  !> Return .true. if two lines (lon1, lat1)-(lon2, lat2) and
  !! (lon3, lat3)-(lon4, lat4) are the same.
  function is_the_same_line(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4)&
       & result(res)
    implicit none
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4 [deg]
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 [deg]

    logical :: res
    !! below are quadraple precision.
    real(kind=QP_K) :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd
    real(kind=QP_K) :: matA(3,3), matI(3,3), datA
    real(kind=QP_K) :: EPSLN=1.0q-32

    call lonlat2xyz(lon1, lat1, xa, ya, za)
    call lonlat2xyz(lon2, lat2, xb, yb, zb)
    call lonlat2xyz(lon3, lat3, xc, yc, zc)
    call lonlat2xyz(lon4, lat4, xd, yd, zd)
    matA(1,1) = xb-xa ; matA(1,2) = xc-xd ; matA(1,3) = xd
    matA(2,1) = yb-ya ; matA(2,2) = yc-yd ; matA(2,3) = yd
    matA(3,1) = zb-za ; matA(3,2) = zc-zd ; matA(3,3) = zd

    call cal_inverse_matrix(matA, matI, datA)

    res=(abs(datA) < EPSLN )

  end function is_the_same_line

  !!===================================================================================
  !> Return .true. if two line segments (lon1, lat1)-(lon2, lat2) and
  !! (lon3,lat3)-(lon4, lat4) intersect.
  !!
  !! \note renamed from is_cross_line
  !! 
  function is_intersecting_lines(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4 ) result(res)
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4 [deg]
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 [deg]

    logical :: res

    real(kind=QP_K) :: alpha, beta
    real(kind=QP_K) :: EPS1=1.q-12
    real(kind=QP_K) :: EPS2=1.q-45

    res=.false.
    call check_intersect( alpha, beta,&
         & lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4)

    if ((EPS1 >= (alpha-1.0_QP_K)).and.(alpha >= -EPS2)) then
      if ((EPS1 >= (beta-1.0_QP_K)).and.(beta >= -EPS2)) then
        res=.true.
      end if
    end if
    
    return
  end function is_intersecting_lines



  !===================================================================================
  !> Return .true. if given point(in (lon,lat)) is in the rectangle
  !! given pair of lon and lat.
  !!
  !! Rewriten by T.Inoue 2014/08/20 from is_in_latlon()
  logical function is_in_lonlat(lon1, lat1, lon2, lat2, lon, lat)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon !< latitude of tested point [deg]
    real(kind=DP_K), intent(IN) :: lat !< longitude of tested point [deg]

    real(kind=DP_K) :: lonZ
    real(kind=DP_K) :: l1, l2, ll
    logical :: lon_inner, lat_inner

    lonZ=lon1

    l1=mod(lon1-lonZ,360.d0); if ( l1 < 0.d0 ) l1 = l1 + 360.d0
    l2=mod(lon2-lonZ,360.d0); if ( l2 < 0.d0 ) l2 = l2 + 360.d0
    ll=mod(lon -lonZ,360.d0); if ( ll < 0.d0 ) ll = ll + 360.d0

    lon_inner = ( l1 <= ll .and. ll < l2 )
    lat_inner = ( lat1 <= lat .and. lat < lat2 )

    is_in_lonlat = lon_inner .and. lat_inner

  end function is_in_lonlat




  !===================================================================================
  !> Return .true. if given point is in the spherical triangle formed on by three points.
  function is_inner(lon1, lat1, lon2, lat2, lon3, lat3, lon, lat) result(res)
    implicit none
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of triangle point 1 [deg]
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of triangle point 1 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of triangle point 2 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of triangle point 2 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of triangle point 3 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of triangle point 3 [deg]
    real(kind=DP_K), intent(IN) :: lon !< longitude of tested point [deg]
    real(kind=DP_K), intent(IN) :: lat !< latitude of tested point [deg]

    logical :: res

    real(kind=DP_K) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
    real(kind=DP_K) :: x, y, z
    real(kind=DP_K) :: detA
    real(kind=DP_K) :: alpha, beta, gamma

    res = .false.

    call lonlat2xyz(lat1, lon1, x1, y1, z1)
    call lonlat2xyz(lat2, lon2, x2, y2, z2)
    call lonlat2xyz(lat3, lon3, x3, y3, z3)
    call lonlat2xyz(lat,  lon,  x,  y,  z )

!!$    detA = x1*y2*z3+x3*y1*z2+x2*y3*z1-x3*y2*z1-x1*y3*z2-x2*y1*z3
    detA = sign(1.d0,&
         & (x1*y2*z3 +x3*y1*z2 +x2*y3*z1 -x3*y2*z1 -x1*y3*z2 -x2*y1*z3) )

!!$  alpha = ((y2*z3-y3*z2)*tx+(z2*x3-z3*x2)*ty+(x2*y3-x3*y2)*tz)/detA
!!$  beta = ((y3*z1-y1*z3)*tx+(z3*x1-z1*x3)*ty+(x3*y1-x1*y3)*tz)/detA
!!$  gamma = ((y1*z2-y2*z1)*tx+(z1*x2-z2*x1)*ty+(x1*y2-x2*y1)*tz)/detA

    alpha = sign(1.d0,&
         & ((y2*z3-y3*z2)*x +(z2*x3-z3*x2)*y +(x2*y3-x3*y2)*z) )
    if (alpha*detA<0.d0) return

    beta  = sign(1.d0,&
         & ((y3*z1-y1*z3)*x +(z3*x1-z1*x3)*y +(x3*y1-x1*y3)*z) )
    if (beta*detA<0.d0) return

    gamma = sign(1.d0,&
         & ((y1*z2-y2*z1)*x +(z1*x2-z2*x1)*y +(x1*y2-x2*y1)*z) )
    if (gamma*detA<0.d0) return

    res = .true.

  end function is_inner

!!$========================================================================
!!$ Private routines
!!$========================================================================

  !!========================================================================
  !> [Private] calc `alpha` and `beta` used in line intersection checking.
  !!
  subroutine check_intersect( alpha, beta, &
       & lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4,&
       & x1, y1, z1, x2, y2, z2 )
    real(kind=QP_K),intent(out) :: alpha, beta
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 
    real(kind=QP_K),intent(out),optional :: x1, y1, z1, x2, y2, z2

    real(kind=QP_K) :: xa, ya, za, xb, yb, zb
    real(kind=QP_K) :: xc, yc, zc, xd, yd, zd
    real(kind=QP_K) :: matA(3,3), matI(3,3), datA
    real(kind=QP_K) :: t

    call lonlat2xyz(lon1, lat1, xa, ya, za)
    call lonlat2xyz(lon2, lat2, xb, yb, zb)
    call lonlat2xyz(lon3, lat3, xc, yc, zc)
    call lonlat2xyz(lon4, lat4, xd, yd, zd)
    matA(1,1) = xb-xa ; matA(1,2) = xc-xd ; matA(1,3) = xd
    matA(2,1) = yb-ya ; matA(2,2) = yc-yd ; matA(2,3) = yd
    matA(3,1) = zb-za ; matA(3,2) = zc-zd ; matA(3,3) = zd

    call cal_inverse_matrix(matA, matI, datA)

    alpha = (xb*matI(1,1)+yb*matI(1,2)+zb*matI(1,3))/datA
    t = xb*matI(3,1)+yb*matI(3,2)+zb*matI(3,3)
    beta = (xb*matI(2,1)+yb*matI(2,2)+zb*matI(2,3))/t

    if ( present( x1 ) ) x1 = xa
    if ( present( y1 ) ) y1 = ya
    if ( present( z1 ) ) z1 = za
    if ( present( x2 ) ) x2 = xb
    if ( present( y2 ) ) y2 = yb
    if ( present( z2 ) ) z2 = zb

    return
  end subroutine check_intersect


  !===================================================================================
  !> Calc inverse matrix (for double precision)
  !!
  !! Matrix A,I are defined as follows;
  !!
  !!      a11 a12 a13
  !!      a21 a22 a23
  !!      a31 a32 a33
  subroutine cal_inverse_matrix_double(A, I, datA)
    implicit none
    real(kind=DP_K), intent(IN) :: A(3, 3)  !< 3x3 matrix
    real(kind=DP_K), intent(OUT) :: I(3, 3) !< inverse of matrix A
    real(kind=DP_K), intent(OUT) :: datA    !< a determinant of A

    datA =  A(1,1)*A(2,2)*A(3,3) &
         & +A(2,1)*A(3,2)*A(1,3) &
         & +A(3,1)*A(1,2)*A(2,3) &
         & -A(1,1)*A(3,2)*A(2,3) &
         & -A(3,1)*A(2,2)*A(1,3) &
         & -A(2,1)*A(1,2)*A(3,3)

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
  !> Calc inverse matrix (for quadruple precision)
  !!
  !! Matrix A,I are defined as follows;
  !!
  !!      a11 a12 a13
  !!      a21 a22 a23
  !!      a31 a32 a33
  subroutine cal_inverse_matrix_quad(A, I, datA)
    implicit none
    real(kind=QP_K), intent(IN) :: A(3, 3)  !< 3x3 matrix
    real(kind=QP_K), intent(OUT) :: I(3, 3) !< inverse of matrix A
    real(kind=QP_K), intent(OUT) :: datA    !< a determinant of A

    datA =  A(1,1)*A(2,2)*A(3,3) &
         & +A(2,1)*A(3,2)*A(1,3) &
         & +A(3,1)*A(1,2)*A(2,3) &
         & -A(1,1)*A(3,2)*A(2,3) &
         & -A(3,1)*A(2,2)*A(1,3) &
         & -A(2,1)*A(1,2)*A(3,3)

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



  !===================================================================================
  !> Calc haversine;
  !!
  !! \f[\text{hav} x \equiv \frac{1-\cos x}{2}\f]
  !!
  !! \attention `theta` here is in __RADIAN__
  function haversine(theta) result(res)
    implicit none
    real(kind=DP_K), intent(IN) :: theta !< degree
    real(kind=DP_K) :: res

    res = (1.d0-cos(theta))*0.5d0

  end function haversine




end module jcf_spherical_lib
