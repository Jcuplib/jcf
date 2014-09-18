!> Wrapper module for backward compatibility.
!!
!! This module is now obsolete and use jcf_spherical_lib instead.
!!
!! \copyright
!!   The MIT License (MIT),  (c) 2013-2014 RIST.
!!   See the [LICENSE](@ref license) file.
!!   
!! 
!! \author Takashi ARAKAWA <arakawa@rist.jp>
!!
!! \version
!!   - 2013/XX/XX T.Arakawa  Initial public version
!!   - 2014/08/29 T.Inoue    Refactored.
!!   .
!!
module jcf_sphere_lib
  use jcf_spherical_lib
  implicit none

  private !! is default
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


!!$  interface cal_inverse_matrix
!!$    module procedure cal_inverse_matrix_double, cal_inverse_matrix_quad
!!$  end interface

  !> genelic subroutine to convert (lat,lon) to cartesian coordinate (X,Y,Z)
  interface latlon2xyz
    module procedure latlon2xyz_double, latlon2xyz_double_quad, latlon2xyz_quad
  end interface

  !> genelic subroutine to convert cartesian coordinate (X,Y,Z) to (lat,lon)
  interface xyz2latlon
    module procedure xyz2latlon_double, xyz2latlon_quad_double, xyz2latlon_quad
  end interface

contains

!===================================================================================
!> initialize this library
  subroutine init_sphere_lib()
    implicit none

    write(0,*)'calling init_sphere_lib, this is OBSOLETE!'
    call init_spherical_lib()
  end subroutine init_sphere_lib


  !===================================================================================
  !> \deprecated Now call get_area_of_spherical_polygon().
  function cal_great_circle_area(num_of_point, lat, lon) result(res)
    real(kind=DP_K) :: res
    integer, intent(IN) :: num_of_point !< number of point
    real(kind=DP_K), intent(IN) :: lat(:)  !< array of latitude
    real(kind=DP_K), intent(IN) :: lon(:)  !< array of longitude

    write(0,*)'calling cal_great_circle_area, this is OBSOLETE!'
    res=get_area_of_spherical_polygon(num_of_point, lon, lat )
  end function cal_great_circle_area



  !===================================================================================
  !> \deprecated Now call get_area_of_spherical_triangle()
  real(kind=DP_K) function cal_great_circle_triangle_area(lat1, lon1, lat2, lon2, lat3, lon3)
    real(kind=DP_K), intent(IN) :: lat1
    real(kind=DP_K), intent(IN) :: lon1
    real(kind=DP_K), intent(IN) :: lat2
    real(kind=DP_K), intent(IN) :: lon2
    real(kind=DP_K), intent(IN) :: lat3
    real(kind=DP_K), intent(IN) :: lon3

    write(0,*)'calling cal_great_circle_triangle_area, this is OBSOLETE!'

    cal_great_circle_triangle_area=get_area_of_spherical_triangle(&
         & lon1, lat1, lon2, lat2, lon3, lat3 )

  end function cal_great_circle_triangle_area

  !===================================================================================
  !> \deprecated Now call get_length_on_great_circle()
  !!
  !!
  real(kind=DP_K) function cal_great_circle_dist(lat1, lon1, lat2, lon2)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1, lon1, lat2, lon2

    write(0,*)'calling cal_great_circle_dist, this is OBSOLETE!'
    cal_great_circle_dist=get_length_on_great_circle(&
         & lon1/D2R, lat1/D2R, lon2/D2R, lat2/D2R)

  end function cal_great_circle_dist


  !===================================================================================
  !> \deprecated Now call get_length_on_great_circle().
  real(kind=DP_K) function get_length(lat1, lon1, lat2, lon2)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2
    write(0,*)'calling get_length, this is OBSOLETE!'
    get_length=get_length_on_great_circle( lon1, lat1, lon2, lat2, normalize=.true.)
  end function get_length

  !===================================================================================
  !> \deprecated Now call is_on_the_line().
  function is_on_line(lat1, lon1, lat2, lon2, latc, lonc) result(res)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]      
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: latc !< latitude of target point [deg]
    real(kind=DP_K), intent(IN) :: lonc !< longitude of target point [deg]
    logical :: res

    write(0,*)'calling is_on_line, this is OBSOLETE!'
    res=is_on_the_line(lon1, lat1, lon2, lat2, lonc, latc)
  end function is_on_line


  !===================================================================================
  !> \deprecated Now call is_the_same_point().
  function is_same_point(lat1, lon1, lat2, lon2) result (res)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    logical :: res 

    write(0,*)'calling is_same_point, this is OBSOLETE!'
    res=is_the_same_point(lon1, lat1, lon2, lat2)
  end function is_same_point
    


  !===================================================================================
  !> \deprecated Now call is_the_same_line().
  function is_same_line(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4) result (res)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 [deg]
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4 [deg]
    logical :: res

    write(0,*)'calling is_same_point, this is OBSOLETE!'
    res=is_the_same_line(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4)

  end function is_same_line


  !===================================================================================
  !> \deprecated Now call is_intersecting_lines() or get_intersection_point().
  logical function is_cross_line(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, cross_lat, cross_lon)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4
    real(kind=DP_K), intent(OUT),optional :: cross_lat !< latitude of cross point
    real(kind=DP_K), intent(OUT),optional :: cross_lon !< longitude of cross point 

    write(0,*)'calling is_cross_line, this is OBSOLETE!'
    call get_intersection_point(cross_lon, cross_lat,&
       & lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) 
  end function is_cross_line



  !===================================================================================
  !> \deprecated Now call get_center_of_spherical_rectangle()
  subroutine cal_great_circle_center_rect(&
       & lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4, clat, clon)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1 [deg]
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2 [deg]
    real(kind=DP_K), intent(IN) :: lat3 !< latitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lon3 !< longitude of point3 [deg]
    real(kind=DP_K), intent(IN) :: lat4 !< latitude of point4 [deg]
    real(kind=DP_K), intent(IN) :: lon4 !< longitude of point4 [deg]
    real(kind=DP_K), intent(OUT) :: clon !< latitude of the center point [deg]
    real(kind=DP_K), intent(OUT) :: clat !< longitude of the center point [deg]


    write(0,*)'calling cal_great_circle_center_rect, this is OBSOLETE!'
    call get_center_of_spherical_rectangle(&
         & clon, clat, &
         & lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4)
  end subroutine cal_great_circle_center_rect


  !!========================================================================
  !> \deprecated Now call is_in_lonlat().
  logical function is_in_latlon(lon1, lat1, lon2, lat2, tlon, tlat)
    implicit none
    real(kind=DP_K), intent(IN) :: lat1 !< latitude of point1 
    real(kind=DP_K), intent(IN) :: lon1 !< longitude of point1
    real(kind=DP_K), intent(IN) :: lat2 !< latitude of point2 
    real(kind=DP_K), intent(IN) :: lon2 !< longitude of point2
    real(kind=DP_K), intent(IN) :: tlon !< latitude of tested point
    real(kind=DP_K), intent(IN) :: tlat !< longitude of tested point

    write(0,*)'calling is_in_latlon, this is OBSOLETE!'
    is_in_latlon=is_in_lonlat(lon1,lat1,lon2,lat2,tlon,tlat)
  end function is_in_latlon




!===================================================================================
!> \depricated converte (lat,lon) to (x, y, z) in the cartesian coordinate.
subroutine latlon2xyz_double(lat, lon, x, y, z)
  implicit none
  real(kind=8), intent(IN) :: lat !< latitude [deg]
  real(kind=8), intent(IN) :: lon !< longitude [deg]
  real(kind=8), intent(OUT) :: x !< x in cartesian
  real(kind=8), intent(OUT) :: y !< y in cartesian
  real(kind=8), intent(OUT) :: z !< z in cartesian

    write(0,*)'calling latlon2xyz_double, this is OBSOLETE!'
    call lonlat2xyz(lon,lat,x,y,z)
end subroutine latlon2xyz_double



!===================================================================================
!> \depricated converte (lat,lon) to (x, y, z) in the cartesian coordinate.
subroutine latlon2xyz_double_quad(lat, lon, x, y, z)
  implicit none
  real(kind=8), intent(IN) :: lat !< latitude [deg]
  real(kind=8), intent(IN) :: lon !< longitude [deg]
  real(kind=16), intent(OUT) :: x !< x in cartesian
  real(kind=16), intent(OUT) :: y !< y in cartesian
  real(kind=16), intent(OUT) :: z !< z in cartesian

    write(0,*)'calling latlon2xyz_double_quad, this is OBSOLETE!'
    call lonlat2xyz(lon,lat,x,y,z)
end subroutine latlon2xyz_double_quad


!===================================================================================
!> \deprecated converte (lat,lon) to (x, y, z) in the cartesian coordinate.
subroutine latlon2xyz_quad(lat, lon, x, y, z)
  implicit none
  real(kind=16), intent(IN) :: lat !< latitude [deg]
  real(kind=16), intent(IN) :: lon !< longitude [deg]
  real(kind=16), intent(OUT) :: x  !< x in cartesian
  real(kind=16), intent(OUT) :: y  !< y in cartesian
  real(kind=16), intent(OUT) :: z  !< z in cartesian

    write(0,*)'calling latlon2xyz_quad, this is OBSOLETE!'
    call lonlat2xyz(lon,lat,x,y,z)
end subroutine latlon2xyz_quad




!===================================================================================
!> \deprecated converte (x, y, z) in the cartesian coordinate to (lat,lon).
subroutine xyz2latlon_double(x, y, z, lat, lon)
  implicit none
  real(kind=8), intent(IN) :: x    !< x in cartesian
  real(kind=8), intent(IN) :: y    !< y in cartesian
  real(kind=8), intent(IN) :: z    !< z in cartesian
  real(kind=8), intent(OUT) :: lat !< latitude [deg]
  real(kind=8), intent(OUT) :: lon !< longitude [deg]


  write(0,*)'calling xyz2latlon_double, this is OBSOLETE!'
  call xyz2lonlat(x,y,z,lon,lat)
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

  write(0,*)'calling xyz2latlon_quad_double, this is OBSOLETE!'
  call xyz2lonlat(x,y,z,lon,lat)
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

  write(0,*)'calling xyz2latlon_quad, this is OBSOLETE!'
  call xyz2lonlat(x,y,z,lon,lat)
end subroutine xyz2latlon_quad



!!$========================================================================
!!$========================================================================
!!$========================================================================
!!$========================================================================
































end module jcf_sphere_lib
