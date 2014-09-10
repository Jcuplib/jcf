module mod_coef
  use jcf_mesh_base, only : polygon_type
  implicit none
  private

  !--------------------------------   public  ----------------------------------!

  public :: cal_monte_carlo_points ! subroutine (polygon, depth)
  public :: cal_monte_carlo_coef ! subroutine (polygon)

  !--------------------------------   private  ---------------------------------!

  real(kind=8) :: ER = 6370000.d0 ! earth radius
  !real(kind=8) :: ER = 1.d0 ! earth radius
  real(kind=8) :: PI, D2R
  integer, parameter :: MAX_POINT = 32

contains

  !=======+=========+=========+=========+=========+=========+=========+=========+
  !> \todo
  !! Intel compiler on ha8000 complains that:
  !!
  !!      fortcom: Error: mod_coef.f90, line 54: A pointer dummy argument
  !!      with the INTENT(IN) attribute shall not appear as an actual argument
  !!      if the associated dummy argument has the INTENT(OUT) or INTENT(INOUT)
  !!      attribute.  [POLYGON]
  !!          call cal_monte_carlo_points_north_pole(polygon, depth)
  !!
  !! But, there's no pointer with INTENT(IN) in this module.
  !! This seems to be a bug of compiler.
  !!
  !! To avoid this, remove INTENT(IN) tentatively.
  !!
  subroutine cal_monte_carlo_points(polygon, depth)
    implicit none
!!$    type(polygon_type), pointer, intent(INOUT) :: polygon ! my polygon
    type(polygon_type), pointer :: polygon ! my polygon
    integer, intent(IN) :: depth
    integer :: i, j, p
    real(kind=8) :: clon, clat
    real(kind=8), allocatable :: lon(:), lat(:)
    integer :: num_of_point
    real(kind=8), allocatable :: xi(:), yi(:)

    num_of_point = polygon%num_of_point

    !write(0,*) "cal_coef 1 ", num_of_point

    allocate(lon(num_of_point))
    allocate(lat(num_of_point))

    ! set point position
    do i = 1, num_of_point
      lon(i) = polygon%point(i)%ptr%x
      lat(i) = polygon%point(i)%ptr%y
      !write(0,*) i, lon(i), lat(i)
    end do

    clon = polygon%data_point%x
    clat = polygon%data_point%y

    if (clat > maxval(lat)) then ! north pole
      write(0,*) "north pole polygon"
      call cal_monte_carlo_points_north_pole(polygon, depth)
      return
    end if

    if (clat < minval(lat)) then ! south pole
      write(0,*) "south pole polygon"
      call cal_monte_carlo_points_south_pole(polygon, depth)
      return
    end if

    call cal_monte_carlo_points_normal(polygon, depth)

    deallocate(lon, lat)

  end subroutine cal_monte_carlo_points

  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine cal_monte_carlo_points_north_pole(polygon, depth)
    implicit none
    type(polygon_type), pointer, intent(INOUT) :: polygon ! my polygon
    integer, intent(IN) :: depth
    integer :: i, j, p
    type(polygon_type), pointer :: hex_polygon, target_polygon
    real(kind=8) :: clon, clat
    real(kind=8), allocatable, dimension(:) :: lon, lat, y_pos, x, y
    integer :: num_of_point
    real(kind=8), allocatable :: xi(:), yi(:)
    real(kind=8), allocatable :: loni(:), lati(:)
    integer :: num_of_inner_point
    integer :: num_of_total_point
    real(kind=8) :: PI, D2R, len

    PI = atan(1.d0)*4.d0
    D2R =PI/180.d0

    num_of_point = polygon%num_of_point

    !write(0,*) "cal_coef north pole 1 ", num_of_point

    allocate(lon(num_of_point))
    allocate(lat(num_of_point))
    allocate(x(num_of_point))
    allocate(y(num_of_point))

    ! set point position
    do i = 1, num_of_point
      lon(i) = polygon%point(i)%ptr%x
      lat(i) = polygon%point(i)%ptr%y
      !write(0,*) i, lon(i), lat(i)
    end do

    clon = polygon%data_point%x
    clat = polygon%data_point%y

    do i = 1, num_of_point
      len = (90.d0-lat(i))*2*PI*ER/360.d0
      x(i) = len*cos(lon(i)*D2R)
      y(i) = len*sin(lon(i)*D2R)

      !write(0,*) i, x(i), y(i)
    end do


    num_of_inner_point = (2**depth)*((2**depth)-1)/2 ! number of inner points in a triangle area

    num_of_total_point = num_of_point*num_of_inner_point ! number of inner points in the polygon

    !write(0,*) "cal_coef 3 ", num_of_inner_point, num_of_total_point

    allocate(xi(num_of_total_point))
    allocate(yi(num_of_total_point))
    allocate(loni(num_of_total_point))
    allocate(lati(num_of_total_point))

    do i = 1, num_of_point-1
      call cal_inner_points(0.d0, 0.d0, x(i), y(i), x(i+1), y(i+1), &
           xi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), &
           yi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), depth) 
    end do

    i = num_of_point
    call cal_inner_points(0.d0, 0.d0, x(i), y(i), x(1), y(1), &
         xi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), &
         yi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), depth) 

    !write(0,*) "cal_coef 4 ", num_of_inner_point, num_of_total_point

    do i = 1, num_of_total_point

      len = sqrt(xi(i)*xi(i)+yi(i)*yi(i))
      lati(i) = 90.d0-len*180.d0/(PI*ER)
      loni(i) = atan2(yi(i), xi(i))*180/PI
      if (yi(i)<0) then
        loni(i) = 360.d0+loni(i)
      end if

      !write(0,*) "x   y   ", xi(i), yi(i)
      !write(0,*) "lon lat ", polygon%monte_carlo_point(i)%point%x, polygon%monte_carlo_point(i)%point%y
    end do

    call cal_num_of_inner_point(polygon, loni, lati, num_of_total_point, 1.d0)

    deallocate(lon, lat)
    deallocate(x, y)

    deallocate(xi, yi)

    deallocate(loni, lati)


  end subroutine cal_monte_carlo_points_north_pole

  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine cal_monte_carlo_points_south_pole(polygon, depth)
    implicit none
    type(polygon_type), pointer, intent(INOUT) :: polygon ! my polygon
    integer, intent(IN) :: depth
    integer :: i, j, p
    type(polygon_type), pointer :: hex_polygon, target_polygon
    real(kind=8) :: clon, clat
    real(kind=8), allocatable, dimension(:) :: lon, lat, y_pos, x, y
    integer :: num_of_point
    real(kind=8), allocatable :: xi(:), yi(:), loni(:), lati(:)
    integer :: num_of_inner_point
    integer :: num_of_total_point
    real(kind=8) :: PI, D2R, len

    PI = atan(1.d0)*4.d0
    D2R =PI/180.d0

    num_of_point = polygon%num_of_point

    !write(0,*) "cal_coef north pole 1 ", num_of_point

    allocate(lon(num_of_point))
    allocate(lat(num_of_point))
    allocate(x(num_of_point))
    allocate(y(num_of_point))

    ! set point position
    do i = 1, num_of_point
      lon(i) = polygon%point(i)%ptr%x
      lat(i) = polygon%point(i)%ptr%y
      !write(0,*) i, lon(i), lat(i)
    end do

    clon = polygon%data_point%x
    clat = polygon%data_point%y

    do i = 1, num_of_point
      len = (lat(i)+90.d0)*2*PI*ER/360.d0
      x(i) = len*cos(lon(i)*D2R)
      y(i) = len*sin(lon(i)*D2R)

      !write(0,*) i, x(i), y(i)
    end do

    num_of_inner_point = (2**depth)*((2**depth)-1)/2 ! number of inner points in a triangle area

    num_of_total_point = num_of_point*num_of_inner_point ! number of inner points in the polygon

    !write(0,*) "cal_coef 3 ", num_of_inner_point, num_of_total_point

    allocate(xi(num_of_total_point))
    allocate(yi(num_of_total_point))
    allocate(loni(num_of_total_point))
    allocate(lati(num_of_total_point))

    do i = 1, num_of_point-1
      call cal_inner_points(0.d0, 0.d0, x(i), y(i), x(i+1), y(i+1), &
           xi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), &
           yi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), depth) 
    end do

    i = num_of_point
    call cal_inner_points(0.d0, 0.d0, x(i), y(i), x(1), y(1), &
         xi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), &
         yi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), depth) 

    !write(0,*) "cal_coef 4 ", num_of_inner_point, num_of_total_point

    do i = 1, num_of_total_point
      len = sqrt(xi(i)*xi(i)+yi(i)*yi(i))
      lati(i) = -90.d0+len*180.d0/(PI*ER)
      loni(i) = atan2(yi(i), xi(i))*180/PI
      if (yi(i)<0) then
        loni(i) = 360.d0+loni(i)
      end if

      !write(0,*) "x   y   ", i, xi(i), yi(i)
      !write(0,*) "lon lat ", polygon%monte_carlo_point(i)%point%x, polygon%monte_carlo_point(i)%point%y
    end do

    call cal_num_of_inner_point(polygon,loni, lati, num_of_total_point, 1.d0)

    deallocate(xi, yi)
    deallocate(loni, lati)
    deallocate(lon, lat)
    deallocate(x, y)

  end subroutine cal_monte_carlo_points_south_pole

  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine cal_monte_carlo_points_normal(polygon, depth)
    implicit none
    type(polygon_type), pointer, intent(INOUT) :: polygon ! my polygon
    integer, intent(IN) :: depth
    integer :: i, j, p
    type(polygon_type), pointer :: hex_polygon, target_polygon
    real(kind=8) :: clon, clat
    real(kind=8), allocatable, dimension(:) :: lon, lat, y_pos, x, y
    integer :: num_of_point
    real(kind=8), allocatable :: xi(:), yi(:)
    real(kind=8), allocatable :: loni(:), lati(:)
    integer :: num_of_inner_point
    integer :: num_of_total_point

    num_of_point = polygon%num_of_point

    !write(0,*) "cal_coef 1 ", num_of_point

    allocate(lon(num_of_point))
    allocate(lat(num_of_point))
    allocate(x(num_of_point))
    allocate(y(num_of_point))

    ! set point position
    do i = 1, num_of_point
      lon(i) = polygon%point(i)%ptr%x
      lat(i) = polygon%point(i)%ptr%y
      !write(0,*) i, lon(i), lat(i)
    end do

    clon = polygon%data_point%x
    clat = polygon%data_point%y

    if (minval(lon) == 0.d0) then
      if (maxval(lon) >= 180.d0) then
        do i = 1, num_of_point
          if (lon(i) == 0.d0) lon(i) = 360.d0
        end do
      end if
    end if

    if (maxval(lon) == 360.d0) then
      if (minval(lon) <= 180.d0) then
        do i = 1, num_of_point
          if (lon(i) == 360.d0) lon(i) = 0.d0
        end do
      end if
    end if

    if (maxval(lon)-minval(lon) >= 180.d0) then
      do i = 1, num_of_point
        if (lon(i) > 180.d0) lon(i) = lon(i) - 360.d0
      end do
    end if


    !do i = 1, num_of_point
    !  write(0,*) i, lon(i), lat(i)
    !end do

    !write(0,*) "cal_coef 2 ", clon, clat

    call cal_local_xy(clon, clat, lon, lat, x, y)

    num_of_inner_point = (2**depth)*((2**depth)-1)/2 ! number of inner points in a triangle area

    num_of_total_point = num_of_point*num_of_inner_point ! number of inner points in the polygon

    !write(0,*) "cal_coef 3 ", num_of_inner_point, num_of_total_point


    allocate(xi(num_of_total_point))
    allocate(yi(num_of_total_point))
    allocate(loni(num_of_total_point))
    allocate(lati(num_of_total_point))

    do i = 1, num_of_point-1
      call cal_inner_points(0.d0, 0.d0, x(i), y(i), x(i+1), y(i+1), &
           xi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), &
           yi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), depth) 
    end do

    i = num_of_point
    call cal_inner_points(0.d0, 0.d0, x(i), y(i), x(1), y(1), &
         xi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), &
         yi(num_of_inner_point*(i-1)+1:num_of_inner_point*i), depth) 

    !write(0,*) "cal_coef 4 ", num_of_inner_point, num_of_total_point

    do i = 1, num_of_total_point

      call cal_global_lat_lon(clon, clat, xi(i), yi(i), loni(i), lati(i))

      !write(0,*) i, xi(i), yi(i)
      !write(0,*)  polygon%monte_carlo_point(i)%point%x, polygon%monte_carlo_point(i)%point%y
    end do

    call cal_num_of_inner_point(polygon, loni, lati, num_of_total_point, 1.d0)

    deallocate(xi, yi)
    deallocate(loni, lati)
    deallocate(x, y)
    deallocate(lon, lat)

  end subroutine cal_monte_carlo_points_normal

  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine cal_num_of_inner_point(polygon, loni, lati, num_of_monte_carlo_point, weight)
    use jcf_mesh_base, only : get_target_polygon_by_num, &
         set_target_polygon_coef1_by_num, set_target_polygon_coef1_by_index
    use jcf_sphere_lib, only : is_inner, is_in_latlon
    implicit none
    type(polygon_type), pointer :: polygon
    real(kind=8), intent(IN) :: loni(:), lati(:)
    integer, intent(IN) :: num_of_monte_carlo_point
    real(kind=8), intent(IN) :: weight
    type(polygon_type), pointer :: target_polygon
    real(kind=8) :: x1, y1, x2, y2, x3, y3
    integer, allocatable :: area(:)
    integer :: i, j, p
    logical :: inner_flag 

    allocate(area(polygon%num_of_target))

    area(:) = 0

    !do j = 1, polygon%num_of_target
    !  target_polygon => get_target_polygon_by_num(polygon, j)
    !  write(0,*) target_polygon%index
    !  write(0,*) target_polygon%point(1)%ptr%x, target_polygon%point(1)%ptr%y
    !  write(0,*) target_polygon%point(3)%ptr%x, target_polygon%point(3)%ptr%y
    !
    !end do

    do i = 1, num_of_monte_carlo_point

      inner_flag = .false.

      target_loop: do j = 1, polygon%num_of_target

        target_polygon => get_target_polygon_by_num(polygon, j)

        if (target_polygon%is_latlon) then
          x1 = target_polygon%point(1)%ptr%x
          y1 = target_polygon%point(1)%ptr%y
          x2 = target_polygon%point(3)%ptr%x
          y2 = target_polygon%point(3)%ptr%y
          !write(0,*) x1, y1, x2, y2

          if (x1 == 360.d0) x1 = 0.d0
          if (x2 == 0.d0) x2 = 360.d0
          if (is_in_latlon(x1, y1, x2, y2, loni(i), lati(i))) then
            area(j) = area(j) + 1
            !write(0,*) j, area(j)
            inner_flag = .true.
            exit target_loop
          end if
        else
          x1 = target_polygon%point(1)%ptr%x
          y1 = target_polygon%point(1)%ptr%y
          do p = 2, target_polygon%num_of_point-1
            x2 = target_polygon%point(p)%ptr%x
            y2 = target_polygon%point(p)%ptr%y
            x3 = target_polygon%point(p+1)%ptr%x
            y3 = target_polygon%point(p+1)%ptr%y

            !write(432, *) x1, y1
            !write(432, *) x2, y2
            !write(432, *) x3, y3
            !write(432, *) is_inner(x1, y1, x2, y2, x3, y3, point_pos_x, point_pos_y)
            !write(432, *)

            if (is_inner(x1, y1, x2, y2, x3, y3, loni(i), lati(i))) then
              area(j) = area(j) + 1
              inner_flag = .true.
              !write(0,*) j, area(j)
              exit target_loop
            end if

          end do
        end if

      end do target_loop

      if ( .not. inner_flag ) then
        write(0, *) "no target polygon ::", polygon%index, i
        !do j = 1, polygon%num_of_point
        !   write(0,*) polygon%point(j)%ptr%x, polygon%point(j)%ptr%y
        !end do
        !write(0, *) loni, lati
        !write(0, *) "   num_of_target : ", polygon%num_of_target
        !do j = 1, polygon%num_of_target
        !  target_polygon => get_target_polygon_by_num(polygon, j)
        !  write(0,*) "target :: ", j, target_polygon%index
        !  write(0,*) target_polygon%point(1)%ptr%x, target_polygon%point(1)%ptr%y
        !  write(0,*) target_polygon%point(2)%ptr%x, target_polygon%point(2)%ptr%y
        !  write(0,*) target_polygon%point(3)%ptr%x, target_polygon%point(3)%ptr%y
        !  write(0,*) target_polygon%point(4)%ptr%x, target_polygon%point(4)%ptr%y
        !  write(0,*)
        !end do
      end if

    end do

    do i = 1, polygon%num_of_target
      call set_target_polygon_coef1_by_num(polygon, i, area(i)*weight)
      target_polygon => get_target_polygon_by_num(polygon, i)
      call set_target_polygon_coef1_by_index(target_polygon, polygon%index, area(i)*weight)
      !write(0,*) "cal_num_of_inner_point ", i, target_polygon%index, polygon%index, area(i)
    end do

    deallocate(area)

  end subroutine cal_num_of_inner_point

  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine cal_monte_carlo_coef(polygon)
    use jcf_mesh_base, only : get_target_polygon_coef1_by_num, set_target_polygon_coef1_by_num, &
         get_target_polygon_by_num
    implicit none
    type(polygon_type), pointer :: polygon
    integer :: i
    real(kind=8) :: total_num, coef
    type(polygon_type), pointer :: target_polygon

    total_num = 0

    do i = 1, polygon%num_of_target
      target_polygon => get_target_polygon_by_num(polygon, i)
      if (target_polygon%mask) then
        total_num = total_num + get_target_polygon_coef1_by_num(polygon, i)
      end if
    end do

    if (total_num == 0) then
      polygon%mask = .false.
      return
    end if

    do i = 1, polygon%num_of_target
      coef = get_target_polygon_coef1_by_num(polygon, i)/total_num
      call set_target_polygon_coef1_by_num(polygon, i, coef)
    end do

  end subroutine cal_monte_carlo_coef


  !===================================================================================

  subroutine cal_locaL_xy(clon, clat, lon, lat, x, y)
    implicit none
    real(kind=8), intent(IN) :: clon, clat
    real(kind=8), intent(IN) :: lon(:), lat(:)
    real(kind=8), intent(OUT) :: x(:), y(:)
    integer :: i

    PI = atan(1.d0)*4.d0
    D2R =PI/180.d0

    do i = 1, size(lat)
      x(i) = 2*PI*ER*(lon(i)-clon)/360.d0*cos(lat(i)*D2R)
      y(i) = 2*PI*ER*(lat(i)-clat)/360.d0
      !write(0,*) "cal_local_xy ", lon(i), lat(i), x(i), y(i)
    end do

  end subroutine cal_local_xy


  !===================================================================================

  subroutine cal_global_lat_lon(clon, clat, x, y, lon, lat)
    implicit none
    real(kind=8), intent(IN) :: clon, clat
    real(kind=8), intent(IN) :: x, y
    real(kind=8), intent(OUT) :: lon, lat

    PI = atan(1.d0)*4.d0
    D2R =PI/180.d0

    lat = clat + y*180.d0/(PI*ER)
    lon = clon + x*180.d0/(PI*ER)*1.d0/cos(lat*D2R)

  end subroutine cal_global_lat_lon


  !===================================================================================

  subroutine cal_inner_points(xo, yo, xa, ya, xb, yb, xi, yi, depth)
    implicit none
    real(kind=8), intent(IN) :: xo, yo, xa, ya, xb, yb
    real(kind=8), intent(INOUT) :: xi(:), yi(:)
    integer, intent(IN) :: depth
    real(kind=8) :: x1, y1, x2, y2
    real(kind=8) :: dx1, dy1, dx2, dy2
    integer :: num_of_line, num_of_points, total_points
    integer :: i

    num_of_line = 2**depth

    total_points = num_of_line*(num_of_line-1)/2

    dx1 =2*(xa-xo)/(2*num_of_line-1)
    dx2 =2*(xb-xo)/(2*num_of_line-1)
    dy1 =2*(ya-yo)/(2*num_of_line-1)
    dy2 =2*(yb-yo)/(2*num_of_line-1)

    do i = 1, num_of_line - 1
      num_of_points = i
      x1 = xo+dx1*i
      y1 = xo+dy1*i
      x2 = xo+dx2*i
      y2 = yo+dy2*i
      !write(0,*) "cal_inner_points ", x1, y1, x2, y2
      call cal_points(x1, y1, x2, y2, xi(i*(i-1)/2+1:i*(i+1)/2), yi(i*(i-1)/2+1:i*(i+1)/2), num_of_points)    
    end do

  end subroutine cal_inner_points


  !===================================================================================

  subroutine cal_points(x1, y1, x2, y2, x, y, num_of_points)
    implicit none
    real(kind=8), intent(IN) :: x1, y1, x2, y2
    real(kind=8), intent(INOUT) :: x(:), y(:)
    integer, intent(IN) :: num_of_points
    real(kind=8) :: dx, dy
    integer :: i


    dx = (x2-x1)/num_of_points
    dy = (y2-y1)/num_of_points

    do i = 1, num_of_points
      x(i) = x1 + dx*(i-1)
      y(i) = y1 + dy*(i-1)
      !write(0,*) "cal_points ", x(i), y(i)
    end do

  end subroutine cal_points

  !===================================================================================

end module mod_coef
