!=======+=========+=========+=========+=========+=========+=========+=========+
!> This module difines API for constructing mesh structure, calculating mapping table
!! and result output.
!! 
!! Making mesh structure is the first stage of calculating mapping table and 
!! interpolation coefficient. 
!! An example code for constructiong mesh structure is listed below.
!! In this example, two mesh, mesh1 and mesh2 is defined. It is assumed that both mesh have
!! global latitude-longitude grid. The number of grid point (data point) of mesh1 and mesh2
!! are (NX1, NY1) and (NX2, NY2) respectively. Total number of vertexes on global lat-lon grid
!! is (NX+1)*(NY+1).
!!
!! \attention
!! Only the case of mesh1 is listed in the example code below.
!! 
!! \attention
!! The example of next polygon setting listed below neglects boundary condision. See user's guide 
!! for more detail.
!! 
!! \code
!!  type(mesh_type) :: mesh1, mesh2
!!   
!!  ! initialize jcf_sphere_lib
!!  call init_shere_lib() 
!!    
!!  ! set number of polygon, number of center point, number of vertexes
!!   call init_mesh(mesh1, NX1*NY1, NX1*NY1, (NX1+1)*(NY1+1)) 
!!    
!!  ! set shape of every polygon and define that the sides of polygon are along lat-lon or not
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      index = i + NX1*(j-1)
!!      call init_polygon(mesh1, index, 4) ! initialize a polygon of mesh1 as rectangle
!!      call set_latlon(mesh1, index, .true.)
!!    end do
!!  end do
!!    
!!  ! set center location of a polygon indexed by i + NX1*(j-1)
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      call set_data_point_location(mesh1, i+NX1*(j-1), x(i,j), y(i,j)) ! set center location of polygon
!!    end do
!!  end do
!!    
!!  ! set location of vertex indexed by i+(NX1+1)*(j-1) 
!!  do j = 1, NY1+1
!!    do i = 1, NX1+1
!!      call set_volume_point_location(mesh1, i+(NX1+1)*(j-1), vx(i,j), vy(i,j)) ! set vertex location of mesh1
!!    end do
!!  end do
!!   
!!  ! set center point and vertexes to the polygon indexed by i + NX1*(j-1)
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      index = i+NX1*(j-1)
!!      call set_data_point(mesh1, index, i+NX1*(j-1), .true.)
!!      call set_volume_point(mesh1, index, i  +(NX1+1)*(j  ))
!!      call set_volume_point(mesh1, index, i  +(NX1+1)*(j-1))
!!      call set_volume_point(mesh1, index, i+1+(NX1+1)*(j-1))
!!      call set_volume_point(mesh1, index, i+1+(NX1+1)*(j  ))
!!    end do
!!  end do
!!   
!!  ! set next polygon
!!  !             (i,j+1)
!!  !(i-1, j)     (i, j )    (i+1, j)
!!  !             (i,j-1)
!!  ! this example neglected boundary setting, see user's guide for more detail
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      index = i+NX1*(j-1)
!!      call set_next_polygon(mesh1, index, i-1+NX1*(j-1))
!!      call set_next_polygon(mesh1, index, i  +NX1*(j-2))
!!      call set_next_polygon(mesh1, index, i+1+NX1*(j-1))
!!      call set_next_polygon(mesh1, index, i  +NX1*(j  ))
!!    end do
!!  end do
!! \endcode
!!   
!!
!! Second stage is overlapped polygon search. This stage is dividet into two phase,
!! The first phase is to search a polygon of which vertex is included in the target polygon. 
!! The second phase is to search a polygon which is overlapped but no vetex is included in the target polygon.
!!
!! \note
!! Both direction mesh1 -> mesh2 and mesh2 -> mesh1 search is necessary for mapping table calculation.
!!
!! \code  
!!  type(polygon_ptr_type), pointer :: start_polygon(:)
!!  type(polygon_type), pointer :: polygon_ptr
!!  integer :: I, j
!!
!!  allocate(start_polygon(1))
!! 
!!  ! The first phase is to search a polygon of whic vertex is include in the target polygon.
!!  !     _____
!!  !  __|_____|___
!!  ! |  |     |   |
!!  ! |  |_____|   |
!!  ! |____________|   
!!  !    
!!  !    
!!  start_polygon(1)%ptr => get_polygon_ptr(mesh2, 1)
!! 
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      polygon_ptr => get_polygon_ptr(mesh1, i + NX1*(j-1))
!!      call search_polygon(polygon_ptr, start_polygon)
!!    end do
!!  end do
!! 
!!  start_polygon(1)%ptr => get_polygon_ptr(mesh1, 1)
!! 
!!  do j = 1, NY2
!!    do i = 1, NX2
!!      polygon_ptr => get_polygon_ptr(mesh2, i + NX2*(j-1))
!!      call search_polygon(polygon_ptr, start_polygon)
!!    end do
!!  end do
!!
!!  ! The second phase is to search a polygon which is overlapped but no vetex is included in the target polygon.
!!  !     _____
!!  !  __|_____|___
!!  ! |  |     |   |
!!  ! |  |     |   |
!!  ! |__|_____|___|  
!!  !    |_____|
!!  !    
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      polygon_ptr => get_polygon_ptr(mesh1, i + NX1*(j-1))
!!      call search_polygon_by_side(polygon_ptr)
!!    end do
!!  end do
!! 
!!  do j = 1, NY2
!!    do i = 1, NX2
!!      polygon_ptr => get_polygon_ptr(mesh2, i + NX2*(j-1))
!!      call search_polygon_by_side(polygon_ptr)
!!    end do
!!  end do
!! \endcode
!!
!! The third stage, the result is output. 
!!
!! \code
!!  type(polygon_type), pointer :: mesh1_polygon, mesh2_polygon
!!  real(kind=8) :: coef
!!  integer :: counter
!!  integer :: i, j, t
!! 
!!  counter = 0
!!  do j = 1, NY1
!!    do i = 1, NX1
!!      mesh1_polygon => get_polygon_ptr(mesh1, i + NX1*(j-1))
!!      do t = 1, mesh1_polygon%num_of_target
!!        mesh2_polygon => get_target_polygon_by_num(mesh1_polygon, t)
!!        counter = counter + 1
!!        coef = get_target_polygon_coef_by_num(mesh1_polygon, t)
!!        write(FILE_UNIT, *) counter, mesh1_polygon%index, mesh2_polygon%index, coef
!!      end do
!!    end do
!!  end do
!! \endcode
!!
!!
!!
!> @author Takashi ARAKAWA <arakawa@rist.jp>
!
!------------------------------------------------------------------------------------
module jcf_mesh_base
  implicit none

!--------------------------------   public  ----------------------------------!

  public :: mesh_type
  public :: polygon_type
  public :: polygon_ptr_type
  public :: target_polygon_type

  public :: init_mesh
  public :: set_data_point_location
  public :: set_volume_point_location
  
  public :: init_polygon
  public :: set_data ! subroutine (my_mesh, polygon_index, data)
  public :: get_data ! real(kind=8) function (my_mesh, polygon_index)
  public :: set_mask ! subroutine (my_mesh, polygon_index, mask)
  public :: get_mask ! logical function (my_mesh, polygon_index)
  public :: set_latlon ! subroutine (my_mesh, polygon_index, is_latlon)
  public :: is_latlon  ! logical function (my_mesh, polygon_index)
  public :: set_data_point
  public :: set_volume_point
  public :: set_next_polygon
  public :: get_data_point_x ! real(kind=8) function (my_mesh, polygon_index)
  public :: get_data_point_y ! real(kind=8) function (my_mesh, polygon_index)
  public :: get_num_of_point ! integer function (my_mesh, polygon_index)
  public :: get_point_x ! real(kind=8) function (my_mesh, polygon_index, point_num)
  public :: get_point_y ! real(kind=8) function (my_mesh, polygon_index, point_num)
  public :: get_polygon_ptr
  public :: search_polygon         ! subroutine (this_polygon, start_polygon, lat_lon_latitude)
  public :: search_polygon_by_side ! subroutine (this_polygon)
  public :: get_target_polygon_by_num ! polygon_type, pointer function (my_polygon, target_polygon_num)
  public :: set_target_polygon_coef1_by_num ! subroutine (my_polygon, target_polygon_num, coef)
  public :: get_target_polygon_coef1_by_num ! real(kind=8) function (my_polygon, target_polygon_num)
  public :: set_target_polygon_coef1_by_index ! subroutine (my_polygon, target_polygon_index, coef)
  public :: get_target_polygon_coef1_by_index ! real(kind=8) function (my_polygon, target_polygon_index)
  public :: set_target_polygon_coef2_by_num ! subroutine (my_polygon, target_polygon_num, coef)
  public :: get_target_polygon_coef2_by_num ! real(kind=8) function (my_polygon, target_polygon_num)
  public :: set_target_polygon_coef2_by_index ! subroutine (my_polygon, target_polygon_index, coef)
  public :: get_target_polygon_coef2_by_index ! real(kind=8) function (my_polygon, target_polygon_index)
  public :: get_target_polygon_by_index ! polygon_type, pointer function (my_polygon, target_polygon_index)
  public :: reset_monitor
  public :: write_monitor_info
  public :: write_polygon_info


  private :: point_type
  private :: point_ptr_type
  private :: side_type
  private :: monte_carlo_point_type
  private :: monitor_type

  private :: search_next_polygon
  private :: search_polygon_from_point
  private :: is_in_this_polygon
  private :: add_target_polygon
  private :: check_mask_correspondence
  private :: check_polygon_mask

!--------------------------------   private  ---------------------------------!

  !> point_type defines a point.
  type point_type
    integer :: index            !<  ID number
    real(kind=8), pointer :: x  !< location in x direction
    real(kind=8), pointer :: y  !< location in y direction
    type(polygon_type), pointer :: target_polygon !< pointer to the polygon in which the point is included
  end type point_type

  type point_ptr_type
    type(point_type), pointer :: ptr
  end type point_ptr_type

  type side_type 
    type(point_ptr_type) :: point(2)
  end type side_type

  !> array of polygon_type
  type polygon_ptr_type
    type(polygon_type), pointer :: ptr !< pointer to polygon_type.
  end type polygon_ptr_type

  !> target_polygon_type defines list of overlapped polygon. 
  type target_polygon_type
    real(kind=8) :: coef1 !< interpolation coefficient:  `r = r + coef*target\%data`
    real(kind=8) :: coef2 !< interpolation coefficient:  `r = r + coef*target\%data`
    integer :: area !< overlapped area
    type(polygon_type), pointer :: target !< pointer to my polygon.
    type(target_polygon_type), pointer :: before !< pointer to the polygon before me.
    type(target_polygon_type), pointer :: next   !< pointer to the polygon next to me.
  end type target_polygon_type

  type monte_carlo_point_type
    type(point_type) :: point
    integer :: polygon_index
  end type monte_carlo_point_type

  !> polygon_type defines a polygon.  
  type polygon_type
    logical :: search_flag      !< flag whether search is completed or not.
    logical :: mask             !< flag whether this polygon is valid or not.
    logical :: is_latlon        !< flag whether the sides of this polygon are along latitude-longitude.
    integer :: index            !< ID number of this polygon.
    real(kind=8) :: data        !< physical value of this polygon (not used).
    type(point_type), pointer :: data_point !< point where the physical value is defined in this polygon.
    integer :: num_of_point                 !< number of vertex.
    integer :: point_set_counter            !<
    type(point_ptr_type), pointer :: point(:) !< array of vertex
    real(kind=8) :: min_x !< min of X of this area (not used)
    real(kind=8) :: max_x !< max of X of this area (not used)
    real(kind=8) :: min_y !< min of Y of this area (not used)
    real(kind=8) :: max_y !< max of Y of this area (not used)
    real(kind=8) :: min_xp !< min of X of vertexes 
    real(kind=8) :: max_xp !< max of X of vertexes 
    real(kind=8) :: min_yp !< min of Y of vertexes 
    real(kind=8) :: max_yp !< max of Y of vertexes 
    integer :: num_of_next !< number of neighbor polygon
    integer :: next_set_counter
    type(polygon_ptr_type), pointer :: next_polygon(:) !< pointer array of neighbor polygon
    integer :: num_of_target                           !< number of overlapped polygon
    type(target_polygon_type), pointer :: target_polygon !< pointer to the list of overlapped polygon
    !integer :: num_of_monte_carlo_point
    !type(monte_carlo_point_type), pointer :: monte_carlo_point(:)
    logical :: is_included !< flag whether this polygon is included in the target polygon
  end type polygon_type

  type(polygon_type), private, pointer :: current_polygon

  !> mesh_type defines a mesh.
  type mesh_type
    integer :: num_of_data_point !< number of data point ( = num_of_polygon )
    type(point_type), pointer :: data_mesh(:) !< array of grid point
    integer :: num_of_volume_point            !< number of vertexes
    type(point_type), pointer :: volume_mesh(:) !< array of vertexes
    integer :: num_of_polygon                   !< number of polygon ( = num_of_data_point )
    type(polygon_type), pointer :: polygon(:)   !< array of polygon
  end type mesh_type

  type monitor_type
    integer :: search_counter
    integer :: level_max
    integer :: count_max
    integer :: max_count_index
    integer :: count_first
    real(kind=8) :: count_mean
  end type monitor_type

  type(monitor_type),private :: monitor

  integer, parameter :: MAX_POINT = 32 !

contains

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> initialize mesh
subroutine init_mesh(this, num_of_polygon, num_of_data_point, num_of_volume_point) 
  implicit none
  type(mesh_type), intent(INOUT) :: this !< this mesh.
  integer, intent(IN) :: num_of_polygon  !< number of polygon.
  integer, intent(IN) :: num_of_data_point !< number of data point.
  integer, intent(IN) :: num_of_volume_point !< number of volume point.

  this%num_of_data_point = num_of_data_point
  allocate(this%data_mesh(num_of_data_point))
  this%num_of_volume_point = num_of_volume_point
  allocate(this%volume_mesh(num_of_volume_point))
  this%num_of_polygon = num_of_polygon
  allocate(this%polygon(num_of_polygon))

end subroutine init_mesh

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set data point location
subroutine set_data_point_location(this, point_index, x, y)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< this mesh
  integer, intent(IN) :: point_index     !< index of data point
  real(kind=8), intent(IN), target :: x  !< data location
  real(kind=8), intent(IN), target :: y  !< data location

  this%data_mesh(point_index)%x => x
  this%data_mesh(point_index)%y => y
  this%data_mesh(point_index)%index = point_index
  Nullify(this%data_mesh(point_index)%target_polygon)

end subroutine set_data_point_location

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set volume point location
subroutine set_volume_point_location(this, point_index, x, y)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< this mesh
  integer, intent(IN) :: point_index     !< indexof data point
  real(kind=8), intent(IN), target :: x !< data location
  real(kind=8), intent(IN), target :: y !< data location

  this%volume_mesh(point_index)%x => x
  this%volume_mesh(point_index)%y => y
  this%volume_mesh(point_index)%index = point_index
  Nullify(this%volume_mesh(point_index)%target_polygon)

end subroutine set_volume_point_location


!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> initialize polygon
subroutine init_polygon(this, polygon_index, num_of_point)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< polygon_type variable
  integer, intent(IN) :: polygon_index   !< index of polygon
  integer, intent(IN) :: num_of_point    !< number of points of the polygon
  integer :: i

  current_polygon => this%polygon(polygon_index)

  current_polygon%index = polygon_index
  current_polygon%mask = .false.
  current_polygon%is_latlon = .false.
  current_polygon%num_of_point = num_of_point
  current_polygon%num_of_next = num_of_point
  current_polygon%search_flag = .false. !!d
  allocate(current_polygon%point(num_of_point))
  allocate(current_polygon%next_polygon(num_of_point))

  do i = 1, num_of_point
    nullify(current_polygon%point(i)%ptr)
    nullify(current_polygon%next_polygon(i)%ptr)
  end do

  current_polygon%point_set_counter = 0
  current_polygon%next_set_counter = 0

  current_polygon%num_of_target = 0 
  nullify(current_polygon%target_polygon)

  current_polygon%min_x = 99999999.d0
  current_polygon%max_x = -99999999.d0
  current_polygon%min_y = 99999999.d0
  current_polygon%max_y = -99999999.d0
  current_polygon%min_xp = 99999999.d0
  current_polygon%max_xp = -99999999.d0
  current_polygon%min_yp = 99999999.d0
  current_polygon%max_yp = -99999999.d0

end subroutine init_polygon

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set value on the data point 
subroutine set_data(this, polygon_index, dt)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< this mesh
  integer, intent(IN) :: polygon_index   !< index of polygon
  real(kind=8), intent(IN) :: dt         !< value of data point

  this%polygon(polygon_index)%data = dt
  
end subroutine set_data

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> get value on the data point 
real(kind=8) function get_data(this, polygon_index)
  implicit none
  type(mesh_type), intent(IN) :: this !< this mesh
  integer, intent(IN) :: polygon_index !<  index of polygon

  get_data = this%polygon(polygon_index)%data

end function get_data

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set mask
subroutine set_mask(this, polygon_index, mask)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< mesh_type variable
  integer, intent(IN) :: polygon_index   !< index of polygon
  logical, intent(IN) :: mask            !< mask flag

  this%polygon(polygon_index)%mask = mask
  
end subroutine set_mask

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> get mask
logical  function get_mask(this, polygon_index)
  implicit none
  type(mesh_type), intent(IN) :: this !< mesh_type variable
  integer, intent(IN) :: polygon_index !< polygon index

  get_mask = this%polygon(polygon_index)%mask

end function get_mask

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set latlon
subroutine set_latlon(this, polygon_index, is_latlon)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< mesh_type variable
  integer, intent(IN) :: polygon_index   !< index of polygon
  logical, intent(IN) :: is_latlon       !< flag whether th sides are along latitude-longitude

  this%polygon(polygon_index)%is_latlon = is_latlon
  
end subroutine set_latlon

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set latlon
logical  function is_latlon(this, polygon_index)
  implicit none
  type(mesh_type), intent(IN) :: this  !< mesh_type variable
  integer, intent(IN) :: polygon_index !< index of polygon

  is_latlon = this%polygon(polygon_index)%is_latlon

end function is_latlon

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set data point
subroutine set_data_point(this, polygon_index, point_index, mask)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< mesh_type variable
  integer, intent(IN) :: polygon_index   !< index of polygon
  integer, intent(IN) :: point_index     !< index of data point
  logical, intent(IN) :: mask            !< flag whether this polygon is valid or not
  
  this%polygon(polygon_index)%data_point => this%data_mesh(point_index)
  this%polygon(polygon_index)%mask = mask

end subroutine set_data_point

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set volume point
subroutine set_volume_point(this, polygon_index, point_index)
  implicit none
  type(mesh_type), intent(INOUT) :: this !< mesh_type variable
  integer, intent(IN) :: polygon_index   !< index of polygon
  integer, intent(IN) :: point_index     !< index of data point
  real(kind=8) :: x, y
  real(kind=8) :: temp_x

  current_polygon => this%polygon(polygon_index)

  current_polygon%point_set_counter = current_polygon%point_set_counter+1
  if (current_polygon%point_set_counter > current_polygon%num_of_point) then
    write(0,*) "point set counter over !!!"
    stop
  end if

  current_polygon%point(current_polygon%point_set_counter)%ptr => this%volume_mesh(point_index)

  x = this%volume_mesh(point_index)%x
  y = this%volume_mesh(point_index)%y

  if (current_polygon%min_x > x) current_polygon%min_x = x
  if (current_polygon%max_x < x) current_polygon%max_x = x
  if (current_polygon%min_y > y) current_polygon%min_y = y
  if (current_polygon%max_y < y) current_polygon%max_y = y

  if (current_polygon%min_xp > x) current_polygon%min_xp = x
  if (current_polygon%max_xp < x) current_polygon%max_xp = x
  if (current_polygon%min_yp > y) current_polygon%min_yp = y
  if (current_polygon%max_yp < y) current_polygon%max_yp = y

  if (current_polygon%point_set_counter == current_polygon%num_of_point) then ! last volume point
    if (abs(current_polygon%max_x-current_polygon%min_x) >= 180.d0) then
      temp_x = current_polygon%min_x
      current_polygon%min_x = current_polygon%max_x - 360.d0
      current_polygon%max_x = temp_x
      temp_x = current_polygon%min_xp
      current_polygon%min_xp = current_polygon%max_xp - 360.d0
      current_polygon%max_xp = temp_x
    end if
  end if

end subroutine set_volume_point

!=======+=========+=========+=========+=========+=========+=========+=========+
!===================================================================================
!> set next polygon
subroutine set_next_polygon(this, polygon_index, next_polygon_index)
  implicit none
  type(mesh_type), intent(INOUT) :: this    !< mesh_type variable
  integer, intent(IN) :: polygon_index      !< index of target polygon
  integer, intent(IN) :: next_polygon_index !< index of next polygon
  
  if (next_polygon_index <= 0) return

  current_polygon => this%polygon(polygon_index)

  current_polygon%next_set_counter = current_polygon%next_set_counter+1
  if (current_polygon%point_set_counter > current_polygon%num_of_point) then
    write(0,*) "next polygon set counter over !!!"
    stop
  end if

  current_polygon%next_polygon(current_polygon%next_set_counter)%ptr => this%polygon(next_polygon_index)

end subroutine set_next_polygon

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get data point x
real(kind=8) function get_data_point_x(this, polygon_index)
  implicit none
  type(mesh_type), intent(IN) :: this  !< mesh_type variable
  integer, intent(IN) :: polygon_index !< index of polygon  
  get_data_point_x = this%polygon(polygon_index)%data_point%x

end function get_data_point_x

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get data point y
real(kind=8) function get_data_point_y(this, polygon_index)
  implicit none
  type(mesh_type), intent(IN) :: this  !< mesh_type variable
  integer, intent(IN) :: polygon_index !< index of polygon  

  get_data_point_y = this%polygon(polygon_index)%data_point%y

end function get_data_point_y

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get the number of vertex
integer function get_num_of_point(this, polygon_index)
  implicit none
  type(mesh_type), intent(IN) :: this  !< mesh_type variable
  integer, intent(IN) :: polygon_index !< index of polygon

  get_num_of_point = this%polygon(polygon_index)%num_of_point

end function get_num_of_point

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get volume point x
real(kind=8) function get_point_x(this, polygon_index, point_num)
  implicit none
  type(mesh_type), intent(IN) :: this  !< mesh_type variable
  integer, intent(IN) :: polygon_index !< index of polygon
  integer, intent(IN) :: point_num     !< index of vertex   

  get_point_x = this%polygon(polygon_index)%point(point_num)%ptr%x

end function get_point_x

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get volume point y
real(kind=8) function get_point_y(this, polygon_index, point_num)
  implicit none
  type(mesh_type), intent(IN) :: this  !< mesh_type variable
  integer, intent(IN) :: polygon_index !< index of polygon  
  integer, intent(IN) :: point_num     !< index of vertex   

  get_point_y = this%polygon(polygon_index)%point(point_num)%ptr%y

end function get_point_y

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get pointer of the polygon
function get_polygon_ptr(this, polygon_index) result(polygon)
  implicit none
  type(mesh_type), intent(IN) :: this    !< mesh_type variable
  integer, intent(IN) :: polygon_index   !< index of polygon

  type(polygon_type), pointer :: polygon 

  polygon => this%polygon(polygon_index)

end function get_polygon_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+
!> reset performance monitor
subroutine reset_monitor()
  implicit none

  monitor%search_counter = 0
  monitor%level_max = 0
  monitor%count_max = 0
  monitor%count_first = 0
  monitor%count_mean = 0.d0

end subroutine reset_monitor

!=======+=========+=========+=========+=========+=========+=========+=========+
!> increment performance monitor
subroutine inc_search_counter(max_count, max_level, polygon_index)
  implicit none
  integer, intent(IN) :: max_count      !< max of search count
  integer, intent(IN) :: max_level      !< max of search level
  integer, intent(IN) :: polygon_index  !< indexof the polygon.

  if (monitor%count_first <= 0) monitor%count_first = max_count
  if (monitor%level_max < max_level) monitor%level_max = max_level
  if (monitor%count_max < max_count) then
    monitor%count_max = max_count
    monitor%max_count_index = polygon_index
  end if
  monitor%count_mean = monitor%count_mean + max_count
  monitor%search_counter = monitor%search_counter + 1

end subroutine inc_search_counter

!=======+=========+=========+=========+=========+=========+=========+=========+
!> search polygons overlapped on this polygon
subroutine search_polygon(this_polygon, start_polygon, lat_lon_latitude)
  implicit none
  type(polygon_type), pointer :: this_polygon            !< target polygon
  type(polygon_ptr_type), pointer :: start_polygon(:)    !< start polygon
  real(kind=8), optional, intent(IN) :: lat_lon_latitude !< latitude, from here the sides of polygon are along latitude-longitude

  type(polygon_type), pointer :: target_polygon
  logical :: is_lat_lon
  real(kind=8) :: x, y
  logical :: is_included
  integer :: search_count, max_level
  integer :: target_length
  integer :: before_index
  integer :: p, pp
  integer :: ii, jj, ll

      do p = 1, this_polygon%num_of_point   
 
        if (associated(this_polygon%point(p)%ptr%target_polygon)) then
          call add_target_polygon(this_polygon, this_polygon%point(p)%ptr%target_polygon)
          call add_target_polygon(this_polygon%point(p)%ptr%target_polygon, this_polygon)
          start_polygon(1)%ptr => this_polygon%point(p)%ptr%target_polygon
          cycle
        end if

        target_length = 99999
 
        do pp = p-1, 1, -1
          if ((associated(this_polygon%point(pp)%ptr%target_polygon)).and.((p-pp) < target_length)) then
            target_length = p-pp
            start_polygon(1)%ptr => this_polygon%point(pp)%ptr%target_polygon
            exit
          end if
        end do

        do pp = p+1, this_polygon%num_of_point
          if ((associated(this_polygon%point(pp)%ptr%target_polygon)).and.((pp-p) < target_length)) then
            start_polygon(1)%ptr => this_polygon%point(pp)%ptr%target_polygon
            exit
          end if
        end do

        x = this_polygon%point(p)%ptr%x
        y = this_polygon%point(p)%ptr%y
       
       is_lat_lon = .false.
       if (present(lat_lon_latitude)) then
         if (y <= lat_lon_latitude) then
           is_lat_lon = .true.
         end if
       end if

        !write(0,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        !write(0,*) "search polygon start "
        !write(0,*) " search point : ", x, y
        !write(0,*) " start target index : ", start_polygon(1)%ptr%index
        !write(0,*) " is_lat_lon         : ", is_lat_lon 

        search_count = 0
        call search_polygon_from_point(x, y, 1, start_polygon, target_polygon, 1, search_count, max_level, is_lat_lon)

        call inc_search_counter(search_count, max_level, this_polygon%index)

        this_polygon%point(p)%ptr%target_polygon => target_polygon

        !write(0,*) 
        !write(0,*) "target polygon index ", target_polygon%index
        !write(0,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "

        call add_target_polygon(this_polygon, target_polygon)
        call add_target_polygon(target_polygon, this_polygon)
        start_polygon(1)%ptr => target_polygon

      end do

      is_included = .true.
      before_index = 0

      do p = 1, this_polygon%num_of_point
        if ((p >= 2).and.(is_included)) then
          is_included = (before_index == this_polygon%point(p)%ptr%target_polygon%index)
        end if
        before_index = this_polygon%point(p)%ptr%target_polygon%index
      end do

      this_polygon%is_included = is_included

end subroutine search_polygon

!=======+=========+=========+=========+=========+=========+=========+=========+
!> search polygons overlapped on this polygon at the side
subroutine search_polygon_by_side(this_polygon)
  implicit none
  type(polygon_type), pointer :: this_polygon !< target polygon

  type(polygon_type), pointer :: current_polygon, next_polygon
  real(kind=8) :: xs, ys, xe, ye
  real(kind=8) :: x, y
  real(kind=8) :: delta_x, delta_y 
  integer :: search_count
  integer :: p, ps, pe
  integer :: ii, jj, ll


      ploop: do p = 1, this_polygon%num_of_point   

        ps = p
        pe = mod(p, this_polygon%num_of_point) + 1

        xs = this_polygon%point(ps)%ptr%x
        ys = this_polygon%point(ps)%ptr%y
        xe = this_polygon%point(pe)%ptr%x
        ye = this_polygon%point(pe)%ptr%y

        if (abs(xe-xs) < 180.d0) then
          delta_x = (xe-xs)/100.d0
        else
          if (xs>=xe) then
            delta_x = (xe-xs+360.d0)/100.d0
          else
            delta_x = (xe-xs-360.d0)/100.d0
          end if
        end if

        delta_y = (ye-ys)/100.d0

       current_polygon => this_polygon%point(p)%ptr%target_polygon
       if (.not.associated(current_polygon)) cycle ! next ploop

       if (.not.associated(this_polygon%point(ps)%ptr%target_polygon)) cycle
       if (.not.associated(this_polygon%point(pe)%ptr%target_polygon)) cycle

       ! if both ps and pe are included in the same polygon, skip search 
       if (this_polygon%point(ps)%ptr%target_polygon%index == this_polygon%point(pe)%ptr%target_polygon%index) cycle

       search_count = 0

       searchloop: do 

         search_count = search_count + 1
         if (search_count > 1280) then
           write(0,*) "seach_polygon_by_side, search count overflow ", this_polygon%index, p
           stop
         end if

         call search_next_polygon(xs, ys, xe, ye, x, y, current_polygon, next_polygon)
         if (.not.associated(next_polygon)) exit searchloop

         if (this_polygon%point(pe)%ptr%target_polygon%index == next_polygon%index) exit searchloop

         if (.not.next_polygon%is_included) then
           call add_target_polygon(this_polygon, next_polygon)
           call add_target_polygon(next_polygon, this_polygon)
         end if
     
         current_polygon => next_polygon
         xs = x + delta_x
         ys = y + delta_y

       end do searchloop

      end do ploop

end subroutine search_polygon_by_side

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine search_next_polygon(xs, ys, xe, ye, x, y, current_polygon, next_polygon)
  use jcf_sphere_lib, only : is_cross_line
  implicit none
  real(kind=8), intent(IN) :: xs, ys, xe, ye
  real(kind=8), intent(OUT) :: x, y
  type(polygon_type), pointer :: current_polygon
  type(polygon_type), pointer :: next_polygon
  real(kind=8) :: x1, y1, x2, y2
  integer :: i1, i2

  do i1 = 1, current_polygon%num_of_point
    i2 = mod(i1, current_polygon%num_of_point) + 1
    x1 = current_polygon%point(i1)%ptr%x
    y1 = current_polygon%point(i1)%ptr%y
    x2 = current_polygon%point(i2)%ptr%x
    y2 = current_polygon%point(i2)%ptr%y
    if (is_cross_line(ys, xs, ye, xe, y1, x1, y2, x2, y, x)) then
      next_polygon => current_polygon%next_polygon(i1)%ptr
      return
    end if
  end do

  nullify(next_polygon)
  x = 9999.d0
  y = 9999.d0

end subroutine search_next_polygon

!=======+=========+=========+=========+=========+=========+=========+=========+

recursive subroutine search_polygon_from_point(mx, my, p_num, polygons, target_polygon, level, &
                                               search_count, max_level, is_lat_lon)
  implicit none
  real(kind=8), intent(IN) :: mx, my
  integer,      intent(IN) :: p_num
  type(polygon_ptr_type), pointer :: polygons(:)
  type(polygon_type), pointer :: target_polygon
  integer, intent(IN) :: level
  integer, intent(INOUT) :: search_count
  integer, intent(INOUT) :: max_level
  logical, intent(IN) :: is_lat_lon

  integer :: next_p_num
  type(polygon_ptr_type), pointer :: next_polygons(:)
  integer :: i, j
  integer :: idx
  integer :: ii,jj,ll

  if (level>1550) then
  !if (level>5) then
    write(0,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(0,*) "                    search_polygon: level over",level  
    write(0,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    stop
  end if

  !write(0,*) "search_polygon start:",level,p_num

  next_p_num = 0
  do i = 1, p_num

    !iidx = polygons(i)%ptr%index
    !ii = mod(idx-1, 32)+1
    !ll = int((idx-1)/(32*32))+1
    !jj = int((idx-32*32*(ll-1)-1)/32)+1

    !write(0,*) "  search target :",i, polygons(i)%ptr%index, ii,jj,ll
    search_count = search_count+1

    !write(0,*)'dbg:search_polygon_from_point:'
    !write(0,*)polygons(i)%ptr%point(1)%ptr%x,polygons(i)%ptr%point(1)%ptr%y
    !write(0,*)polygons(i)%ptr%point(2)%ptr%x,polygons(i)%ptr%point(2)%ptr%y
    !write(0,*)polygons(i)%ptr%point(3)%ptr%x,polygons(i)%ptr%point(3)%ptr%y
    !write(0,*)polygons(i)%ptr%point(4)%ptr%x,polygons(i)%ptr%point(4)%ptr%y
    !write(0,*)polygons(i)%ptr%point(5)%ptr%x,polygons(i)%ptr%point(5)%ptr%y


    if (is_in_this_polygon(mx, my, polygons(i)%ptr, is_lat_lon)) then
      max_level = level
      target_polygon => polygons(i)%ptr
      return
    else
      polygons(i)%ptr%search_flag = .true.
      do j = 1, polygons(i)%ptr%num_of_next
        if (associated(polygons(i)%ptr%next_polygon(j)%ptr)) then
          !!if (polygons(i)%ptr%next_polygon(j)%ptr%search_flag==.false.) then
          if (.not. polygons(i)%ptr%next_polygon(j)%ptr%search_flag) then
            next_p_num=next_p_num+1
          end if
        end if
      end do
    end if
  end do

  !write(0,*) "  search_polygon next_p_num:",next_p_num

  allocate(next_polygons(next_p_num))

  next_p_num = 0
  do i = 1, p_num
      do j = 1, polygons(i)%ptr%num_of_next
        if (associated(polygons(i)%ptr%next_polygon(j)%ptr)) then
           !write(0,*) "   search_polygon next_polygon info ", &
           !               polygons(i)%ptr%index, &
           !               polygons(i)%ptr%next_polygon(j)%ptr%search_flag, &
           !               polygons(i)%ptr%next_polygon(j)%ptr%index, &
           !               polygons(i)%ptr%next_polygon(j)%ptr%num_of_next
          !! if (polygons(i)%ptr%next_polygon(j)%ptr%search_flag==.false.) then
          if (.not. polygons(i)%ptr%next_polygon(j)%ptr%search_flag) then
            next_p_num=next_p_num+1
            next_polygons(next_p_num) = polygons(i)%ptr%next_polygon(j)
            polygons(i)%ptr%next_polygon(j)%ptr%search_flag = .true.
          end if
        end if
      end do
  end do

  call search_polygon_from_point(mx, my, next_p_num, next_polygons, target_polygon, level+1, &
                                 search_count, max_level, is_lat_lon)

  do i = 1, p_num
    polygons(i)%ptr%search_flag = .false.
  end do
  do i = 1, next_p_num
    next_polygons(i)%ptr%search_flag = .false.
  end do

  deallocate(next_polygons)
  
end subroutine search_polygon_from_point

!=======+=========+=========+=========+=========+=========+=========+=========+

logical function is_in_this_polygon(x, y, polygon, is_lat_lon)
  use jcf_sphere_lib, only : is_inner
  implicit none
  real(kind=8), intent(IN) :: x, y
  type(polygon_type), intent(IN) :: polygon
  logical, intent(IN) :: is_lat_lon
  real(kind=8) :: x0, y0, x1, y1, x2, y2
  integer :: i

  is_in_this_polygon = .false.
  
  if (is_lat_lon) then
    !write(0,*) "is_in_this_polygon ", x, polygon%max_xp, polygon%min_xp
    !write(0,*) "                   ", y, polygon%max_yp, polygon%min_yp
    !if ((x < polygon%min_xp).or.(x > polygon%max_xp)) return
    if ((polygon%max_xp-polygon%min_xp) < mod(x-polygon%min_xp, 360.d0)) return
    if ((y < polygon%min_yp).or.(y > polygon%max_yp)) return
    is_in_this_polygon = .true.
    return
  end if
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! min max is required !!!!!!!!!!!!!!


  x0 = polygon%point(1)%ptr%x
  y0 = polygon%point(1)%ptr%y

  !write(0,*) "   is_in_this_polygon ", polygon%num_of_point-2
  !write(0,*) "       x,  y : ",x,y
  !write(0,*) "       x0, y0 : ",x0,y0

  do i = 1, polygon%num_of_point-2
    x1 = polygon%point(i+1)%ptr%x
    y1 = polygon%point(i+1)%ptr%y
    x2 = polygon%point(i+2)%ptr%x
    y2 = polygon%point(i+2)%ptr%y
  !write(0,*) "       x1, y1 : ",x1,y1
  !write(0,*) "       x2, y2 : ",x2,y2
    if (is_inner(x0,y0,x1,y1,x2,y2,x,y)) then
      is_in_this_polygon = .true.
      return
    end if
  end do

end function is_in_this_polygon

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine add_target_polygon(my_polygon, target_polygon)
  implicit none
  type(polygon_type), pointer :: my_polygon
  type(polygon_type), pointer :: target_polygon
  type(target_polygon_type), pointer :: current_target

  !write(0,*) "add target polygon ", my_polygon%index, target_polygon%index

  ! check same polygon
  if (associated(my_polygon%target_polygon)) then
    current_target => my_polygon%target_polygon
    do while(associated(current_target))
      if (current_target%target%index == target_polygon%index) return ! same polygon
      current_target => current_target%next
    end do
  end if


  my_polygon%num_of_target = my_polygon%num_of_target + 1
  !write(0,*) "inclement num_of_target ", my_polygon%num_of_target

  if (.not.associated(my_polygon%target_polygon)) then
    allocate(my_polygon%target_polygon)
    my_polygon%target_polygon%target => target_polygon 
    nullify(my_polygon%target_polygon%before)
    nullify(my_polygon%target_polygon%next) 
  else
    current_target => my_polygon%target_polygon
    do while(associated(current_target%next)) 
      current_target => current_target%next
    end do
    allocate(current_target%next)
    current_target%next%before => current_target
    current_target => current_target%next
    current_target%target => target_polygon
    nullify(current_target%next)
  end if

end subroutine add_target_polygon

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get pointer of target polygon
function get_target_polygon_by_num(my_polygon, target_polygon_num) result(polygon)
  implicit none
  type(polygon_type), pointer :: my_polygon !< my poygon
  integer, intent(IN) :: target_polygon_num !< target polygon number

  type(polygon_type), pointer :: polygon
  type(target_polygon_type), pointer :: current_target
  integer :: counter
  nullify(polygon)

  if (my_polygon%num_of_target < target_polygon_num) return

  counter = 0 
  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    counter = counter + 1
    if (counter == target_polygon_num) then
      polygon => current_target%target
      return
    end if
    current_target => current_target%next
  end do
  
end function get_target_polygon_by_num

!=======+=========+=========+=========+=========+=========+=========+=========+
!> set coefficient of target polygon
subroutine set_target_polygon_coef1_by_num(my_polygon, target_polygon_num, coef) 
  implicit none
  type(polygon_type), pointer :: my_polygon !< my polygon
  integer, intent(IN) :: target_polygon_num !< target polygon number
  real(kind=8), intent(IN) :: coef          !< coefficient of target polygon

  type(target_polygon_type), pointer :: current_target
  integer :: counter

  if (my_polygon%num_of_target < target_polygon_num) return

  counter = 0 
  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    counter = counter + 1
    if (counter == target_polygon_num) then
      current_target%coef1 = coef
      return
    end if
    current_target => current_target%next
  end do
  
end subroutine set_target_polygon_coef1_by_num

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get coefficient of target polygon
real(kind=8) function get_target_polygon_coef1_by_num(my_polygon, target_polygon_num) 
  implicit none
  type(polygon_type), pointer :: my_polygon !< my polygon
  integer, intent(IN) :: target_polygon_num !< target polygon number

  type(target_polygon_type), pointer :: current_target
  integer :: counter

  if (my_polygon%num_of_target < target_polygon_num) return

  counter = 0 
  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    counter = counter + 1
    if (counter == target_polygon_num) then
      get_target_polygon_coef1_by_num = current_target%coef1
      return
    end if
    current_target => current_target%next
  end do
  
end function get_target_polygon_coef1_by_num

!=======+=========+=========+=========+=========+=========+=========+=========+
!> set coefficient of target polygon
subroutine set_target_polygon_coef1_by_index(my_polygon, target_polygon_index, coef) 
  implicit none
  type(polygon_type), pointer :: my_polygon   !< my polygon
  integer, intent(IN) :: target_polygon_index !< index of target polygon
  real(kind=8), intent(IN) :: coef            !< coefficient of target polygon

  type(target_polygon_type), pointer :: current_target

  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    if (current_target%target%index == target_polygon_index) then
      current_target%coef1 = coef
      return
    end if
    current_target => current_target%next
  end do
  
end subroutine set_target_polygon_coef1_by_index

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get coefficient of target polygon
real(kind=8) function get_target_polygon_coef1_by_index(my_polygon, target_polygon_index) 
  implicit none
  type(polygon_type), pointer :: my_polygon   !< my polygon
  integer, intent(IN) :: target_polygon_index !< index of target polygon

  type(target_polygon_type), pointer :: current_target

  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    if (current_target%target%index == target_polygon_index) then
      get_target_polygon_coef1_by_index = current_target%coef1
      return
    end if
    current_target => current_target%next
  end do
  
end function get_target_polygon_coef1_by_index

!=======+=========+=========+=========+=========+=========+=========+=========+
!> set coefficient of target polygon
subroutine set_target_polygon_coef2_by_num(my_polygon, target_polygon_num, coef) 
  implicit none
  type(polygon_type), pointer :: my_polygon !< my polygon
  integer, intent(IN) :: target_polygon_num !< target polygon number
  real(kind=8), intent(IN) :: coef          !< coefficient of target polygon

  type(target_polygon_type), pointer :: current_target
  integer :: counter

  if (my_polygon%num_of_target < target_polygon_num) return

  counter = 0 
  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    counter = counter + 1
    if (counter == target_polygon_num) then
      current_target%coef2 = coef
      return
    end if
    current_target => current_target%next
  end do
  
end subroutine set_target_polygon_coef2_by_num

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get coefficient of target polygon
real(kind=8) function get_target_polygon_coef2_by_num(my_polygon, target_polygon_num) 
  implicit none
  type(polygon_type), pointer :: my_polygon !< my polygon
  integer, intent(IN) :: target_polygon_num !< target polygon number

  type(target_polygon_type), pointer :: current_target
  integer :: counter

  if (my_polygon%num_of_target < target_polygon_num) return

  counter = 0 
  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    counter = counter + 1
    if (counter == target_polygon_num) then
      get_target_polygon_coef2_by_num = current_target%coef2
      return
    end if
    current_target => current_target%next
  end do
  
end function get_target_polygon_coef2_by_num

!=======+=========+=========+=========+=========+=========+=========+=========+
!> set coefficient of target polygon
subroutine set_target_polygon_coef2_by_index(my_polygon, target_polygon_index, coef) 
  implicit none
  type(polygon_type), pointer :: my_polygon   !< my polygon
  integer, intent(IN) :: target_polygon_index !< index of target polygon
  real(kind=8), intent(IN) :: coef            !< coefficient of target polygon

  type(target_polygon_type), pointer :: current_target

  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    if (current_target%target%index == target_polygon_index) then
      current_target%coef2 = coef
      return
    end if
    current_target => current_target%next
  end do
  
end subroutine set_target_polygon_coef2_by_index

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get coefficient of target polygon
real(kind=8) function get_target_polygon_coef2_by_index(my_polygon, target_polygon_index) 
  implicit none
  type(polygon_type), pointer :: my_polygon   !< my polygon
  integer, intent(IN) :: target_polygon_index !< index of target polygon

  type(target_polygon_type), pointer :: current_target

  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    if (current_target%target%index == target_polygon_index) then
      get_target_polygon_coef2_by_index = current_target%coef2
      return
    end if
    current_target => current_target%next
  end do
  
end function get_target_polygon_coef2_by_index

!=======+=========+=========+=========+=========+=========+=========+=========+
!> get pointer of target polygon
function get_target_polygon_by_index(my_polygon, target_polygon_index) result(polygon)
  implicit none
  type(polygon_type), pointer :: my_polygon   !< my polygon
  integer, intent(IN) :: target_polygon_index !< index of target polygon

  type(polygon_type), pointer :: polygon
  type(target_polygon_type), pointer :: current_target

  nullify(polygon)
  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    if (current_target%target%index == target_polygon_index) then
      polygon => current_target%target
      return
    end if
    current_target => current_target%next
  end do
  
end function get_target_polygon_by_index


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine check_mask_correspondence(this)
  implicit none
  type(mesh_type), intent(IN) :: this  
  integer :: i

  do i = 1, this%num_of_polygon
    if (this%polygon(i)%mask) then ! recv polygon is valid
      if (.not.check_polygon_mask(this%polygon(i))) then
        write(0,*) "mask correspondence check error "
        stop        
      end if 
    end if
  end do

end subroutine check_mask_correspondence

!=======+=========+=========+=========+=========+=========+=========+=========+

logical function check_polygon_mask(my_polygon)
  implicit none
  type(polygon_type), intent(IN) :: my_polygon
  type(target_polygon_type), pointer :: current_target

  current_target => my_polygon%target_polygon
  do while(associated(current_target))
    if (current_target%target%mask) then
      check_polygon_mask = .true.
      return     
    end if
  end do

  check_polygon_mask = .false.

end function check_polygon_mask

!=======+=========+=========+=========+=========+=========+=========+=========+
!> output polygon information to STDERR
subroutine write_polygon_info(my_polygon)
  implicit none
  type(polygon_type), pointer :: my_polygon     !< my polygon

  type(polygon_type), pointer :: target_polygon
  integer :: p
  integer :: ii, jj

  write(0,*)
  write(0,*) "number of my polygon target ", my_polygon%num_of_target
  do p = 1, my_polygon%num_of_target
    target_polygon => get_target_polygon_by_num(my_polygon, p)
    ii = mod(target_polygon%index-1,360)+1
    jj = int(target_polygon%index/360) + 1

    write(0,*) "my target polygon ", p, ii,jj, target_polygon%index, target_polygon%is_included
  end do

end subroutine write_polygon_info

!=======+=========+=========+=========+=========+=========+=========+=========+
!> wite performance information to file FID
subroutine write_monitor_info(FID)
  implicit none
  integer, intent(IN) :: FID    !< file id

  write(FID,*) "     monitor info " 
  write(FID,'(A,I15)') "----- number of search   : ",monitor%search_counter
  write(FID,'(A,I15)') "----- first search count : ", monitor%count_first
  write(FID,'(A,I15)') "----- max search level   : ", monitor%level_max
  write(FID,'(A,I15)') "----- max search count   : ", monitor%count_max
  write(FID,'(A,I15)') "----- max polygon index  : ", monitor%max_count_index
  write(FID,'(A,F15.8)') "----- mean search count  : ", monitor%count_mean/monitor%search_counter

end subroutine write_monitor_info

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jcf_mesh_base
