!------------------------------------------------------------------------------------
!> @author
!> Arakawa T.
!
!> @brief
!> read io grid

module nicamio_grid
  use io_grid, only : io
  use nicam_grid, only: nicam
  implicit none
  private

  public :: init_nicamio ! subroutine ()
  public :: read_grid_info   ! subroutine ()
  public :: search_target_polygons ! subroutine ()
  public :: cal_coefficient_io_to_nicam ! subrouitne ()
  public :: cal_coefficient_nicam_to_io_area ! subrouitne ()
  public :: cal_coefficient_nicam_to_io_monte_carlo ! subrouitne ()
!!$  public :: interpolate_io_to_nicam ! subroutine ()
  public :: interpolate_nicam_to_io ! subroutine ()

  integer :: NICAM_NX = 32
  integer :: NICAM_NY = 32
  integer :: NICAM_NZ = 40
  integer :: NICAM_NL = 10
  integer :: IO_NX = 144 !180
  integer :: IO_NY =  72 !90
  integer :: IO_XDIV = 1
  integer :: IO_YDIV = 1  

contains

  !===================================================================================
  !> @breaf
  !> initialize nicam, io grid

  subroutine init_nicamio(cnf_file)
    use nicam_grid, only : nicam_init_grid => init_grid, &
         nicam_get_grid_size => get_grid_size
    use io_grid, only : io_init_grid => init_grid, &
         io_get_grid_size => get_grid_size
    implicit none
    character(len=*), intent(IN) :: cnf_file
    character(len=128) :: llmap_dir, llmap_base

    call nicam_init_grid(cnf_file)
    call nicam_get_grid_size(NICAM_NX, NICAM_NY, NICAM_NL)


    call io_init_grid(cnf_file)
    call io_get_grid_size(IO_NX, IO_NY, IO_XDIV, IO_YDIV)

  end subroutine init_nicamio

  !===================================================================================
  !> @breaf
  !> read nicam io grid information

  subroutine read_grid_info()
    use nicam_grid, only : nicam_read_grid => read_grid
    implicit none

    call nicam_read_grid()

  end subroutine read_grid_info

  !===================================================================================
  !> @breaf
  !> search target polygons from NICAM grid to IO grid

  subroutine search_target_polygons()
    implicit none

    call point_search_nicam_to_io()
    call point_search_io_to_nicam()
    call line_search_nicam_to_io()
    call line_search_io_to_nicam()

    !write(0,*) "search_target_polygon stop"
    !stop

  end subroutine search_target_polygons

  !===================================================================================

  subroutine point_search_nicam_to_io()
    use jcf_mesh_base, only : polygon_type, polygon_ptr_type, search_polygon, get_polygon_ptr, &
         write_polygon_info, reset_monitor, write_monitor_info
    use nicam_grid, only : get_polygon_index, index2ijl, NORTH_POLE, SOUTH_POLE
    use io_grid, only : cal_lat_lon_to_index
    use mod_logfile, only : LOG_FID
    implicit none
    type(polygon_ptr_type), pointer :: start_polygon(:)
    type(polygon_type), pointer :: io_polygon
    type(polygon_type), pointer :: nicam_polygon
    integer :: start_lat, start_lon, start_np, start_sp
    real(kind=8) :: nicam_lat, nicam_lon
    integer :: i, j, k

    call reset_monitor()

    allocate(start_polygon(1))

    nicam_polygon => get_polygon_ptr(nicam, get_polygon_index(1,1,1))
    nicam_lon = nicam_polygon%point(1)%ptr%x
    nicam_lat = nicam_polygon%point(1)%ptr%y

    call cal_lat_lon_to_index(nicam_lat, nicam_lon, start_lat, start_lon)

    start_polygon(1)%ptr => get_polygon_ptr(io, (start_lat-1)*IO_NX+start_lon)

    do k = 1, NICAM_NL
      do j = 1, NICAM_NY
        do i = 1, NICAM_NX

          nicam_polygon => get_polygon_ptr(nicam, get_polygon_index(i,j,k))

          call search_polygon(nicam_polygon, start_polygon)

        end do
      end do
    end do


    nicam_polygon => get_polygon_ptr(nicam, NORTH_POLE)
    nicam_lon = nicam_polygon%point(1)%ptr%x
    nicam_lat = nicam_polygon%point(1)%ptr%y

    call cal_lat_lon_to_index(nicam_lat, nicam_lon, start_lat, start_lon)

    start_polygon(1)%ptr => get_polygon_ptr(io, (start_lat-1)*IO_NX+start_lon)

    call search_polygon(nicam_polygon, start_polygon)

    nicam_polygon => get_polygon_ptr(nicam, SOUTH_POLE)
    nicam_lon = nicam_polygon%point(1)%ptr%x
    nicam_lat = nicam_polygon%point(1)%ptr%y

    call cal_lat_lon_to_index(nicam_lat, nicam_lon, start_lat, start_lon)

    start_polygon(1)%ptr => get_polygon_ptr(io, (start_lat-1)*IO_NX+start_lon)

    call search_polygon(nicam_polygon, start_polygon)



    write(LOG_FID, *) "Msg : Sub[point_search_nicam_to_io]/Mod[nicamio_grid]"
    write(LOG_FID, *) "----- point search nicam to io finish"
    call write_monitor_info(LOG_FID)
    write(LOG_FID, *)

    ! check write
    !io_polygon => get_polygon_ptr(io, 360*(154-1)+308)
    !call write_polygon_info(io_polygon)

    deallocate(start_polygon)

  end subroutine point_search_nicam_to_io

  !===================================================================================

  subroutine point_search_io_to_nicam()
    use jcf_mesh_base, only : polygon_type, polygon_ptr_type, search_polygon, get_polygon_ptr, &
         write_polygon_info, reset_monitor, write_monitor_info
    use nicam_grid, only : SOUTH_POLE
    use io_grid, only : get_polygon_index
    use mod_logfile, only : LOG_FID
    implicit none
    type(polygon_ptr_type), pointer :: start_polygon(:)
    type(polygon_type), pointer :: io_polygon
    type(polygon_type), pointer :: nicam_polygon
    integer :: i, j

    call reset_monitor()

    allocate(start_polygon(1))

    start_polygon(1)%ptr => get_polygon_ptr(nicam, SOUTH_POLE)

    do j = 1,IO_NY
      do i = 1, IO_NX

        io_polygon => get_polygon_ptr(io, get_polygon_index(i,j))

        call search_polygon(io_polygon, start_polygon)

      end do
    end do




    write(LOG_FID, *) "Msg : Sub[point_search_io_to_nicam]/Mod[nicamio_grid]"
    write(LOG_FID, *) "----- point search io to nicam finish"
    call write_monitor_info(LOG_FID)
    write(LOG_FID, *)


    ! check write
    !nicam_polygon => get_polygon_ptr(nicam, 101)
    !call write_polygon_info(nicam_polygon)

    deallocate(start_polygon)

  end subroutine point_search_io_to_nicam


  !===================================================================================

  subroutine line_search_nicam_to_io()
    use jcf_mesh_base, only : polygon_type, polygon_ptr_type, search_polygon_by_side, get_polygon_ptr, &
         write_polygon_info, reset_monitor, write_monitor_info
    use nicam_grid, only : get_polygon_index, index2ijl
    use mod_logfile, only : LOG_FID
    implicit none
    type(polygon_type), pointer :: io_polygon
    type(polygon_type), pointer :: nicam_polygon
    integer :: i, j, k

    do k = 1, NICAM_NL
      do j = 1, NICAM_NY
        do i = 1, NICAM_NX

          nicam_polygon => get_polygon_ptr(nicam, get_polygon_index(i,j,k))

          call search_polygon_by_side(nicam_polygon)

        end do
      end do
    end do


    write(LOG_FID, *) "Msg : Sub[line_search_nicam_to_io]/Mod[nicamio_grid]"
    write(LOG_FID, *) "----- line search nicam to io finish"
    write(LOG_FID, *)

  end subroutine line_search_nicam_to_io

  !===================================================================================

  subroutine line_search_io_to_nicam()
    use jcf_mesh_base, only : polygon_type, polygon_ptr_type, search_polygon_by_side, get_polygon_ptr, &
         write_polygon_info, reset_monitor, write_monitor_info
    use io_grid, only : get_polygon_index
    use mod_logfile, only : LOG_FID
    implicit none
    type(polygon_type), pointer :: io_polygon
    integer :: i, j

    do j = 1, IO_NY
      do i = 1, IO_NX

        io_polygon => get_polygon_ptr(io, get_polygon_index(i,j))
        call search_polygon_by_side(io_polygon)

      end do
    end do


    write(LOG_FID, *) "Msg : Sub[line_search_io_to_nicam]/Mod[nicamio_grid]"
    write(LOG_FID, *) "----- line search io to nicam finish"
    write(LOG_FID, *)

  end subroutine line_search_io_to_nicam

  !===================================================================================

  subroutine write_mapping_table_io_to_nicam(file_name)
    use jcf_mesh_base, only : polygon_type, polygon_ptr_type, get_polygon_ptr, &
         get_target_polygon_by_num, get_target_polygon_coef1_by_num
    use nicam_grid, only : get_polygon_index, NORTH_POLE, SOUTH_POLE
    implicit none
    character(len=*), intent(IN) :: file_name
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    integer :: i, j, k, t
    integer :: counter
    real(kind=8) :: coef
    integer, parameter :: FILE_UNIT = 111

    open(unit=FILE_UNIT, FILE=trim(file_name), access = "SEQUENTIAL")

    counter = 0
    do k = 1, NICAM_NL
      do j = 1, NICAM_NY
        do i = 1, NICAM_NX

          nicam_polygon => get_polygon_ptr(nicam, get_polygon_index(i,j,k))
          if (nicam_polygon%mask) then
            do t = 1, nicam_polygon%num_of_target
              io_polygon => get_target_polygon_by_num(nicam_polygon, t)
              if (io_polygon%mask) then
                counter = counter + 1
                coef = get_target_polygon_coef1_by_num(nicam_polygon, t)
                write(FILE_UNIT, '(I8,I8,F15.8)') nicam_polygon%index, io_polygon%index, coef
              end if
            end do
          end if

        end do
      end do
    end do

    nicam_polygon => get_polygon_ptr(nicam, NORTH_POLE)
    if (nicam_polygon%mask) then
      do t = 1, nicam_polygon%num_of_target
        io_polygon => get_target_polygon_by_num(nicam_polygon, t)
        if (io_polygon%mask) then
          counter = counter + 1
          coef = get_target_polygon_coef1_by_num(nicam_polygon, t)
          write(FILE_UNIT, '(I8,I8,F15.8)') nicam_polygon%index, io_polygon%index, coef
        end if
      end do
    end if

    nicam_polygon => get_polygon_ptr(nicam, SOUTH_POLE)
    if (nicam_polygon%mask) then
      do t = 1, nicam_polygon%num_of_target
        io_polygon => get_target_polygon_by_num(nicam_polygon, t)
        if (io_polygon%mask) then
          counter = counter + 1
          coef = get_target_polygon_coef1_by_num(nicam_polygon, t)
          write(FILE_UNIT, '(I8,I8,F15.8)') nicam_polygon%index, io_polygon%index, coef
        end if
      end do
    end if

    close(FILE_UNIT)

  end subroutine write_mapping_table_io_to_nicam

  !===================================================================================

  subroutine write_mapping_table_nicam_to_io(file_name)
    use jcf_mesh_base, only : polygon_type, polygon_ptr_type, get_polygon_ptr, &
         get_target_polygon_by_num, get_target_polygon_coef1_by_num
    use io_grid, only : get_polygon_index
    implicit none
    character(len=*), intent(IN) :: file_name
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    integer :: i, j, k, t
    integer :: counter
    real(kind=8) :: coef
    integer, parameter :: FILE_UNIT = 111

    open(unit=FILE_UNIT, FILE=trim(file_name), access = "SEQUENTIAL")

    counter = 0
    do j = 1, IO_NY
      do i = 1, IO_NX

        io_polygon => get_polygon_ptr(io, get_polygon_index(i,j))
        !if (io_polygon%mask) then
        do t = 1, io_polygon%num_of_target
          nicam_polygon => get_target_polygon_by_num(io_polygon, t)
          !if (nicam_polygon%mask) then
          counter = counter + 1
          coef = get_target_polygon_coef1_by_num(io_polygon, t)
          !write(FILE_UNIT, '(I8,I8,I8,F15.8)') counter, io_polygon%index, nicam_polygon%index, coef
          write(FILE_UNIT, '(I8,I8,F15.7)') io_polygon%index, nicam_polygon%index, coef
          !end if
        end do
        !end if
      end do
    end do

    close(FILE_UNIT)

  end subroutine write_mapping_table_nicam_to_io

  !===================================================================================

  subroutine cal_coefficient_io_to_nicam(table_file_name, coef_check_mode, coef_cal_type, side_div_num)
    use jcf_sphere_lib, only : cal_great_circle_area
    use jcf_mesh_base, only : polygon_type, get_polygon_ptr, get_target_polygon_by_num, &
         get_target_polygon_coef1_by_num
    use jcf_coef_base, only :  &
         cal_coefficient
    use nicam_grid, only : get_polygon_index, NORTH_POLE, SOUTH_POLE
    implicit none
    character(len=*), intent(IN) :: table_file_name
    integer, intent(IN) :: coef_check_mode
    character(len=*), intent(IN) :: coef_cal_type
    integer, intent(IN) :: side_div_num
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    real(kind=8) :: coef
    integer :: i, j, k
    integer :: p
    logical :: is_check_coef
    !real(kind=8) :: test_lat(3), test_lon(3)
    !real(kind=8) :: area

    !test_lat(1) = -90.d0 ; test_lon(1) = 0.d0
    !test_lat(2) = 0.d0 ; test_lon(2) = 0.d0
    !test_lat(3) = 0.d0 ; test_lon(3) = 90.d0

    !area = cal_great_circle_area(3, test_lat, test_lon)
    !write(0,*) "interpolate_io_to_nicam, test area = ", area

    is_check_coef = .false.
    if (coef_check_mode == 1) is_check_coef = .true.

    write(0,*) "coef cal start"
    do k = 1, NICAM_NL
      do j = 1, NICAM_NY
        do i = 1, NICAM_NX
          nicam_polygon => get_polygon_ptr(nicam, get_polygon_index(i,j,k))
          call cal_coefficient(nicam_polygon, is_check_coef)
          !call cal_coefficient_without_mask(nicam_polygon)
        end do
      end do
    end do

    nicam_polygon => get_polygon_ptr(nicam, NORTH_POLE)
    call cal_coefficient(nicam_polygon, is_check_coef)
    nicam_polygon => get_polygon_ptr(nicam, SOUTH_POLE)
    call cal_coefficient(nicam_polygon, is_checK_coef)


    write(0,*) "coef cal end"

    call write_mapping_table_io_to_nicam(trim(table_file_name))

  end subroutine cal_coefficient_io_to_nicam

  !===================================================================================

  subroutine interpolate_io_to_nicam()
    use jcf_sphere_lib, only : cal_great_circle_area
    use jcf_mesh_base, only : polygon_type, get_polygon_ptr, get_target_polygon_by_num, &
         get_target_polygon_coef1_by_num
    use nicam_grid, only : get_polygon_index
    implicit none
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    real(kind=8) :: coef
    integer :: i, j, k
    integer :: p

    do k = 1, NICAM_NL
      do j = 1, NICAM_NY
        do i = 1, NICAM_NX
          nicam_polygon => get_polygon_ptr(nicam, get_polygon_index(i,j,k))
          nicam_polygon%data = 0.d0
          do p = 1, nicam_polygon%num_of_target
            io_polygon => get_target_polygon_by_num(nicam_polygon, p)          
            coef = get_target_polygon_coef1_by_num(nicam_polygon, p)
            nicam_polygon%data = nicam_polygon%data + io_polygon%data*coef
          end do
        end do
      end do
    end do

  end subroutine interpolate_io_to_nicam

  !===================================================================================

  subroutine cal_coefficient_nicam_to_io_area(table_file_name, coef_check_mode)
    use jcf_sphere_lib, only : cal_great_circle_area
    use jcf_mesh_base, only : polygon_type, get_polygon_ptr, get_target_polygon_by_num, &
         get_target_polygon_coef1_by_num
    use jcf_coef_base, only :  cal_coefficient
    use io_grid, only : get_polygon_index
    use mod_logfile, only : LOG_FID
    implicit none
    character(len=*), intent(IN) :: table_file_name
    integer, intent(IN) :: coef_check_mode
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    integer :: i, j
    logical :: is_check_coef

    is_check_coef = .false.
    if (coef_check_mode == 1) is_check_coef = .true.

    write(0,*) "coef cal start"
    do j = 1, IO_NY
      do i = 1, IO_NX
        io_polygon => get_polygon_ptr(io, get_polygon_index(i,j))
        call cal_coefficient(io_polygon, is_check_coef)
      end do
    end do

    write(0,*) "coef cal end"

    call write_mapping_table_nicam_to_io(trim(table_file_name))

    write(LOG_FID, *) "Msg : Sub[cal_coefficient_nicam_to_io_area]/Mod[nicamio_grid]"
    write(LOG_FID, *) "----- area ratio coefficient calculation nicam to io finish"
    write(LOG_FID, *) "----- mapping table file name : "//trim(table_file_name)
    write(LOG_FID, *)

  end subroutine cal_coefficient_nicam_to_io_area

  !===================================================================================

  subroutine cal_coefficient_nicam_to_io_monte_carlo&
       &(table_file_name, coef_check_mode, side_div_num)
    use jcf_sphere_lib, only : cal_great_circle_area
    use jcf_mesh_base, only : polygon_type, get_polygon_ptr, get_target_polygon_by_num, &
         get_target_polygon_coef1_by_num
    use mod_coef, only :  cal_monte_carlo_points, cal_monte_carlo_coef
    use nicam_grid, only : get_polygon_index_nicam => get_polygon_index, NORTH_POLE, SOUTH_POLE
    use io_grid, only : get_polygon_index_io => get_polygon_index
    use mod_logfile, only : LOG_FID
    implicit none
    character(len=*), intent(IN) :: table_file_name
    integer, intent(IN) :: coef_check_mode
    integer, intent(IN) :: side_div_num
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    integer :: i, j, k
    logical :: is_check_coef

    write(LOG_FID, *) "Msg : Sub[cal_coefficient_nicam_to_io_monte_carlo]/Mod[nicoco_grid]"

    is_check_coef = .false.
    if (coef_check_mode == 1) is_check_coef = .true.

    write(0,*) "coef cal start"
    write(0,*) "calculation type : monte carlo"
    write(0,*) "number of inner points : ", (2**side_div_num)*(2**side_div_num-1)/2

    write(LOG_FID, *) "----- monte carlo point calculation start"
    write(LOG_FID, *) "----- number of inner points : ", (2**side_div_num)*(2**side_div_num-1)/2

    nicam_polygon => get_polygon_ptr(nicam, SOUTH_POLE)
    call cal_monte_carlo_points(nicam_polygon, side_div_num)
    nicam_polygon => get_polygon_ptr(nicam, NORTH_POLE)
    call cal_monte_carlo_points(nicam_polygon, side_div_num)

    do k = 1, NICAM_NL
      do j = 1, NICAM_NY
        do i = 1, NICAM_NX
          nicam_polygon => get_polygon_ptr(nicam, get_polygon_index_nicam(i,j,k))
          call cal_monte_carlo_points(nicam_polygon, side_div_num)
        end do
      end do
    end do

    write(LOG_FID, *) "----- monte carlo point calculation finish"
    write(LOG_FID, *) "----- coefficient calculation start"


    write(0,*) "coef cal second step"

    do j = 1, IO_NY
      do i = 1, IO_NX
        io_polygon => get_polygon_ptr(io, get_polygon_index_io(i,j))
        call cal_monte_carlo_coef(io_polygon)
      end do
    end do

    write(0,*) "coef cal end"

    call write_mapping_table_nicam_to_io(trim(table_file_name))

    write(LOG_FID, *) "----- coefficient calculation finish"
    write(LOG_FID, *) "----- mapping table file name : "//trim(table_file_name)
    write(LOG_FID, *)

  end subroutine cal_coefficient_nicam_to_io_monte_carlo

  !===================================================================================

  subroutine interpolate_nicam_to_io()
    use jcf_sphere_lib, only : cal_great_circle_area
    use jcf_mesh_base, only : polygon_type, get_polygon_ptr, get_target_polygon_by_num, &
         get_target_polygon_coef1_by_num
    use io_grid, only : get_polygon_index
    implicit none
    type(polygon_type), pointer :: nicam_polygon
    type(polygon_type), pointer :: io_polygon
    real(kind=8) :: coef
    integer :: i, j, k
    integer :: p

    do j = 1, IO_NY
      do i = 1, IO_NX
        io_polygon => get_polygon_ptr(io, get_polygon_index(i,j))
        io_polygon%data = 0.d0
        do p = 1, io_polygon%num_of_target
          nicam_polygon => get_target_polygon_by_num(io_polygon, p)          
          coef = get_target_polygon_coef1_by_num(io_polygon, p)
          io_polygon%data = io_polygon%data + nicam_polygon%data*coef
        end do
      end do
    end do

  end subroutine interpolate_nicam_to_io

  !===================================================================================

end module nicamio_grid

