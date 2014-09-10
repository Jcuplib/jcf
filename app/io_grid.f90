!------------------------------------------------------------------------------------
!> @author
!! Arakawa T.
!
!> @brief read io grid and output GMT file
!! 
!! io mesh definition:
!!
!!                np4
!!    
!!           vp1-------vp4
!!            |         | 
!!       np1  |    dp   |   np3
!!            |         |
!!           vp2------- vp3
!!       J
!!       |        np2
!!        -> I
!------------------------------------------------------------------------------------

module io_grid
  use jcf_mesh_base, only : mesh_type
  implicit none
  private

  public :: io

  public :: init_grid           ! initialize io grid
  public :: get_grid_size       ! subroutine (NNX, NNY, XDIV, YDIV)
  public :: get_polygon_index   ! integer funciton (i, j)
  public :: cal_lat_lon_to_index ! subroutine(nicam_lat, nicam_lon, lat_index, lon_index)

  !integer :: NX = 144 !< number of longitudinal grid points
  !integer :: NY =  72 !< number of laitudinal grid points
  integer :: NX = 2304 !< number of longitudinal grid points
  integer :: NY = 1152 !< number of laitudinal grid points
  integer :: XDIV = 1
  integer :: YDIV = 1

  integer, parameter :: FILE_UNIT = 88

  real(kind=8), pointer :: xt(:,:) !< data grid x (NX, NY)
  real(kind=8), pointer :: yt(:,:) !< data grid y (NX, NY)

  real(kind=8), pointer :: vx(:,:) !< volume grid x (NX+1, NY+1)
  real(kind=8), pointer :: vy(:,:) !< volume grid y (NX+1, NY+1)

  type(mesh_type) :: io

contains

  !===================================================================================
  !> Initialize io grid
  !!
  subroutine init_grid(cnf_file)
    use jcf_sphere_lib, only : init_sphere_lib
    use jcf_mesh_base, only : init_mesh, init_polygon, set_latlon
    implicit none
    character(len=*), intent(IN) :: cnf_file !< namelist file name.

    character(len=128) :: llmap_base
    character(len=128) :: llmap_info
    integer :: i, j

    call read_ico2ll_cnf(cnf_file, llmap_base, XDIV, YDIV)
    llmap_info = trim(llmap_base)//".info"
    call read_lat_lon_info(trim(llmap_info), NX, NY)

    allocate(xt(NX, NY), yt(NX, NY))
    allocate(vx(NX+1, NY+1), vy(NX+1, NY+1))

    call init_sphere_lib()

    call init_mesh(io, NX*NY, NX*NY, (NX+1)*(NY+1))

    do j = 1, NY
      do i = 1, NX
        call init_polygon(io, get_polygon_index(i,j), 4) ! rectangle initialize
        call set_latlon(io, get_polygon_index(i,j), .true.) ! set lat lon grid flag
      end do
    end do

    call set_data_grid()
    call cal_volume_grid()
    call set_io_grid()

  end subroutine init_grid


  !*=======+=========+=========+=========+=========+=========+=========+=========+
  !> read namelist &iogrd
  subroutine read_ico2ll_cnf(file_name, mapbase, xd, yd)
    implicit none
    character(len=*), intent(IN) :: file_name
    character(len=*), intent(OUT) :: mapbase
    integer, intent(OUT) :: xd, yd

    character(128) :: llmap_base = './llmap'
    integer :: xdiv = 1
    integer :: ydiv = 1

    integer,parameter :: ctl_fid = 130
    integer :: ierr

    namelist / iogrd /&
         llmap_base,         &
         xdiv,               & 
         ydiv


    open(CTL_FID,             &
         file=trim(file_name),   &
         form='formatted',    &
         status='old',        &
         iostat=ierr)
    if(ierr/=0) then
      write(*,*) 'Cannot open PARAMETER file!'
      stop
    end if
    rewind(ctl_fid)
    read(ctl_fid,nml=iogrd,end=999)
999 close(ctl_fid)
    write(*,nml=iogrd)

    xd = xdiv
    yd = ydiv
    mapbase = trim(llmap_base)

  end subroutine read_ico2ll_cnf


  !*=======+=========+=========+=========+=========+=========+=========+=========+
  !> Read lan_lon info from llmap.info file.
  subroutine read_lat_lon_info(file_name, imax, jmax)
    implicit none
    character(len=*), intent(IN) :: file_name !< llmap.info file name
    integer, intent(INOUT) :: imax !< imax specified in llmap.info
    integer, intent(INOUT) :: jmax !< jmax specified in llmap.info

    real(8),allocatable :: lon(:),lat(:)
    integer :: ierr

    integer,parameter :: llmap_fid = 131

    open(llmap_fid,file=trim(file_name),&
         form='unformatted',status='old' ,iostat=ierr)
    if(ierr/=0) then
      write(*,*) 'Cannot open llmap info file!'
      stop
    end if
    read(llmap_fid) imax
    allocate(lon(imax))
    read(llmap_fid) lon(:)
    read(llmap_fid) jmax
    allocate(lat(jmax))
    read(llmap_fid) lat(:)
    close(llmap_fid)

    deallocate(lon, lat)

  end subroutine read_lat_lon_info


  !===================================================================================
  !> Get grid size defined in llmap.
  subroutine get_grid_size(nnx, nny, nxdiv, nydiv)
    implicit none
    integer, intent(OUT) :: nnx, nny, nxdiv, nydiv

    nnx = NX
    nny = NY
    nxdiv = xdiv
    nydiv = ydiv

  end subroutine get_grid_size


  !===================================================================================
  !> Read io grid data file
  integer function get_data_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_index = i + NX*(j-1)

  end function get_data_index


  !===================================================================================
  !> Read io grid data file
  integer function get_volume_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_index = i + (NX+1)*(j-1)

  end function get_volume_index

  !===================================================================================
  !> Read io grid data file
  integer function get_polygon_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    if (i>NX) then
      write(9,*) "io get_polygon_index: i index error"
      stop
    end if
    if (j>NY) then
      write(9,*) "io get_polygon_index: j index error"
      stop
    end if

    get_polygon_index = i + NX*(j-1)

  end function get_polygon_index

  !===================================================================================
  !> read io grid data file
  subroutine set_data_grid()
    implicit none
    integer :: i, j
    real(kind=8) :: dx, dy

    dx = 360.d0/NX
    dy = 180.d0/NY

    do j = 1, NY
      do i = 1, NX
        xt(i,j) = dx*0.5+(i-1)*dx      
        yt(i,j) = dy*0.5+(j-1)*dy-90.d0
      end do
    end do


  end subroutine set_data_grid

  !===================================================================================
  !>
  subroutine cal_volume_grid()
    use jcf_sphere_lib, only : cal_great_circle_center_rect
    implicit none
    integer :: i,j

    do j = 2, NY
      do i = 2, NX
        vx(i,j) = (xt(i-1,j)+xt(i,j))*0.5d0
        vy(i,j) = (yt(i,j-1)+yt(i,j))*0.5d0
      end do
    end do

    do i = 1, NX+1
      vx(i,1) = vx(i,2)
      vx(i,NY+1) = vx(i,NY)
      vy(i,1) = -90.d0
      vy(i,NY+1) = 90.d0
    end do

    do j = 1, NY+1
      vx(1,j) = 0.d0
      vx(NX+1,j) = 360.d0
      vy(1,j) = vy(2,j)
      vy(NX+1,j) = vy(NY,j)
    end do

  end subroutine cal_volume_grid

  !===================================================================================
  !>
  subroutine set_io_grid()
    use jcf_mesh_base, only : set_data_point_location, set_volume_point_location, &
         set_data_point, set_volume_point, set_next_polygon
    implicit none
    integer :: polygon_index
    integer :: i, j

    do j = 1, NY
      do i = 1, NX
        call set_data_point_location(io, get_data_index(i,j), xt(i,j), yt(i,j))
      end do
    end do

    do j = 1, NY+1
      do i = 1, NX+1
        if (j<=NY) then
          if ((vy(i,j) == vy(i,j+1))) then
            write(0,*) "vy setting error"
            write(0,*) i, j, vy(i,j)
            stop
          end if
        end if
        call set_volume_point_location(io, get_volume_index(i,j), vx(i,j), vy(i,j))
      end do
    end do

    do j = 1, NY
      do i = 1, NX
        polygon_index = get_polygon_index(i,j)
        call set_data_point(io, polygon_index, get_data_index(i,j), .true.)
        call set_volume_point(io, polygon_index, get_volume_index(i  ,j+1))
        call set_volume_point(io, polygon_index, get_volume_index(i  ,j  ))
        call set_volume_point(io, polygon_index, get_volume_index(i+1,j  ))
        call set_volume_point(io, polygon_index, get_volume_index(i+1,j+1))
      end do
    end do

    do j = 2, NY-1
      do i = 2, NX-1
        polygon_index = get_polygon_index(i,j)
        call set_next_polygon(io, polygon_index, get_polygon_index(i-1,j  ))
        call set_next_polygon(io, polygon_index, get_polygon_index(i  ,j-1))
        call set_next_polygon(io, polygon_index, get_polygon_index(i+1,j  ))
        call set_next_polygon(io, polygon_index, get_polygon_index(i  ,j+1))
      end do
    end do

    ! SW corner
    polygon_index = get_polygon_index(1,1)
    call set_next_polygon(io, polygon_index, get_polygon_index(NX,1))
    call set_next_polygon(io, polygon_index, get_polygon_index(2,1))
    call set_next_polygon(io, polygon_index, get_polygon_index(1,2))

    ! South side
    do i = 2, NX-1
      polygon_index = get_polygon_index(i,1)
      call set_next_polygon(io, polygon_index, get_polygon_index(i-1,1))
      call set_next_polygon(io, polygon_index, get_polygon_index(i+1,1))
      call set_next_polygon(io, polygon_index, get_polygon_index(i,2))
    end do

    ! SE corner
    polygon_index = get_polygon_index(NX,1)
    call set_next_polygon(io, polygon_index, get_polygon_index(NX-1,1))
    call set_next_polygon(io, polygon_index, get_polygon_index(1,1))
    call set_next_polygon(io, polygon_index, get_polygon_index(NX,2))

    ! West side
    do j = 2, NY-1
      polygon_index = get_polygon_index(1,j)
      call set_next_polygon(io, polygon_index, get_polygon_index(NX,j))
      call set_next_polygon(io, polygon_index, get_polygon_index(1,j-1))
      call set_next_polygon(io, polygon_index, get_polygon_index(2,j))
      call set_next_polygon(io, polygon_index, get_polygon_index(1,j+1))
    end do

    ! East side
    do j = 2, NY-1
      polygon_index = get_polygon_index(NX,j)
      call set_next_polygon(io, polygon_index, get_polygon_index(NX-1,j))
      call set_next_polygon(io, polygon_index, get_polygon_index(NX,j-1))
      call set_next_polygon(io, polygon_index, get_polygon_index(1,j))
      call set_next_polygon(io, polygon_index, get_polygon_index(NX,j+1))
    end do

    ! NW corner
    polygon_index = get_polygon_index(1,NY)
    call set_next_polygon(io, polygon_index, get_polygon_index(NX,NY))
    call set_next_polygon(io, polygon_index, get_polygon_index(1,NY-1))
    call set_next_polygon(io, polygon_index, get_polygon_index(2,NY))

    ! NE corner
    polygon_index = get_polygon_index(NX,NY)
    call set_next_polygon(io, polygon_index, get_polygon_index(NX-1,NY))
    call set_next_polygon(io, polygon_index, get_polygon_index(NX,NY-1))
    call set_next_polygon(io, polygon_index, get_polygon_index(1,NY))

    ! North side
    do i = 2, NX-1
      polygon_index = get_polygon_index(i,NY)
      call set_next_polygon(io, polygon_index, get_polygon_index(i-1,NY))
      call set_next_polygon(io, polygon_index, get_polygon_index(i,NY-1))
      call set_next_polygon(io, polygon_index, get_polygon_index(i+1,NY))
    end do

  end subroutine set_io_grid


  !===================================================================================
  !>
  real(kind=8) function get_data(i, j)
    use jcf_mesh_base, only : get_mesh_data => get_data
    implicit none
    integer, intent(IN) :: i, j

    get_data = get_mesh_data(io, get_polygon_index(i,j))

  end function get_data


  !===================================================================================
  !>
  real(kind=8) function get_data_grid_x(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_grid_x = xt(i,j)

  end function get_data_grid_x


  !===================================================================================
  !>
  real(kind=8) function get_data_grid_y(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_grid_y = yt(i,j)

  end function get_data_grid_y


  !===================================================================================
  !>
  real(kind=8) function get_volume_grid_x(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_grid_x = vx(i,j)

  end function get_volume_grid_x


  !===================================================================================
  !>
  real(kind=8) function get_volume_grid_y(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_grid_y = vy(i,j)

  end function get_volume_grid_y


  !===================================================================================
  !>
  subroutine cal_lat_lon_to_index(lat, lon, lat_index, lon_index)
    implicit none
    real(kind=8), intent(IN) :: lat, lon
    integer, intent(OUT) :: lat_index, lon_index

    lon_index = mod(int(lon/360.d0*NX), NX)+1
    if (lon_index < 1) lon_index = 1
    if (lon_index > NX) lon_index = NX

    lat_index = int((lat+90)/180.d0*NY)+1
    if (lat_index < 1) lat_index = 1
    if (lat_index > NY ) lat_index = NY

  end subroutine cal_lat_lon_to_index

  !===================================================================================

end module io_grid
