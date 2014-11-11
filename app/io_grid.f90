!> Read and maintain IO grid.
!!
!! This module manages/maintains type(mesh_type) instance for IO component of Coupled-NICAM.
!!
!! IO grid component uses regular spacing lat-lon coordinates.
!!
!!
!! 
!! io mesh definition
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
!!
!!
!! - dp : data point
!! - vpX: volume point number X 
!! - npX: next(neighbour) polygon X's data point
!! .
!!
!! \author
!!   - Arakawa T.
!------------------------------------------------------------------------------------

module io_grid
  use jcf_mesh_base, only : mesh_type
  implicit none
  private

  public :: io

  public :: init_grid !< initialize io grid
  public :: read_grid

  public :: get_data ! real(kind=8) function (NX, NY)
  public :: get_data_grid_x ! real(kind=8) function (NX, NY)
  public :: get_data_grid_y ! real(kind=8) function (NX, NY)
  public :: get_volume_grid_x ! real(kind=8) function (NX+1, NY+1)
  public :: get_volume_grid_y ! real(kind=8) function (NX+1, NY+1)

  public :: get_grid_size ! subroutine (NNX, NNY, XDIV, YDIV)

  public :: get_polygon_index ! integer funciton (i, j)
  public :: get_polygon_tuple 

  public :: get_polygon_from_idx

  public :: cal_lat_lon_to_index ! subroutine(nicam_lat, nicam_lon, lat_index, lon_index)

  integer, parameter :: FNLEN = 1024

  type(mesh_type) :: io !< type(mesh_type) instance for IO grid.

  integer :: NX  !< number of longitudinal grid points
  integer :: NY  !< number of laitudinal grid points

  real(kind=8), pointer :: xt(:,:) !< data grid x (NX, NY)
  real(kind=8), pointer :: yt(:,:) !< data grid y (NX, NY)

  real(kind=8), pointer :: vx(:,:) !< volume grid x (NX+1, NY+1)
  real(kind=8), pointer :: vy(:,:) !< volume grid y (NX+1, NY+1)


contains

  !===================================================================================
  !> initialize io grid
  subroutine init_grid(cnf_file, cnf_fid)
    use jcf_spherical_lib, only : init_spherical_lib
    use jcf_mesh_base, only : init_mesh, init_polygon, set_latlon
    use jcf_misc, only: &
         & jcf_avail_fid, &
         & LOG_fid => jcf_log_fid
    implicit none
    character(len=*), intent(IN), optional :: cnf_file !< namelist file name.
    integer         , intent(IN), optional :: cnf_fid  !< namelist file lun.

    character(len=FNLEN) :: llmapfile

    integer :: i, j
    integer :: ierr

    logical :: open_cnf_myself
    integer :: fid
    character(len=FNLEN) :: fname

    !! Open cnf(namelist) file if cnf_file is given, else just rewind.
    if ( present(cnf_fid) ) then
      fid = cnf_fid
      open_cnf_myself = .false.
    else if ( present(cnf_file) ) then
      call open_cnf_file()
      open_cnf_myself = .true.
    else
      write(0,*)'Err:init_grid/io_grid: Must specify ether fid or fname.'
      call exit(1)
    end if

    call read_cnf()

    if ( open_cnf_myself ) then
      close(fid)
    end if

    call read_lat_lon_info()

    write(LOG_fid, '(A)')    "*** init_grid/io_grid:"
    write(LOG_fid, '(A,A)')  "  llmapfile:",trim(llmapfile)
    write(LOG_fid, '(A,I5)') "     NX : ", NX
    write(LOG_fid, '(A,I5)') "     NY : ", NY
    write(LOG_fid, *)

    allocate(  xt(NX, NY)    , &
         &     yt(NX, NY)    , &
         &     vx(NX+1, NY+1), &
         &     vy(NX+1, NY+1)  )

    call init_spherical_lib()
    call init_mesh(io, NX*NY, NX*NY, (NX+1)*(NY+1))

    do j = 1, NY
      do i = 1, NX
        call init_polygon(io, get_polygon_index(i,j), 4) ! rectangle initialize
        call set_latlon(io, get_polygon_index(i,j), .true.) ! set lat lon grid flag
      end do
    end do

  contains
    !!------------------------------------------------------------------------
    subroutine open_cnf_file()
      fname = cnf_file
      fid   = jcf_avail_fid()
      open(fid,             &
           file=fname,   &
           form='formatted',    &
           status='old',        &
           iostat=ierr)
      if(ierr/=0) then
        write(0,*) 'Err:init_grid/io_grid: Cannot open cnf file:'//trim(fname)
        call exit(1)
      end if
    end subroutine open_cnf_file
    !!------------------------------------------------------------------------
    subroutine read_cnf()
      implicit none

      character(128) :: llmap_dir = './'
      character(128) :: llmap_base = 'llmap'

      namelist / iogrd /&
           llmap_dir,          &
           llmap_base

      ierr = 0
      rewind(fid)
      read(fid,nml=iogrd,iostat=ierr) !! \todo error handling.
!!$    write(LOG_fid, nml=iogrd )
      if ( ierr < 0 ) then
        write(LOG_fid,*) '*** iogrd is not specified. use default.'
      elseif( ierr > 0 ) then
        write(*,      *) 'xxx Not appropriate names in namelist iogrd. STOP.'
        write(LOG_fid,*) 'xxx Not appropriate names in namelist iogrd. STOP.'
        stop
      endif

      llmapfile = trim(llmap_dir) // "/" // trim(llmap_base) // ".info"

    end subroutine read_cnf
    !!------------------------------------------------------------------------
    subroutine read_lat_lon_info()
      implicit none

      real(8),allocatable :: lon(:),lat(:)
      integer :: llmap_fid = 20141017
      integer :: ierr

      real(8) :: PI,R2D

      PI = atan(1.d0)*4.d0
      R2D = 180.d0/PI

      open(llmap_fid,file=trim(llmapfile), form='unformatted',status='old' )

      read(llmap_fid) NX
      allocate(lon(NX))
      read(llmap_fid) lon(:)
!!$  read(llmap_fid) !! skip lon
      read(llmap_fid) NY
      allocate(lat(NY))
      read(llmap_fid) lat(:)
!!$  read(llmap_fid) !! skip lat

      close(llmap_fid)

!!$  write(0,*)'dbg:read_lat_lon_info:lon(deg):'
!!$  write(0,'(8F12.5)')lon*R2D
!!$  write(0,*)'dbg:read_lat_lon_info:lat(deg):'
!!$  write(0,'(8F12.5)')lat*R2D

      deallocate(lon, lat)

    end subroutine read_lat_lon_info

  end subroutine init_grid


  !!===================================================================================
  !> read io_grid data file and setup mesh_type instance.
  subroutine read_grid()
    implicit none
    !-----------------------------------------------------------------------------

    call set_data_grid()

    call cal_volume_grid()
    call set_grid()


    return
  end subroutine read_grid

  !===================================================================================
  !> set io grid coord.
  !!
  !! For IO grid, coord of gridpoints are calc'ed here, not read from file.
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

!!$  write(0,*)'dbg:set_data_grid:xt'
!!$  write(0,'(8F12.5)')xt
!!$  write(0,*)'dbg:set_data_grid:yt'
!!$  write(0,'(8F12.5)')yt


  end subroutine set_data_grid

  !===================================================================================
  !> Calculate volume grid.
  subroutine cal_volume_grid()
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

!!$  write(0,*)'dbg:cal_volume_grid:vx'
!!$  write(0,'(8F12.5)')vx
!!$  write(0,*)'dbg:cal_volume_grid:vy'
!!$  write(0,'(8F12.5)')vy


  end subroutine cal_volume_grid


  !===================================================================================
  !> Set coord of each data/volume point of type(mesh_type) io
  subroutine set_grid()
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

  end subroutine set_grid

  !===================================================================================
  !> return size of grid in X/Y direction
  subroutine get_grid_size(nnx, nny)
    implicit none
    integer, intent(OUT) :: nnx, nny

    nnx = NX
    nny = NY

  end subroutine get_grid_size

  !===================================================================================
  !> Return serial index of (i,j)th data point.
  integer function get_data_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_index = i + NX*(j-1)

  end function get_data_index

  !===================================================================================
  !> Return serial index of (i,j)th volume point.
  integer function get_volume_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_index = i + (NX+1)*(j-1)

  end function get_volume_index

  !===================================================================================
  !> Return data value of (i,j)th polygon
  real(kind=8) function get_data(i, j)
    use jcf_mesh_base, only : get_mesh_data => get_data
    implicit none
    integer, intent(IN) :: i, j

    get_data = get_mesh_data(io, get_polygon_index(i,j))

  end function get_data

  !===================================================================================

  real(kind=8) function get_data_grid_x(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_grid_x = xt(i,j)

  end function get_data_grid_x

  !===================================================================================

  real(kind=8) function get_data_grid_y(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_grid_y = yt(i,j)

  end function get_data_grid_y

  !===================================================================================

  real(kind=8) function get_volume_grid_x(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_grid_x = vx(i,j)

  end function get_volume_grid_x

  !===================================================================================

  real(kind=8) function get_volume_grid_y(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_grid_y = vy(i,j)

  end function get_volume_grid_y

  !===================================================================================

  !===================================================================================
  !> Return serial index of (i,j)th polygon.
  integer function get_polygon_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    if ( i>NX )then
      write(0,*) "Err:get_polygon_index/io_grid: i index error:",i
      call exit(1)
    end if
    if (j>NY) then
      write(0,*) "Err:get_polygon_index/io_grid: j index error:",j
      call exit(1)
    end if

    get_polygon_index = i + NX*(j-1)

  end function get_polygon_index

  !!===================================================================================
  !> Return (i,j) tuple of idx'th polygon
  subroutine get_polygon_tuple( idx, i, j )
    implicit none
    integer, intent(IN) :: idx
    integer, intent(OUT) :: i, j

!!$d  write(0,*)'dbg:get_polygon_tuple:idx=',idx
!!$d  write(0,*)'dbg:get_polygon_tuple:NX= ',NX

    i = mod(idx-1, NX)+1
    j = (idx-i)/NX+1
  end subroutine get_polygon_tuple

  !!===================================================================================
  !> return idx'th polygon(array of coord of each vertecies).
  subroutine get_polygon_from_idx( idx, poly, ierr )
    use jcf_mesh_base, only : get_point_x, get_point_y
    integer,intent(in) :: idx
    real(kind=8),allocatable,intent(out) :: poly(:,:)
    integer,intent(out) :: ierr

    integer,parameter :: np = 4 !< number of vertices of this polygon.
    integer :: i, j, n

    call get_polygon_tuple( idx, i, j )
    if ( i > NX .or. j > NY ) then
      write(0,*)'Error:io_grid:invalid idx:',idx
      call exit(1)
    end if

    allocate( poly(2,np) )

    do n=1,np
      poly(1,n)=get_point_x(io,idx,n)
      poly(2,n)=get_point_y(io,idx,n)
    end do

    return

  end subroutine get_polygon_from_idx




  !===================================================================================

  subroutine cal_lat_lon_to_index( ii,jj,lon,lat )
    implicit none
    integer, intent(OUT) :: ii
    integer, intent(OUT) :: jj
    real(kind=8), intent(IN) :: lon
    real(kind=8), intent(IN) :: lat

    ii = mod(int(lon/360.d0*NX), NX)+1
    if (ii < 1) ii = 1
    if (ii > NX) ii = NX

    jj = int((lat+90)/180.d0*NY)+1
    if (jj < 1) jj = 1
    if (jj > NY ) jj = NY

  end subroutine cal_lat_lon_to_index

  !===================================================================================

!!$
!!$ Below are not maintained, but left as an example of use of jcf.
!!$

  !===================================================================================

  subroutine write_data_grid(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    use jcf_misc, only: jcf_avail_fid
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j
    integer :: fid

    fid = jcf_avail_fid()
    open(unit = fid, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do j = 1, NY
      do i = 1, NX
        write(fid, *) get_data_point_x(io, get_polygon_index(i,j)), get_data_point_y(io, get_polygon_index(i,j))
      end do
    end do

    close(fid)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_data_grid

  !===================================================================================

  subroutine write_data_grid_index(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    use jcf_misc, only: jcf_avail_fid
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j
    integer :: font_size, font_angle, font_no, grid_index
    real(kind=8) :: lon, lat
    character(len=256) :: data_str

    integer :: fid

    font_size = 4
    font_angle = 0
    font_no = 1

    fid = jcf_avail_fid()
    open(unit = fid, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do j = 1, NY
      do i = 1, NX
        lon = get_data_point_x(io, get_polygon_index(i,j))
        lat = get_data_point_y(io, get_polygon_index(i,j))
        grid_index = 100*i+j
        write(data_str, "(2f10.5,3I2,A,I4)") lon, lat+0.1, font_size, font_angle, font_no, " BC", i
        write(fid, *) trim(data_str)
        write(data_str, "(2f10.5,3I2,A,I4)") lon, lat, font_size, font_angle, font_no, " TC ", j
        write(fid, *) trim(data_str)
      end do
    end do

    close(fid)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_data_grid_index

  !===================================================================================

  subroutine write_volume_grid(file_name)
    use jcf_mesh_base, only : get_point_x, get_point_y
    use jcf_misc, only: jcf_avail_fid
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j
    integer :: fid

    fid = jcf_avail_fid()
    open(unit = fid, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do j = 1, NY
      do i = 1, NX
        write(fid, *) get_point_x(io, get_polygon_index(i,j),1), get_point_y(io, get_polygon_index(i,j),1)
        write(fid, *) get_point_x(io, get_polygon_index(i,j),2), get_point_y(io, get_polygon_index(i,j),2)
        write(fid, *) get_point_x(io, get_polygon_index(i,j),3), get_point_y(io, get_polygon_index(i,j),3)
        write(fid, *) ">"
      end do
    end do

    do i = 1, NX
      write(fid, *) get_point_x(io, get_polygon_index(i,NY),1), get_point_y(io, get_polygon_index(i,NY),1)
    end do
    write(fid, *) get_point_x(io, get_polygon_index(NX,NY),4), get_point_y(io, get_polygon_index(NX,NY),4)
    write(fid, *) ">"

    do j = 1, NY
      write(fid, *) get_point_x(io, get_polygon_index(NX,j),3), get_point_y(io, get_polygon_index(NX,j),3)
    end do
    write(fid, *) get_point_x(io, get_polygon_index(NX,NY),4), get_point_y(io, get_polygon_index(NX,NY),4)


    close(fid)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_volume_grid


end module io_grid
