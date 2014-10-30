!> Read COCO grid and manage type(mesh_type) instance.
!!
!! This module manages/maintains type(mesh_type) instance for COCO
!! ''Tri-polar'' grid.
!!
!! Grid points their latitude north than `TLAT` are assumed to be
!! oblique, and south than this are assumed to be usual lon-lat
!! coordinate (but NOT limited to regular spacing).
!!
!! \note
!! Grid point coordinates are read from some format file. If you
!! change it, modify read_grid().
!!
!!
!! COCO mesh definition is assumed as follows;
!!
!!
!!              np4
!!       
!!          vp1-------vp4
!!          |         | 
!!     np1  |    dp   |   np3
!!          |         |
!!          vp2------- vp3
!!      J
!!      |        np2
!!       -> I
!!    
!! - dp : data point
!! - vpX: volume point number X 
!! - npX: next(neighbour) polygon X's data point
!! .
!!
!! \author
!!   - Arakawa T.
!------------------------------------------------------------------------------------

module coco_grid
  use jcf_mesh_base, only : mesh_type
  implicit none
  private

  public :: coco

  public :: init_grid
  public :: read_grid

  public :: get_data
  public :: get_mask
  public :: get_read_mask
  public :: get_read_mask_value
  public :: get_data_grid_x
  public :: get_data_grid_y
  public :: get_volume_grid_x
  public :: get_volume_grid_y

  public :: get_grid_size

  public :: get_polygon_index
  public :: get_polygon_tuple 

  public :: get_polygon_from_idx

  public :: get_rang
  public :: get_tripolar_lat
  public :: cal_first_guess_tuple


  integer, parameter :: STR_LEN = 1024
  integer, parameter :: FILE_UNIT = 887

  type(mesh_type) :: coco !< type(mesh_type) instance for COCO grid.

  integer :: NX                 !< number of longitudinal grid points
  integer :: NY                 !< number of laitudinal grid points

  real(kind=8) :: TLAT          !< tripole latitude(degN)


  !! contents of coco_grid_file

  real(kind=8), pointer :: xt(:,:)   !< data grid x (NX, NY)
  real(kind=8), pointer :: yt(:,:)   !< data grid y (NX, NY)
  real(kind=8), pointer :: rang(:,:) !< oblique (NX, NY)
  real(kind=8), pointer :: mask(:,:) !< land mask from COCO. (NX, NY)

  real(kind=8), pointer :: vx(:,:) !< volume grid x (NX+1, NY+1)
  real(kind=8), pointer :: vy(:,:) !< volume grid y (NX+1, NY+1)


  character(len=STR_LEN) :: coco_grid_file='GRIDO.m52' !< COCO grid file name
  integer                :: coco_nx=360                !< NX
  integer                :: coco_ny=256                !< NY
  real(kind=8)           :: tripole_lat = 62.73081     !< TLAT
  namelist / cocogrd /  &
       & coco_grid_file,&
       & coco_nx,       &
       & coco_ny,       &
       & tripole_lat

  real(kind=8), parameter,private :: EPSILON=1.d-10

contains

  !!===================================================================================
  !> Initialize type(mesh_type) instance for coco grid.
  !!
  !! If cnf_fid is given, assumes that cnf file is already opened
  !! somewhere.
  !! 
  !! If cnf_file is given, assumes that cnf file is NOT opened yet,
  !! open/read/close this file here.
  !!
  !! None of cnf_fid/cnf_file is given, cause error.
  !! 
  subroutine init_grid(cnf_file, cnf_fid)
    use jcf_spherical_lib, only : &
         & init_spherical_lib
    use jcf_mesh_base, only : &
         & init_mesh,         &
         & init_polygon
    use jcf_misc, only: &
         & jcf_avail_fid, &
         & LOG_FID => jcf_log_fid
    implicit none
    character(len=*), intent(IN), optional :: cnf_file !< namelist file name.
    integer         , intent(IN), optional :: cnf_fid  !< namelist file lun.

    integer :: fid
    character(len=STR_LEN) :: fname

    logical :: open_cnf_myself
    integer :: i, j
    integer :: ierr


    !! Open cnf(namelist) file if cnf_file is given, else just rewind.
    if ( present(cnf_fid) ) then
      fid = cnf_fid
      open_cnf_myself = .false.
    else if ( present(cnf_file) ) then
      call open_cnf_file()
      open_cnf_myself = .true.
    else
      write(0,*)'Err:init_grid/coco_grid: Must specify ether fid or fname.'
      call exit(1)
    end if

    rewind(fid)
    read(fid, nml=cocogrd,iostat=ierr)
    if ( ierr < 0 ) then
      write(LOG_FID,*) '*** cocogrd is not specified. use default.'
    elseif( ierr > 0 ) then
      write(*,      *) 'xxx Not appropriate names in namelist cocogrd. STOP.'
      write(LOG_FID,*) 'xxx Not appropriate names in namelist cocogrd. STOP.'
      stop
    endif
    write(LOG_FID, nml=cocogrd )
    if ( open_cnf_myself ) then
      close(fid)
    end if

    NX = coco_nx
    NY = coco_ny
    TLAT = tripole_lat+EPSILON

    write(LOG_FID, '(A)')       "Msg : Sub[init_grid]/Mod[coco_grid]"
    write(LOG_FID, '(A)')       "----- coco grid initialize "
    write(LOG_FID, '(A,I5)')    " ----- NX : ", NX
    write(LOG_FID, '(A,I5)')    " ----- NY : ", NY
    write(LOG_FID, '(A,F10.5)') " ----- TriLat : ", TLAT
    write(LOG_FID, *)

    allocate( xt  (NX,   NY)  ,&
         &    yt  (NX,   NY)  ,&
         &    rang(NX,   NY)  ,&
         &    mask(NX,   NY)  ,&
         &    vx  (NX+1, NY+1),&
         &    vy  (NX+1, NY+1) )

    call init_spherical_lib()
    call init_mesh(coco, NX*NY, NX*NY, (NX+1)*(NY+1))

    do j = 1, NY
      do i = 1, NX
        call init_polygon(coco, get_polygon_index(i,j), 4) ! rectangle initialize
      end do
    end do

    return

  contains
    subroutine open_cnf_file()
      fname = cnf_file
      fid   = jcf_avail_fid()
      open(fid,             &
           file=fname,   &
           form='formatted',    &
           status='old',        &
           iostat=ierr)
      if(ierr/=0) then
        write(0,*) 'Err:init_grid/coco_grid: Cannot open cnf file:'//trim(fname)
        call exit(1)
      end if
    end subroutine open_cnf_file

  end subroutine init_grid


  !!===================================================================================
  !> read coco_grid data file and setup mesh_type instance.
  subroutine read_grid()
    implicit none

    call input_grid( trim(coco_grid_file) )

    call cal_volume_grid
    call set_grid

    return
  end subroutine read_grid


  !!===================================================================================
  !> read coco_grid data file.
  subroutine input_grid(fname)
    use jcf_misc, only: &
         & jcf_avail_fid,&
         & LOG_FID => jcf_log_fid
    implicit none
    character(len=*), intent(in) :: fname

    integer :: fid, ierr
    integer :: l

    fid = jcf_avail_fid()
    open(  unit   = fid,           &
         & file   = trim(fname),   &
         & access = "SEQUENTIAL",  &
         & form   = "UNFORMATTED", &
         & status = "OLD",         &
         & iostat = ierr           )

    if ( ierr /= 0 ) then
      write(0,*) 'Err:input_grid/coco_grid: Cannot open coco_grid file!'//trim(fname)
      call exit(1)
    endif

    read(fid) xt
    read(fid) yt
    read(fid)   !! sx
    read(fid)   !! sy
    read(fid) mask
    read(fid) rang

    close(fid)

    write(LOG_FID,*) '*** coco_grid file was read successfully:'//trim(coco_grid_file)

    return
  end subroutine input_grid


  !!===================================================================================
  !> Calculate volume grid from data_grid
  subroutine cal_volume_grid()
    use jcf_spherical_lib, only : get_center_of_spherical_rectangle
    implicit none

    real(kind=8) :: dx(0:NX+1, 0:NY+1) , dy(0:NX+1, 0:NY+1)
    real(kind=8) :: dx1, dx2, dx3, dx4
    real(kind=8) :: dy1, dy2, dy3, dy4
    integer :: i,j

    do j = 1, NY
      do i = 1, NX
        dx(i,j) = xt(i,j)
        dy(i,j) = yt(i,j)
      end do
    end do

    do j = 1, NY
      dx(0, j) = dx(NX, j)
      dy(0, j) = dy(NX, j)
      dx(NX+1, j) = dx(1, j)
      dy(NX+1, j) = dy(1, j)
    end do

    do i = 1, NX
      dx(i, 0) = dx(i, 1)
      dy(i, 0) = -90.d0
    end do

    do i = 1, NX/2
      dx(i, NY+1) = dx(NX-i+1, NY)
      dy(i, NY+1) = dy(NX-i+1, NY)
    end do

    do i = NX/2+1, NX
      dx(i, NY+1) = dx(NX-i+1, NY)
      dy(i, NY+1) = dy(NX-i+1, NY)
    end do

    dx(0, 0) = dx(0, 1)
    dy(0, 0) = -90.d0

    dx(NX+1, 0) = dx(NX+1, 1)
    dy(NX+1, 0) = -90.d0

    dx(0, NY+1) = dx(NX, NY+1)
    dy(0, NY+1) = dy(NX, NY+1)
    dx(NX+1, NY+1) = dx(1, NY+1)
    dy(NX+1, NY+1) = dy(1, NY+1)

    do j = 1, NY+1
      do i = 1, NX+1
        dx1 = dx(i-1, j-1)
        dx2 = dx(i  , j-1)
        dx3 = dx(i-1, j)
        dx4 = dx(i  , j)
        dy1 = dy(i-1, j-1)
        dy2 = dy(i  , j-1)
        dy3 = dy(i-1, j)
        dy4 = dy(i  , j)

        if ((max(dx1, dx2, dx3, dx4)-min(dx1, dx2, dx3, dx4))>180.d0) then ! cross 0 line
          if (dx1 > 180.d0) dx1 = dx1 - 360.d0
          if (dx2 > 180.d0) dx2 = dx2 - 360.d0
          if (dx3 > 180.d0) dx3 = dx3 - 360.d0
          if (dx4 > 180.d0) dx4 = dx4 - 360.d0
        end if


        vx(i,j) = (dx1+dx2+dx3+dx4)/4.d0

        vx(i,j) = mod(vx(i,j), 360.d0)
        if (vx(i,j) < 0.d0) vx(i,j) = vx(i,j) + 360.d0

        vy(i,j) = (dy1+dy2+dy3+dy4)/4.d0

        if (vy(i,j) > TLAT ) then ! Tripole grid
          call get_center_of_spherical_rectangle(       &
               & vx(i,j), vy(i,j),                      &
               & dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4 )
        end if
      end do
    end do

    ! south pole grid
    do i = 1, NX+1
      vy(i,1) = -90.d0
    end do

    ! north pole grid
    i = NX/4+1 ; j = NY+1
    vy(i,j) = 90.d0
    i = NX/4*3+1 ; j = NY+1
    vy(i,j) = 90.d0

    return
  end subroutine cal_volume_grid


  !!===================================================================================
  !> Set coord of each data/volume point of type(mesh_type).
  subroutine set_grid()
    use jcf_mesh_base, only : &
         & set_data_point_location,&
         & set_volume_point_location,&
         & set_data_point,&
         & set_volume_point,&
         & set_next_polygon
    implicit none

    integer :: polygon_index
    integer :: i, j

    !! set coord of (i,j)th data point as (xt(i,j),yt(i,j))
    do j = 1, NY
      do i = 1, NX
        call set_data_point_location(coco, get_data_index(i,j), xt(i,j), yt(i,j))
      end do
    end do

    !! set coord of (i,j)th volume point as (vx(i,j),vy(i,j))
    !! \todo get error check out of do-loop.
    do j = 1, NY+1
      do i = 1, NX+1
        if (j<=NY) then
          if ((vy(i,j) == vy(i,j+1))) then
            write(0,*) "Err:set_grid/coco_grid: vy setting error:"
            write(0,*) i, j, vy(i,j)
            call exit(1)
          end if
        end if
        call set_volume_point_location(coco, get_volume_index(i,j), vx(i,j), vy(i,j))
      end do
    end do

    !! set (i,j)th polygon's data/volume point.
    do j = 1, NY
      do i = 1, NX
        polygon_index = get_polygon_index(i,j)
        call set_data_point(coco, polygon_index, get_data_index(i,j), (mask(i,j) == 1.d0))
        call set_volume_point(coco, polygon_index, get_volume_index(i  ,j+1))
        call set_volume_point(coco, polygon_index, get_volume_index(i  ,j  ))
        call set_volume_point(coco, polygon_index, get_volume_index(i+1,j  ))
        call set_volume_point(coco, polygon_index, get_volume_index(i+1,j+1))
      end do
    end do

    !! set (i,j)th polygon's neighbour polygon.
    do j = 2, NY-1
      do i = 2, NX-1
        polygon_index = get_polygon_index(i,j)
        call set_next_polygon(coco, polygon_index, get_polygon_index(i-1,j  ))
        call set_next_polygon(coco, polygon_index, get_polygon_index(i  ,j-1))
        call set_next_polygon(coco, polygon_index, get_polygon_index(i+1,j  ))
        call set_next_polygon(coco, polygon_index, get_polygon_index(i  ,j+1))
      end do
    end do

    ! SW corner
    polygon_index = get_polygon_index(1,1)
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX,1))
    call set_next_polygon(coco, polygon_index, get_polygon_index(2,1))
    call set_next_polygon(coco, polygon_index, get_polygon_index(1,2))

    ! South side
    do i = 2, NX-1
      polygon_index = get_polygon_index(i,1)
      call set_next_polygon(coco, polygon_index, get_polygon_index(i-1,1))
      call set_next_polygon(coco, polygon_index, get_polygon_index(i+1,1))
      call set_next_polygon(coco, polygon_index, get_polygon_index(i,2))
    end do

    ! SE corner
    polygon_index = get_polygon_index(NX,1)
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX-1,1))
    call set_next_polygon(coco, polygon_index, get_polygon_index(1,1))
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX,2))

    ! West side
    do j = 2, NY-1
      polygon_index = get_polygon_index(1,j)
      call set_next_polygon(coco, polygon_index, get_polygon_index(NX,j))
      call set_next_polygon(coco, polygon_index, get_polygon_index(1,j-1))
      call set_next_polygon(coco, polygon_index, get_polygon_index(2,j))
      call set_next_polygon(coco, polygon_index, get_polygon_index(1,j+1))
    end do

    ! NW corner
    polygon_index = get_polygon_index(1,NY)
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX,NY))
    call set_next_polygon(coco, polygon_index, get_polygon_index(1,NY-1))
    call set_next_polygon(coco, polygon_index, get_polygon_index(2,NY))
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX,NY))

    ! East side
    do j = 2, NY-1
      polygon_index = get_polygon_index(NX,j)
      call set_next_polygon(coco, polygon_index, get_polygon_index(NX-1,j))
      call set_next_polygon(coco, polygon_index, get_polygon_index(NX,j-1))
      call set_next_polygon(coco, polygon_index, get_polygon_index(1,j))
      call set_next_polygon(coco, polygon_index, get_polygon_index(NX,j+1))
    end do

    ! NE corner
    polygon_index = get_polygon_index(NX,NY)
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX-1,NY))
    call set_next_polygon(coco, polygon_index, get_polygon_index(NX,NY-1))
    call set_next_polygon(coco, polygon_index, get_polygon_index(1,NY))
    call set_next_polygon(coco, polygon_index, get_polygon_index(1,NY))

    ! North side
    do i = 2, NX-1
      polygon_index = get_polygon_index(i,NY)
      call set_next_polygon(coco, polygon_index, get_polygon_index(i-1,NY))
      call set_next_polygon(coco, polygon_index, get_polygon_index(i,NY-1))
      call set_next_polygon(coco, polygon_index, get_polygon_index(i+1,NY))
      call set_next_polygon(coco, polygon_index, get_polygon_index(NX-i+1,NY))
    end do

    return
  end subroutine set_grid


  !!===================================================================================
  !> return size of grid in X/Y direction
  subroutine get_grid_size(lon_num,lat_num)
    implicit none
    integer, intent(out) :: lon_num !< num of grid in X direction
    integer, intent(out) :: lat_num !< num of grid in Y direction

    lon_num = NX
    lat_num = NY

    return
  end subroutine get_grid_size

  !!===================================================================================
  !> return `TLAT` Tri-polar transition latitude.
  subroutine get_tripolar_lat(tri_lat)
    implicit none
    real(kind=8), intent(out) :: tri_lat !< Tri-polar transition latitude [degN]

    tri_lat = TLAT

    return
  end subroutine get_tripolar_lat

  !!===================================================================================
  !> Return serial index of (i,j)th data point.
  integer function get_data_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_index = i + NX*(j-1)

  end function get_data_index

  !!===================================================================================
  !> Return serial index of (i,j)th volume point.
  integer function get_volume_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_index = i + (NX+1)*(j-1)

  end function get_volume_index

  !!===================================================================================
  !> Return data value of (i,j)th polygon
  real(kind=8) function get_data(i, j)
    use jcf_mesh_base, only : get_mesh_data => get_data
    implicit none
    integer, intent(IN) :: i, j

    get_data = get_mesh_data(coco, get_polygon_index(i,j))

  end function get_data

  !!===================================================================================
  !> 
  logical function get_read_mask(i, j)
    implicit none
    integer, intent(IN) :: i, j

    get_read_mask = (mask(i,j) > 0)

  end function get_read_mask

  !!===================================================================================
  !>
  logical function get_mask(i, j)
    use jcf_mesh_base, only : get_mesh_mask => get_mask
    implicit none
    integer, intent(IN) :: i, j

    get_mask = get_mesh_mask(coco, get_data_index(i,j))

  end function get_mask


  !!===================================================================================
  !>
  real(kind=8) function get_read_mask_value(i, j)
    implicit none
    integer, intent(IN) :: i, j

    get_read_mask_value = mask(i,j)

  end function get_read_mask_value


  !!===================================================================================
  !> Return x-coord of (i,j)th data point
  real(kind=8) function get_data_grid_x(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_grid_x = xt(i,j)

  end function get_data_grid_x

  !!===================================================================================
  !> Return y-coord of (i,j)th data point
  real(kind=8) function get_data_grid_y(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_data_grid_y = yt(i,j)

  end function get_data_grid_y

  !!===================================================================================
  !> Return rang(oblique) of (i,j)th data point
  real(kind=8) function get_rang(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_rang = rang(i,j)

  end function get_rang

  !!===================================================================================
  !> Return x-coord of (i,j)th volume point
  real(kind=8) function get_volume_grid_x(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_grid_x = vx(i,j)

  end function get_volume_grid_x

  !!===================================================================================
  !> Return y-coord of (i,j)th volume point
  real(kind=8) function get_volume_grid_y(i,j)
    implicit none
    integer, intent(IN) :: i, j

    get_volume_grid_y = vy(i,j)

  end function get_volume_grid_y


  !!===================================================================================
  !> return idx'th polygon(array of coord of each vertecies).
  subroutine get_polygon_from_idx( idx, poly, ierr )
    use jcf_mesh_base, only : get_point_x, get_point_y
    integer,intent(in) :: idx
    real(kind=8),allocatable,intent(out) :: poly(:,:)
    integer,intent(out) :: ierr

    integer :: i,j,l
    integer,parameter :: k = 1

    integer :: n,np
    real(kind=8) :: lon, lat
    character(len=*),parameter :: form='(F0.10,",",F0.10)'


    ierr = 0

    call get_polygon_tuple( idx, i, j )
    if ( i > NX .or. j > NY ) then
      write(0,*)'Error:coco_grid:invalid idx:',idx
      call exit(1)
    end if


    np = 4
    allocate( poly(2,np) )

    do n=1,np
      poly(1,n)=get_point_x(coco,idx,n)
      poly(2,n)=get_point_y(coco,idx,n)
    end do

    return

  end subroutine get_polygon_from_idx

  !!===================================================================================
  !> Return serial index of (i,j)th polygon.
  integer function get_polygon_index(i,j)
    implicit none
    integer, intent(IN) :: i, j

    if ( i>NX )then
      write(0,*) "Err:get_polygon_index/coco_grid: i index error:",i
      call exit(1)
    end if
    if (j>NY) then
      write(0,*) "Err:get_polygon_index/coco_grid: j index error:",j
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

  !===================================================================================


  !===================================================================================
  !> In order to get first guess for polygon search, get (i,j) from given(lon,lat).
  !!
  !! **Caution**  This routine does NOT consider tri-polar region.
  !! 
  subroutine cal_first_guess_tuple(ii,jj,lon,lat)
    implicit none
    integer, intent(OUT) :: ii
    integer, intent(OUT) :: jj
    real(kind=8), intent(IN) :: lon
    real(kind=8), intent(IN) :: lat


    ii = mod(int((lon+300)/360.d0*NX), NX)+1
    if ( ii < 1  ) ii = 1
    if ( ii > NX ) ii = NX

    jj = int((lat+90)/180.d0*NY)+1
    if ( jj < 1  ) jj = 1
    if ( jj > NY ) jj = NY

    return
  end subroutine cal_first_guess_tuple




!!$
!!$ Below are not maintained, but left as an example of use of jcf.
!!$

  !!===================================================================================
  !>
  subroutine write_data_grid(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j

    open(unit = FILE_UNIT, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do j = 1, NY
      do i = 1, NX
        if (mask(i, j) == 1.d0) then
          write(FILE_UNIT, *) get_data_point_x(coco, get_polygon_index(i,j)), get_data_point_y(coco, get_polygon_index(i,j))
        end if
      end do
    end do

    close(FILE_UNIT)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_data_grid

  !!===================================================================================
  !>
  subroutine write_data_grid_index(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j
    integer :: font_size, font_angle, font_no, grid_index
    real(kind=8) :: lon, lat
    character(len=256) :: data_str

    font_size = 4
    font_angle = 0
    font_no = 1

    open(unit = FILE_UNIT, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do j = 1, NY
      do i = 1, NX
        lon = get_data_point_x(coco, get_polygon_index(i,j))
        lat = get_data_point_y(coco, get_polygon_index(i,j))
        grid_index = 100*i+j
        write(data_str, "(2f10.5,3I2,A,I4)") lon, lat+0.1, font_size, font_angle, font_no, " BC", i
        write(FILE_UNIT, *) trim(data_str)
        write(data_str, "(2f10.5,3I2,A,I4)") lon, lat, font_size, font_angle, font_no, " TC ", j
        write(FILE_UNIT, *) trim(data_str)
      end do
    end do

    close(FILE_UNIT)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_data_grid_index

  !!===================================================================================
  !>
  subroutine write_volume_grid(file_name)
    use jcf_mesh_base, only : get_point_x, get_point_y
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j

    open(unit = FILE_UNIT, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do j = 1, NY
      do i = 1, NX
        write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(i,j),1), get_point_y(coco, get_polygon_index(i,j),1)
        write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(i,j),2), get_point_y(coco, get_polygon_index(i,j),2)
        write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(i,j),3), get_point_y(coco, get_polygon_index(i,j),3)
        write(FILE_UNIT, *) ">"
      end do
    end do

    do i = 1, NX
      write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(i,NY),1), get_point_y(coco, get_polygon_index(i,NY),1)
    end do
    write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(NX,NY),4), get_point_y(coco, get_polygon_index(NX,NY),4)
    write(FILE_UNIT, *) ">"

    do j = 1, NY
      write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(NX,j),3), get_point_y(coco, get_polygon_index(NX,j),3)
    end do
    write(FILE_UNIT, *) get_point_x(coco, get_polygon_index(NX,NY),4), get_point_y(coco, get_polygon_index(NX,NY),4)


    close(FILE_UNIT)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_volume_grid

end module coco_grid

