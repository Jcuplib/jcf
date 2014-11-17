!> Read and maintain NICAM grid.
!!
!! This module manages/maintains type(mesh_type) instance for NICAM
!! ''Icosahedoral'' grid.
!!
!! NICAM uses 'A-grid' and all physical quantities are defined on the
!! vertecies of triangles formed by icosahedral grid. So a shape of
!! the control volume is a hexagon, some of singular points a
!! pentagon.
!!
!!
!! NICAM mesh definition in this module is as follows;
!!      
!!           np6          np5
!!                  vp6
!!      
!!            vp1        vp5
!!      
!!       np1        dp        np4
!!      
!!            vp2        vp4
!!      
!!                  vp3
!!            np2         np3
!!      
!!      
!!
!! - dp : data point
!! - vpX: volume point number X 
!! - npX: next(neighbour) polygon X's data point
!! .
!!
!!
!! \note
!! For a unique way of domain decomposition of NICAM, each grid point
!! and control volume(polygon) are indexed by a triplet (i,j,l),
!! rather than a doublet(i,j) like COCO or usual lat-lon coordinate.
!!
!! \author
!! - Arakawa T.

module nicam_grid
  use jcf_mesh_base, only : mesh_type
  implicit none

  public :: nicam

  public :: init_grid
  public :: read_grid

  public :: get_data
!!$  public :: get_mask
  public :: get_data_lon
  public :: get_data_lat
  public :: get_data_grid_x
  public :: get_data_grid_y
  public :: get_volume_grid_x
  public :: get_volume_grid_y

  public :: get_grid_size
  public :: get_num_of_volume_point

  public :: get_data_index
  public :: get_volume_index
  public :: get_polygon_index
  public :: get_polygon_tuple

  public :: get_polygon_from_idx

  integer, parameter :: FNLEN = 1024

  type(mesh_type) :: nicam !< type(mesh_type) instance for NICAM grid.

  integer,private :: NX
  integer,private :: NY
  integer,private :: NL

  integer, parameter,private :: DMD = 10 !< number of outer-most domain of ICO grid.
  integer, parameter,private :: GDUMMY = 1 !< width of halo region
  integer,private :: gall
  integer,private :: gall1d
  integer,private :: gall_pl
  integer,private :: lall
  integer,private :: lall_pl
  integer,private :: lin        ! number of inner domain, lall/DMD
  integer,private :: lin_1d     ! sqrt(lin)

  integer,public :: TOP_POLE
  integer,public :: BOT_POLE

  real(kind=8), allocatable,private :: grid_x     (:,:,:,:)
  real(kind=8), allocatable,private :: grid_xt    (:,:,:,:,:)
  real(kind=8), allocatable,private :: grid_x_pl  (:,:,:,:)
  real(kind=8), allocatable,private :: grid_xt_pl (:,:,:,:)

  integer, allocatable,private :: region_info(:,:,:) ! region_info(Target_rgn:Target_direc, My_Direc, My_Rgn_ID)

  type grid_type
    real(kind=8) :: lon
    real(kind=8) :: lat
!!$    real(kind=8) :: landmask
  end type grid_type

  type(grid_type), pointer,private :: data_grid  (:,:,:)   ! grid center
  type(grid_type), pointer,private :: volume_grid(:,:,:,:) ! grid vertex(up,down)


  integer :: glevel = 5 !< glevel
  integer :: rlevel = 0 !< rlevel
  integer :: vlayer = 1                     !< dummy here
  character(len=FNLEN) :: rgnmngfname = ""  !< region info filename
  character(len=FNLEN) :: hgrid_fname = ""  !< horizontal grid filename
  character(len=FNLEN) :: vgrid_fname = ""  !< dummy here
  character(len=FNLEN) :: topo_fname  = ""  !< dummy here

  namelist / nicamgrd / &
       glevel,                 &
       rlevel,                 &
       vlayer,                 &
       rgnmngfname,            &
       hgrid_fname,            &
       vgrid_fname,            &
       topo_fname


contains

  !===================================================================================
  !> Initialize type(mesh_type) instance for nicam grid.
  !!
  !! Ether `cnf_file` or `cnf_fid` must be supplied.
  !! - `cnf_file` is supplied, open/close this file in this routine,
  !! - `cnf_fid` is supplied, assumes that cnf file is already opened,
  !! - None of `cnf_file`/`cnf_fid` is given, cause error.
  !! .
  subroutine init_grid(cnf_file, cnf_fid)
    use jcf_spherical_lib, only: &
         & init_spherical_lib
    use jcf_mesh_base, only: &
         & init_mesh, &
         & init_polygon
    use jcf_misc, only: &
         & jcf_avail_fid, &
         & LOG_FID => jcf_log_fid
    implicit none
    character(len=*), intent(IN), optional :: cnf_file !< namelist file name.
    integer         , intent(IN), optional :: cnf_fid  !< namelist file lun.

    integer :: fid
    character(len=FNLEN) :: fname

    logical :: open_cnf_myself
    integer :: i, j, l
    integer :: ierr
    integer :: dout, din


    ! Open cnf(namelist) file if cnf_file is given, else just rewind.
    if ( present(cnf_fid) ) then
      fid = cnf_fid
      open_cnf_myself = .false.
    else if ( present(cnf_file) ) then
      call open_cnf_file()
      open_cnf_myself = .true.
    else
      write(0,'(A)')'ERR:init_grid/nicam_grid: Must specify ether fid or fname.'
      call exit(1)
    end if

    rewind(fid)
    read(fid,nml=nicamgrd,iostat=ierr)
    if ( ierr < 0 ) then
      write(LOG_FID,'(A)') '*** init_grid/nicam_grid: nicamgrd is not specified, use default.'
    elseif( ierr > 0 ) then
      write(0,      '(A)') 'ERR:init_grid/nicam_grid: nicamgrd read error.'
      call exit(1)
    end if
!!$    write(LOG_FID,nml=nicamgrd)
    if ( open_cnf_myself ) then
      close(fid)
    end if

    NX      = 2**(glevel-rlevel)
    NY      = NX
    NL      = (2**rlevel)**2 * DMD
    lall    = NL
    lin     = (2**rlevel)**2
    lin_1d  = 2**rlevel

    lall_pl = 2
    gall1d  = NX+2
    gall    = gall1d * gall1d
    gall_pl = 6

    TOP_POLE = NX*NY*lall+1
    BOT_POLE = NX*NY*lall+2

    write(LOG_FID, '(A)')   "*** init_grid/nicam_grid:"
    write(LOG_FID,'(A15,": ",I0)') "glevel", glevel
    write(LOG_FID,'(A15,": ",I0)') "rlevel", rlevel
    write(LOG_FID,'(A15,": ",I0)') "NX", NX
    write(LOG_FID,'(A15,": ",I0)') "NY", NY
    write(LOG_FID,'(A15,": ",I0)') "NL", NL
    write(LOG_FID,'(A15,": ",I0)') "gall", gall
    write(LOG_FID,'(A15,": ",I0)') "lin_1d", lin_1d
    write(LOG_FID,'(A15,": ",I0)') "TOP_POLE", TOP_POLE
    write(LOG_FID,'(A15,": ",I0)') "BOT_POLE", BOT_POLE
    write(LOG_FID,*)

    allocate( region_info(2,4,lall) )

    allocate( grid_x     (gall,1,lall,  3) )
    allocate( grid_xt    (gall,1,lall,2,3) )
    grid_x    (:,:,:,:)   = 0.D0
    grid_xt   (:,:,:,:,:) = 0.D0
    allocate( grid_x_pl  (gall_pl,1,lall_pl,3) )
    allocate( grid_xt_pl (gall_pl,1,lall_pl,3) )
    grid_x_pl (:,:,:,:)   = 0.D0
    grid_xt_pl(:,:,:,:) = 0.D0


    allocate( data_grid  (0:NX+1,0:NY+1,lall+lall_pl  ) )
    allocate( volume_grid(0:NX+1,0:NY+1,lall+lall_pl,2) )

    call init_spherical_lib()
    call init_mesh(nicam, NX*NY*lall+2, NX*NY*lall+2, 2*NX*NY*lall)

    do dout = 1, DMD
      do din  = 1, lin
        l = lin*(dout-1) + din
        do j = 1, NY
          do i = 1, NX
            if ( din == 1 .and. i == 1 .and. j == 1 ) then
              call init_polygon(nicam, get_polygon_index(i,j,l), 5) ! pentagonal initialization
            else
              call init_polygon(nicam, get_polygon_index(i,j,l), 6) ! hexagonal initialization
            end if
          end do
        end do
      end do
    end do

    call init_polygon(nicam, TOP_POLE, 5)
    call init_polygon(nicam, BOT_POLE, 5)

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
        write(0,'(A)') 'ERR:init_grid/coco_grid: Cannot open cnf file:'//trim(fname)
        call exit(1)
      end if
    end subroutine open_cnf_file



  end subroutine init_grid


  !!===================================================================================
  !> Read nicam_grid data file and setup mesh_type instance.
  !!
  !! Read in rgninfo file and hgrid files.
  subroutine read_grid()
    use jcf_misc, only: &
         & LOG_FID => jcf_log_fid
    implicit none

    call read_rgn_info  ( trim(rgnmngfname) )
    write(LOG_FID,'(A)') '*** read_grid/nicam_grid:rgn_info file was read successfully: '//trim(rgnmngfname)
    call input_hgrid    ( trim(hgrid_fname) )
    write(LOG_FID,'(A)') '*** read_grid/nicam_grid:hgrid files were read successfully: '//trim(hgrid_fname)
    call cal_data_grid
    call set_grid

    return
  end subroutine read_grid


  !!===================================================================================
  !> Return value of (i,j,l)'th polygon (data_point)
  real(kind=8) function get_data(i, j, l)
    use jcf_mesh_base, only : get_mesh_data => get_data
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_data = get_mesh_data(nicam, get_polygon_index(i,j,l))

  end function get_data


  !!===================================================================================
  !> Return mask of (i,j,l)'th polygon
  !!
  !! __Tentatively disabled__.
  !!
  logical function get_mask()
!!$    use jcf_mesh_base, only : get_mesh_mask => get_mask
!!$    implicit none
!!$    integer, intent(IN) :: i, j, k
!!$
!!$    get_mask = get_mesh_mask(nicam, get_polygon_index(i,j,k))
!!$
  end function get_mask


  !!===================================================================================
  !> Return longitude(deg) of (i,j,l)'th data_point
  real(kind=8) function get_data_lon(i,j,l)
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_data_lon = data_grid(i,j,l)%lon

  end function get_data_lon


  !!===================================================================================
  !> Return latitude(deg) of (i,j,l)'th data_point
  real(kind=8) function get_data_lat(i,j,l)
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_data_lat = data_grid(i,j,l)%lat

  end function get_data_lat


  !!===================================================================================
  !> Return x-coord of (i,j,l)'th polygon's data_point
  !!
  !! \todo Is there any difference with get_data_lon() above ??
  real(kind=8) function get_data_grid_x(i,j,l)
    use jcf_mesh_base, only : get_data_point_x
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_data_grid_x = get_data_point_x(nicam, get_polygon_index(i,j,l))

  end function get_data_grid_x


  !!===================================================================================
  !> Returen y-coord of (i,j,l)'th polygon's data_point
  real(kind=8) function get_data_grid_y(i,j,l)
    use jcf_mesh_base, only : get_data_point_y
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_data_grid_y = get_data_point_y(nicam, get_polygon_index(i,j,l))

  end function get_data_grid_y


  !!===================================================================================
  !> Return x-coord of p'th volume point of (i,j,l)'th polygon.
  real(kind=8) function get_volume_grid_x(i, j, l, p)
    use jcf_mesh_base, only : get_point_x
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple
    integer, intent(IN) :: p !< up/down of volume grid.

    get_volume_grid_x = get_point_x(nicam, get_polygon_index(i, j, l), p)

  end function get_volume_grid_x


  !!===================================================================================
  !> Return y-coord of p'th volume point of (i,j,l)'th polygon.
  real(kind=8) function get_volume_grid_y(i, j, l, p)
    use jcf_mesh_base, only : get_point_y
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple
    integer, intent(IN) :: p !< up/down of volume grid.

    get_volume_grid_y = get_point_y(nicam, get_polygon_index(i, j, l), p)

  end function get_volume_grid_y


  !!===================================================================================
  !> Return total grid size
  !!
  !!
  subroutine get_grid_size(nnx, nny, nnl)
    implicit none
    integer, intent(OUT) :: nnx !< X size
    integer, intent(OUT) :: nny !< Y size
    integer, intent(OUT) :: nnl !< L size

    nnx = NX
    nny = NY
    nnl = NL

    return
  end subroutine get_grid_size


  !!===================================================================================
  !> Return number of volume points of (i,j,l)'th polygon
  integer function get_num_of_volume_point(i, j, l)
    use jcf_mesh_base, only : get_num_of_point
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_num_of_volume_point = get_num_of_point(nicam, get_polygon_index(i, j, l))

  end function get_num_of_volume_point


  !!===================================================================================
  !> Return serialized index of (i,j,l)'th data_point
  integer function get_data_index(i, j, l)
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_data_index = i+NX*(j-1)+NX*NY*(l-1)

  end function get_data_index


  !!===================================================================================
  !> Return serialized index of p'th volume point of (i,j,l)'th polygon.
  integer function get_volume_index(i, j, l, p)
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple
    integer, intent(IN) :: p !< up/down of volume grid.

    get_volume_index = p+2*(i-1)+2*NX*(j-1)+2*NX*NY*(l-1)

  end function get_volume_index


  !!===================================================================================
  !> Return serialized index of (i,j,l)'th polygon
  integer function get_polygon_index(i,j,l)
    implicit none
    integer, intent(IN) :: i !< data point tuple
    integer, intent(IN) :: j !< data point tuple
    integer, intent(IN) :: l !< data point tuple

    get_polygon_index = i+NX*(j-1)+NX*NY*(l-1)

  end function get_polygon_index


  !!===================================================================================
  !> Return tuple (i,j,l) of idx'th polygon
  !!
  !! As you know, this is the reverse function of get_polygon_index() above.
  subroutine get_polygon_tuple(idx, i, j, l)
    implicit none
    integer, intent(IN)  :: idx !< polygon index
    integer, intent(OUT) :: i   !< polygon tuple
    integer, intent(OUT) :: j   !< polygon tuple
    integer, intent(OUT) :: l   !< polygon tuple

    i = mod(idx-1, NX)+1
    l = int((idx-1)/(NX*NX))+1
    j = int((idx-NX*NX*(l-1)-1)/NX)+1

  end subroutine get_polygon_tuple


  !!===================================================================================
  !> Return polygon coords of idx'th polygon.
  !!
  !! Polygon coords `poly(2,NP)` is an array of (x,y) of each polygon vertices.
  !!
  !! \note poly is deallocated here even if pre-allocated before call.
  subroutine get_polygon_from_idx( idx, poly, ierr )
    use jcf_mesh_base, only : get_point_x, get_point_y
    integer,intent(in) :: idx   !< polygon's index
    real(kind=8),allocatable,intent(out) :: poly(:,:) !< coords of polygon
    integer,intent(out) :: ierr !< non-zefo if error.

    integer :: i,j,l
    integer,parameter :: k = 1

    integer :: n,np
    real(kind=8) :: lon, lat

    ierr = 0

    call get_polygon_tuple( idx, i, j, l )
    np = get_num_of_volume_point(i,j,l)

!!$d  write(0,*)'dbg:get_polygon_from_idx:'
!!$d  write(0,*)'dbg:idx,i,j,l:',idx,i,j,l
!!$d  write(0,*)'dbg:num_of_volume_point:',np

    if ( allocated ( poly ) ) then
      deallocate( poly )
    end if
    allocate( poly(2,np) )


    do n=1,np
      poly(1,n)=get_point_x(nicam,idx,n)
      poly(2,n)=get_point_y(nicam,idx,n)
    end do

    return

  end subroutine get_polygon_from_idx

!!$
!!$
!!$ Below are private routines.
!!$
!!$
  !!===================================================================================
  !! read rgn_info file.
  subroutine read_rgn_info(fname)
    use jcf_misc, only: &
         & jcf_avail_fid
    implicit none

    character(len=*), intent(in) :: fname

    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn

    integer :: rgnid, sw(2), nw(2), ne(2), se(2)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se

    integer :: fid, ierr
    integer :: l

    fid = jcf_avail_fid()
    open( unit   = fid,         &
         file   = trim(fname), &
         form   = 'formatted', &
         status = 'old')

    rewind(fid)
    read(fid,nml=rgn_info)
    do l = 1, num_of_rgn
      read(fid,nml=rgn_link_info)

      region_info(1,1,rgnid) = sw(1)
      region_info(2,1,rgnid) = sw(2)
      region_info(1,2,rgnid) = nw(1)
      region_info(2,2,rgnid) = nw(2)
      region_info(1,3,rgnid) = ne(1)
      region_info(2,3,rgnid) = ne(2)
      region_info(1,4,rgnid) = se(1)
      region_info(2,4,rgnid) = se(2)
    end do

    close(fid)

    return
  end subroutine read_rgn_info



  !!===================================================================================
  !! readin hgrid files
  subroutine input_hgrid(basename)
    use jcf_misc, only: &
         & jcf_avail_fid, &
         & LOG_FID => jcf_log_fid
    implicit none

    character(len=*), intent(in) :: basename

    character(len=128) :: fname

    integer :: fid, ierr
    integer :: l

    integer :: i,j,suf
    suf(i,j) = gall1d * ((j)-1) + (i)

    logical,parameter :: bgrid_dump = .true.

    do l = 1, lall
      call consist_rgn_fname(fname,basename,l)
      fid = jcf_avail_fid()
      open( unit   = fid,           &
           file   = trim(fname),   &
           form   = 'unformatted', &
           status = 'old')

      read(fid) gall1d

      read(fid) grid_x(:,1:1,l,1)
      read(fid) grid_x(:,1:1,l,2)
      read(fid) grid_x(:,1:1,l,3)

      read(fid) grid_xt(:,1:1,l,1:2,1)
      read(fid) grid_xt(:,1:1,l,1:2,2)
      read(fid) grid_xt(:,1:1,l,1:2,3)

      close(fid)
    end do

    ! fill unused grid
    grid_x (suf(gall1d,1),:,:,:)   = grid_x (suf(gall1d,2),:,:,:)
    grid_x (suf(1,gall1d),:,:,:)   = grid_x (suf(2,gall1d),:,:,:)

    grid_xt(suf(gall1d,1),:,:,:,:) = grid_xt(suf(gall1d,2),:,:,:,:)
    grid_xt(suf(1,gall1d),:,:,:,:) = grid_xt(suf(2,gall1d),:,:,:,:)


    !! TOP_POLE, BOTTOM_POLE
    fname=trim(basename)//'.pl'
    fid =jcf_avail_fid()
    open(  unit   = fid,           &
         & file   = fname,         &
         & status = 'old',         &
         & form   = 'unformatted', &
         & iostat = ierr           )

    if ( ierr /= 0 ) then
      write(0,'(A)') 'ERR:input_hgrid/nicam_grid: Cannot open hgrid(.pl) file!', trim(fname)
      call exit(1)
    end if

    read(fid) grid_x_pl(:,:,:,:)
    if(bgrid_dump) then
      read(fid) grid_xt_pl(:,:,:,:)
    end if

    close(fid)

    return
  end subroutine input_hgrid


  !===================================================================================

  subroutine cal_data_grid()
!!$  use jcf_sphere_lib, only: xyz2latlon
    use jcf_spherical_lib, only: xyz2lonlat
    implicit none

    real(8) :: x, y, z
    real(8) :: lat, lon
    integer :: i, j, l, ij, ll

    do l = 1, lall
      do j = 0, NY+1
        do i = 0, NX+1
          ij = (NX+2)*j+i+1

          x = grid_x(ij,1,l,1)
          y = grid_x(ij,1,l,2)
          z = grid_x(ij,1,l,3)

          call xyz2lonlat(x,y,z,lon,lat)

          data_grid(i,j,l)%lon = mod(lon,360.D0)
          data_grid(i,j,l)%lat = lat

!!$     data_grid(i,j,l)%landmask = landmask(ij,1,l)

          x = grid_xt(ij,1,l,1,1)
          y = grid_xt(ij,1,l,1,2)
          z = grid_xt(ij,1,l,1,3)

          call xyz2lonlat(x,y,z,lon,lat)

          volume_grid(i,j,l,1)%lon = mod(lon,360.D0)
          volume_grid(i,j,l,1)%lat = lat

          x = grid_xt(ij,1,l,2,1)
          y = grid_xt(ij,1,l,2,2)
          z = grid_xt(ij,1,l,2,3)

          call xyz2lonlat(x,y,z,lon,lat)

          volume_grid(i,j,l,2)%lon = mod(lon,360.D0)
          volume_grid(i,j,l,2)%lat = lat

        end do
      end do
    end do

    !! For Top/Bottom pole
    !! \todo where these ''magic number'' come from ??
    !! \todo 0...5 should be 1...6 ??
    do ll = 1,2
      do i = 0,5
        l = lall + ll
        ij = i+1
        x = grid_x_pl(ij,1,ll,1)
        y = grid_x_pl(ij,1,ll,2)
        z = grid_x_pl(ij,1,ll,3)
        call xyz2lonlat(x,y,z,lon,lat)

        data_grid(i,1,l)%lon = mod(lon,360.D0)
        data_grid(i,1,l)%lat = lat
        !! \todo How about landmask ??
!!$    data_grid(ij,1,l)%landmask = landmask(ij,1,l)

        x = grid_xt_pl(ij,1,ll,1)
        y = grid_xt_pl(ij,1,ll,2)
        z = grid_xt_pl(ij,1,ll,3)
        call xyz2lonlat(x,y,z,lon,lat)

        volume_grid(i,1,l,1)%lon = mod(lon,360.D0)
        volume_grid(i,1,l,1)%lat = lat
        !! for TOP/BOTTOM pole, no ''down'' volume grid.

      end do
    end do

!!$    write(0,*)'dbg:cal_data_grid:data_grid(PL):'
!!$    do l=lall+1,lall+2
!!$      do i=0,5
!!$        write(0,*)data_grid(i,1,l)%lon,data_grid(i,1,l)%lat
!!$      end do
!!$    end do

    return
  end subroutine cal_data_grid


  !===================================================================================

  subroutine set_grid()
    use jcf_mesh_base, only : &
         & set_data_point_location, &
         & set_data_point, &
         & set_volume_point_location, &
         & set_volume_point, &
         & set_next_polygon
    implicit none

    integer :: sw_rgn, nw_rgn, ne_rgn, se_rgn ! region info
    integer :: sw_drc, nw_drc, ne_drc, se_drc ! direction info
    integer :: next_rgn
    integer :: polygon_index

    integer :: i, j, l, dout, din

    do l = 1, lall
      do j = 1, NY
        do i = 1, NX
          call set_data_point_location  ( nicam, get_data_index  (i,j,l),   data_grid(i,j,l)%lon,     data_grid(i,j,l)%lat     )
          call set_volume_point_location( nicam, get_volume_index(i,j,l,1), volume_grid(i,j,l,1)%lon, volume_grid(i,j,l,1)%lat )
          call set_volume_point_location( nicam, get_volume_index(i,j,l,2), volume_grid(i,j,l,2)%lon, volume_grid(i,j,l,2)%lat )
        end do
      end do
    end do

!!$  call set_data_point_location( nicam, TOP_POLE, 0.D0,  90.D0 )
!!$  call set_data_point_location( nicam, BOT_POLE, 0.D0, -90.D0 )
    call set_data_point_location( nicam, TOP_POLE, data_grid(0,1,lall+1)%lon, data_grid(0,1,lall+1)%lat )
    call set_data_point_location( nicam, BOT_POLE, data_grid(0,1,lall+2)%lon, data_grid(0,1,lall+2)%lat )

    do l =1, lall
      do j = 1, NY
        do i = 1, NX
          call set_data_point(nicam, get_polygon_index(i,j,l), get_data_index(i,j,l), .true.)
        end do
      end do
    end do

    call set_data_point(nicam, TOP_POLE, TOP_POLE, .true.)
    call set_data_point(nicam, BOT_POLE, BOT_POLE, .true.)

    do dout = 1, DMD
      do din = 1, lin

        l = lin*(dout-1)+din

        do j = 2, NY
          do i = 2, NX
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,j-1,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,j-1,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j-1,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,j  ,l,1))
          end do
        end do

        sw_rgn = region_info(1, 1, l)
        sw_drc = region_info(2, 1, l)
        nw_rgn = region_info(1, 2, l)
        nw_drc = region_info(2, 2, l)
        ne_rgn = region_info(1, 3, l)
        ne_drc = region_info(2, 3, l)
        se_rgn = region_info(1, 4, l)
        se_drc = region_info(2, 4, l)

        ! SW side
        j = 1
        select case (sw_drc)
        case(3) !
          do i = 2, NX
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,NY ,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,NY ,sw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,NY ,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,j  ,l,1))
          end do
        case(4) ! south pole
          do i = 2, NX
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(NX ,NY-i+2,sw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX ,NY-i+1,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX ,NY-i+1,sw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j    ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j    ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i-1,j    ,l,1))
          end do
        case default
          write(0,'(A)') "set_grid swdrc error"
          stop
        end select

        ! NW side
        i = 1
        select case(nw_drc)
        case(4)
          do j = 2, NY
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,j-1,nw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,j-1,nw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j-1,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,j,nw_rgn,1))
          end do
        case(3) ! north pole
          do j = 2, NY
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(NX-j+1,NY,nw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX-j+2,NY,nw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j-1,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX-j+1,NY,nw_rgn,2))
          end do
        case default
          write(0,'(A)') "set_grid nw_drc error"
          stop
        end select

        i = 1 ; j = 1
        if (din == 1) then
          if (l <= lall/2) then ! north
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,nw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,nw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index( 1,NY,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index( i, j, l, 1))
            call set_volume_point(nicam, polygon_index, get_volume_index( i, j, l, 2))
          else ! south
            polygon_index = get_polygon_index(i,j,l)
            call set_volume_point(nicam, polygon_index, get_volume_index(NX, 1,nw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,sw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index( i, j, l, 1))
            call set_volume_point(nicam, polygon_index, get_volume_index( i, j, l, 2))
          end if
          cycle
        end if

        if ((din>=2).and.(din<=lin_1d)) then ! sw side of outer domain
          select case (sw_drc)
          case(3) ! normal
            polygon_index = get_polygon_index(i,j,l)
            next_rgn = region_info(1, 1, nw_rgn) ! sw region of nw region
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,next_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,next_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i ,NY,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX, 1,nw_rgn,1))
          case(4)
            polygon_index = get_polygon_index(i,j,l)
            next_rgn = region_info(1, 3, sw_rgn) ! ne region of sw region
            call set_volume_point(nicam, polygon_index, get_volume_index(NX, 1,next_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,sw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX, 1, nw_rgn,1))
          case default
            write(0,'(A)') "sw side of outer domain, direc error"
            stop
          end select
          cycle
        end if

        if (mod(din-1, lin_1d) == 0) then ! nw side of outer domain
          select case (nw_drc)
          case(4) ! normal
            polygon_index = get_polygon_index(i,j,l)
            next_rgn = region_info(1, 1, nw_rgn) ! sw region of nw region
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,next_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,next_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i ,NY,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX, 1,nw_rgn,1))
          case(3)
            polygon_index = get_polygon_index(i,j,l)
            next_rgn = region_info(1, 4, nw_rgn) ! se region of nw region
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,nw_rgn,1))
            call set_volume_point(nicam, polygon_index, get_volume_index( 1,NY,next_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i ,NY,sw_rgn,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
            call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
            call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,nw_rgn,2))
          case default
            write(0,'(A)') "nw side of outer domain, direc error"
            stop
          end select
          cycle
        end if

        polygon_index = get_polygon_index(i,j,l)
        next_rgn = region_info(1, 1, nw_rgn) ! sw region of nw region
        call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,next_rgn,2))
        call set_volume_point(nicam, polygon_index, get_volume_index(NX,NY,next_rgn,1))
        call set_volume_point(nicam, polygon_index, get_volume_index(i ,NY,sw_rgn,2))
        call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,1))
        call set_volume_point(nicam, polygon_index, get_volume_index(i  ,j  ,l,2))
        call set_volume_point(nicam, polygon_index, get_volume_index(NX, 1,nw_rgn,1))

      end do
    end do

    call set_volume_point(nicam, TOP_POLE, get_volume_index(1,NY,lin*1-lin_1d+1,2))
    call set_volume_point(nicam, TOP_POLE, get_volume_index(1,NY,lin*2-lin_1d+1,2))
    call set_volume_point(nicam, TOP_POLE, get_volume_index(1,NY,lin*3-lin_1d+1,2))
    call set_volume_point(nicam, TOP_POLE, get_volume_index(1,NY,lin*4-lin_1d+1,2))
    call set_volume_point(nicam, TOP_POLE, get_volume_index(1,NY,lin*5-lin_1d+1,2))

    call set_volume_point(nicam, BOT_POLE, get_volume_index(NX,1,lin*5+lin_1d,1))
    call set_volume_point(nicam, BOT_POLE, get_volume_index(NX,1,lin*6+lin_1d,1))
    call set_volume_point(nicam, BOT_POLE, get_volume_index(NX,1,lin*7+lin_1d,1))
    call set_volume_point(nicam, BOT_POLE, get_volume_index(NX,1,lin*8+lin_1d,1))
    call set_volume_point(nicam, BOT_POLE, get_volume_index(NX,1,lin*9+lin_1d,1))

    ! set next polygon
    call set_next()


    ! Find out which polygon contains real North/South pole ??
    call set_polar_polygon(nicam)


  end subroutine set_grid








  !===================================================================================

  subroutine set_next()
    use jcf_mesh_base, only : set_next_polygon
    use jcf_mesh_base, only : get_polygon_ptr, polygon_type

    implicit none

    integer :: sw_rgn, nw_rgn, ne_rgn, se_rgn ! region info
    integer :: sw_drc, nw_drc, ne_drc, se_drc ! direction info
    integer :: polygon_index

    integer :: i,j,l, dout, din

    do dout = 1, DMD
      do din = 1, lin

        l = lin*(dout-1)+din

        ! set inner area
        do j = 2, NY-1
          do i = 2, NX-1
            polygon_index = get_polygon_index(i,j,l)
            call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,j-1,l))
            call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,j-1,l))
            call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,j  ,l))
            call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,j+1,l))
            call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,j+1,l))
            call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,j  ,l))
          end do
        end do

        sw_rgn = region_info(1, 1, l)
        sw_drc = region_info(2, 1, l)
        nw_rgn = region_info(1, 2, l)
        nw_drc = region_info(2, 2, l)
        ne_rgn = region_info(1, 3, l)
        ne_drc = region_info(2, 3, l)
        se_rgn = region_info(1, 4, l)
        se_drc = region_info(2, 4, l)

        ! set SW side
        call set_sw_side_polygon(l, sw_rgn, sw_drc)

        ! set NW side
        call set_nw_side_polygon(l, nw_rgn, nw_drc)

        ! set NE side
        call set_ne_side_polygon(l, ne_rgn, ne_drc)

        ! set SE side
        call set_se_side_polygon(l, se_rgn, se_drc)

        ! set (1,1) polygon
        call set_1_1_polygon(dout, din, nw_rgn, nw_drc, sw_rgn, sw_drc)

        ! set (NX,NY) polygon
        call set_NX_NY_polygon(dout, din, ne_rgn, ne_drc, se_rgn, se_drc)

        ! set (NX,1) polygon
        call set_NX_1_polygon(dout, din, sw_rgn, sw_drc, se_rgn, se_drc)

        ! set (1,NY) polygon
        call set_1_NY_polygon(dout, din, nw_rgn, nw_drc, ne_rgn, ne_drc)

      end do
    end do

    ! north pole
    call set_next_polygon(nicam, TOP_POLE, get_polygon_index(1,NY,lin*2-lin_1d+1))
    call set_next_polygon(nicam, TOP_POLE, get_polygon_index(1,NY,lin*3-lin_1d+1))
    call set_next_polygon(nicam, TOP_POLE, get_polygon_index(1,NY,lin*4-lin_1d+1))
    call set_next_polygon(nicam, TOP_POLE, get_polygon_index(1,NY,lin*5-lin_1d+1))
    call set_next_polygon(nicam, TOP_POLE, get_polygon_index(1,NY,lin*1-lin_1d+1)) ! 2013/05/17 mod

    call set_next_polygon(nicam, BOT_POLE, get_polygon_index(NX,1,lin*5+lin_1d))
    call set_next_polygon(nicam, BOT_POLE, get_polygon_index(NX,1,lin*6+lin_1d))
    call set_next_polygon(nicam, BOT_POLE, get_polygon_index(NX,1,lin*7+lin_1d))
    call set_next_polygon(nicam, BOT_POLE, get_polygon_index(NX,1,lin*8+lin_1d))
    call set_next_polygon(nicam, BOT_POLE, get_polygon_index(NX,1,lin*9+lin_1d))

  end subroutine set_next

  !===================================================================================

  subroutine set_sw_side_polygon(my_rgn, sw_rgn, sw_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: my_rgn, sw_rgn, sw_drc
    integer :: polygon_index
    integer :: i,j,l

    l = my_rgn
    j = 1
    select case (sw_drc)
    case(3) ! north east side
      do i = 2, NX-1
        polygon_index = get_polygon_index(i,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,NY ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,NY ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,1  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,1+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,1+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,1  ,l))
      end do
    case(4) ! south east side
      do i = 2, NX-1
        polygon_index = get_polygon_index(i,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX ,NY-i+2,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX ,NY-i+1,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,1  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,1+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,1+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,1  ,l))
      end do
    case default
      write(0,'(A)') "set_sw_side_polygon direc error"
      stop
    end select

  end subroutine set_sw_side_polygon

  !===================================================================================

  subroutine set_nw_side_polygon(my_rgn, nw_rgn, nw_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: my_rgn, nw_rgn, nw_drc
    integer :: polygon_index
    integer :: i,j,l

    l = my_rgn
    i = 1
    select case(nw_drc)
    case(4) ! south east side
      do j = 2, NY-1
        polygon_index = get_polygon_index(1,j,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,j-1,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1  ,j-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,j  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,j+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1  ,j+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,j,nw_rgn))
      end do
    case(3) ! north east side
      do j = 2, NY-1
        polygon_index = get_polygon_index(1,j,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-j+2,NY,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1  ,j-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,j  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,j+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1  ,j+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-j+1,NY,nw_rgn))
      end do
    case default
      write(0,'(A)') "set_nw_side_polygon direc error"
      stop
    end select

  end subroutine set_nw_side_polygon

  !===================================================================================

  subroutine set_ne_side_polygon(my_rgn, ne_rgn, ne_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: my_rgn, ne_rgn, ne_drc
    integer :: polygon_index
    integer :: i,j,l

    l = my_rgn
    j = NY
    select case(ne_drc)
    case(1) ! south west side
      do i = 2, NX-1
        polygon_index = get_polygon_index(i,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,NY  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,1,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,1,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,NY  ,l))
      end do
    case(2) ! north west side
      do i = 2, NX-1
        polygon_index = get_polygon_index(i,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i  ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i+1,NY  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1,NY-i+1,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1,NY-i+2,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(i-1,NY  ,l))
      end do
    case default
      write(0,'(A)') "set_ne_side_polygon direc error"
      stop
    end select

  end subroutine set_ne_side_polygon

  !===================================================================================

  subroutine set_se_side_polygon(my_rgn, se_rgn, se_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: my_rgn, se_rgn, se_drc
    integer :: polygon_index
    integer :: i,j,l

    l = my_rgn
    i = NX
    select case(se_drc)
    case(2) ! north west side
      do j = 2, NY-1
        polygon_index = get_polygon_index(NX,j,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,j-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,j-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1,j  ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1,j+1,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,j+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,j  ,l))
      end do
    case(1) ! south west side
      do j = 2, NY-1
        polygon_index = get_polygon_index(NX,j,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,j-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,j-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-j+2,1,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-j+1,1,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,j+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,j  ,l))
      end do
    case default
      write(0,'(A)') "set_se_side_polygon direc error"
      stop
    end select

  end subroutine set_se_side_polygon

  !===================================================================================

  subroutine set_1_1_polygon(dout, din, nw_rgn, nw_drc, sw_rgn, sw_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: dout, din, nw_rgn, nw_drc, sw_rgn, sw_drc
    integer :: polygon_index, next_rgn
    integer :: l

    l = lin*(dout-1)+din
    if (din == 1) then
      if (dout <= DMD/2) then ! north
        polygon_index = get_polygon_index(1,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,NY,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 1,NY,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,1  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,1+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1  ,1+1,l))
      else ! south
        polygon_index = get_polygon_index(1,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX, 1,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,NY,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,1  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1+1,1+1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1  ,1+1,l))
      end if
      return
    end if

    if ((din>=2).and.(din<=lin_1d)) then ! sw side of outer domain
      select case (sw_drc)
      case(3) ! normal
        goto 1000
      case(4)
        polygon_index = get_polygon_index(1,1,l)
        next_rgn = region_info(1, 3, sw_rgn) ! ne region of sw region
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX, 1,next_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,NY,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 2, 1, l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 1, 1, l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 1, 2, l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX, 1, nw_rgn))
      case default
        write(0,'(A)') "sw side of outer domain, direc error"
        stop
      end select
      return
    end if


    if (mod(din-1, lin_1d) == 0) then ! nw side of outer domain
      select case (nw_drc)
      case(4) ! normal
        goto 1000
      case(3)
        polygon_index = get_polygon_index(1,1,l)
        next_rgn = region_info(1, 4, nw_rgn) ! se region of nw region
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 1,NY,next_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 1,NY,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 1, 2, l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 2, 2, l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index( 2, 1, l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,NY,nw_rgn))
      case default
        write(0,'(A)') "nw side of outer domain, direc error"
        stop
      end select
      return
    end if

1000 continue

    polygon_index = get_polygon_index(1,1,l)
    next_rgn = region_info(1, 1, nw_rgn) ! sw region of nw region
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX,NY,next_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index( 1,NY,sw_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index( 2, 1, l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index( 2, 2, l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index( 1, 2, l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX, 1,nw_rgn))

  end subroutine set_1_1_polygon

  !===================================================================================

  subroutine set_NX_NY_polygon(dout, din, ne_rgn, ne_drc, se_rgn, se_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: dout, din, ne_rgn, ne_drc, se_rgn, se_drc
    integer :: polygon_index
    integer :: l, next_rgn

    l = lin*(dout-1)+din

    if (din == lin) then
      if (dout <= DMD/2) then ! north
        polygon_index = get_polygon_index(NX,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,2   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,l))
      else ! south
        polygon_index = get_polygon_index(NX,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,l))
      end if
      return
    end if

    if ((din >= lin-lin_1d+1).and.(din <= lin)) then ! ne side of outer domain
      select case(ne_drc)
      case(1) ! normal
        goto 1000
      case(2)
        polygon_index = get_polygon_index(NX,NY,l)
        next_rgn = region_info(1, 1, ne_rgn) ! sw region of ne region
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,2   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,l))
      case default
        write(0,'(A)') "nw side of outer domain, direc error"
        stop
      end select
      return
    end if

    if (mod(din, lin_1d) == 0) then ! se side of outer domain
      select case(se_drc)
      case(2) ! normal
        goto 1000
      case(1)
        polygon_index = get_polygon_index(NX,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,l))
      case default
        write(0,'(A)') "nw side of outer domain, direc error"
        stop
      end select

      return
    end if


1000 continue

    polygon_index = get_polygon_index(NX,NY,l)
    next_rgn = region_info(1, 4, ne_rgn) ! se region of ne region
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY-1,l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,se_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,next_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,ne_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,l))

  end subroutine set_NX_NY_polygon

  !===================================================================================

  subroutine set_NX_1_polygon(dout, din, sw_rgn, sw_drc, se_rgn, se_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: dout, din, sw_rgn, sw_drc, se_rgn, se_drc
    integer :: polygon_index
    integer :: l

    l = lin*(dout-1)+din

    if (din == lin_1d) then ! south corner
      if ( l <= DMD/2) then ! north
        goto 1000
      else ! south
        polygon_index = get_polygon_index(NX,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, BOT_POLE)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,1   ,l))
      end if
      return
    end if


    if (din < lin_1d) then ! sw side
      select case(sw_drc)
      case(3) ! normal
        goto 1000
      case(4)
        polygon_index = get_polygon_index(NX,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,2   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,1   ,l))
      case default
        write(0,'(A)') "nw side of outer domain, direc error"
        stop
      end select
      return
    end if

    if (mod(din, lin_1d) == 0) then ! se side
      select case(se_drc)
      case(2) ! normal
        goto 1000
      case(1)
        polygon_index = get_polygon_index(NX,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY  ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,1   ,l))
      case default
        write(0,'(A)') "nw side of outer domain, direc error"
        stop
      end select
      return
    end if


1000 continue

    polygon_index = get_polygon_index(NX,1,l)
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,NY  ,sw_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY  ,sw_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,se_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,2   ,se_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,1   ,l))


  end subroutine set_NX_1_polygon

  !===================================================================================

  subroutine set_1_NY_polygon(dout, din, nw_rgn, nw_drc, ne_rgn, ne_drc)
    use jcf_mesh_base, only : set_next_polygon
    implicit none
    integer, intent(IN) :: dout, din, nw_rgn, nw_drc, ne_rgn, ne_drc
    integer :: polygon_index
    integer :: l, next_rgn

    l = lin*(dout-1)+din

    if (din == lin-lin_1d+1) then ! north corner
      if (dout <= DMD/2) then ! north
        polygon_index = get_polygon_index(1,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, TOP_POLE)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,nw_rgn))
      else ! south
        goto 1000
      end if
      return
    end if

    if (mod(din-1, lin_1d) == 0) then ! nw side
      select case(nw_drc)
      case(4) ! normal
        goto 1000
      case(3)
        polygon_index = get_polygon_index(1,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,1   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,nw_rgn))
      case default
        write(0,'(A)') "set_grid l error"
        stop
      end select
      return
    end if

    if (din> lin-lin_1d+1) then ! ne side
      select case(ne_drc)
      case(1) ! normal
        goto 1000
      case(2)
        polygon_index = get_polygon_index(1,NY,l)
        next_rgn = region_info(1, 3, ne_rgn) ! ne region of ne region
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,next_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY  ,nw_rgn))
      case default
        write(0,'(A)') "set_grid l error"
        stop
      end select
      return
    end if

1000 continue

    polygon_index = get_polygon_index(1,NY,l)
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY-1,nw_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY-1,l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,l))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,1   ,ne_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,1   ,ne_rgn))
    call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,NY  ,nw_rgn))

  end subroutine set_1_NY_polygon

  !===================================================================================

  subroutine set_polar_polygon(this)
    use jcf_mesh_base,only : polygon_type, is_in_this_polygon,get_data_point_x,get_data_point_y
    use jcf_misc, only: LOG_FID => jcf_log_fid
    type(mesh_type), intent(INOUT) :: this !! this mesh

    type(polygon_type),pointer :: polygon

    integer :: nn
    logical :: res
    integer :: found_np, found_sp
    integer :: index_np, index_sp

    found_np=0
    found_sp=0
    index_np=-1
    index_sp=-1
    do nn = this%num_of_polygon, 1, -1
!!$    write(0,*)'dbg:set_polar_polygon:nn:',nn
      polygon => this%polygon(nn)
      res = is_in_this_polygon(0.d0, 90.d0,polygon,is_lat_lon=.false.)
      if ( res ) then
        polygon%including_np=.true.
        found_np = found_np+1
        index_np = nn
      end if
      res = is_in_this_polygon(0.d0,-90.d0,polygon,is_lat_lon=.false.)
      if ( res ) then
        polygon%including_sp =.true.
        found_sp = found_sp+1
        index_sp = nn
      end if
    end do

    write(LOG_FID,'(A)')'*** set_polar_polygon/nicam_grid: Searching North/South Pole polygon done.'
    write(LOG_FID,'(5X,A,": ",I0)')'North Pole Index',index_np
    write(LOG_FID,'(5X,A,": ",I0)')'South Pole Index',index_sp

    if     ( found_np < 0 ) then
      write(LOG_FID,'(A)')'WRN:set_polar_polygon/nicam_grid: No North Pole polygon found!'
    elseif ( found_np > 1 ) then
      write(LOG_FID,'(A)')'WRN:set_polar_polygon/nicam_grid: TOO MANY North Pole polygon found!'
    else
      if ( index_np > 0 ) then
        write(LOG_FID,'(A,2F10.5)')'*** set_polar_polygon/nicam_grid:NP point:',&
             & get_data_point_x(nicam,index_np),get_data_point_y(nicam,index_np)
      end if
    end if

    if     ( found_sp < 0 ) then
      write(LOG_FID,'(A)')'WRN:set_polar_polygon/nicam_grid: No South Pole polygon found!'
    elseif ( found_sp > 1 ) then
      write(LOG_FID,'(A)')'WRN:set_polar_polygon/nicam_grid: TOO MANY South Pole polygon found!'
    else
      if ( index_sp > 0 ) then
        write(LOG_FID,'(A,2F10.5)')'*** set_polar_polygon/nicam_grid:SP point:',&
             & get_data_point_x(nicam,index_sp),get_data_point_y(nicam,index_sp)
      end if
    end if


  end subroutine set_polar_polygon


  !!========================================================================
  !! consist filename for NICAM region data file.
  !!
  !! NICAM region data file name consist of:
  !!   `basename`+".rgn"+`XXXXX`
  !! here, XXXXX is a five-digit-region-number.
  !!
  !! \attention
  !! `XXXXX` is rgn_number-1.
  subroutine consist_rgn_fname(filename,basename,rgnid)
    character(len=*),intent(out) :: filename
    character(len=*),intent(in)  :: basename
    integer,         intent(in)  :: rgnid

    character(len=5) :: cnum

    if ( len(filename) .lt. len_trim(basename)+9 ) then
      write(0,'(A)')'ERR:consist_rgn_fname/nicam_grid:insufficient length for filename.'
      call exit(1)
    end if
    
    write(cnum,'(I5.5)')rgnid -1

    filename=trim(basename)//'.rgn'//cnum

    return
  end subroutine consist_rgn_fname




!!$
!!$
!!$ Below are not maintained, but left as an example of use of jcf.
!!$
!!$

  !===================================================================================

  subroutine write_data_grid(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    use jcf_misc, only: jcf_avail_fid
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j, k

    integer :: fid

    fid = jcf_avail_fid()
    open(unit = fid, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do k = 1, lall
      do j = 1, NY
        do i = 1, NX
          write(fid, *) get_data_point_x(nicam, get_polygon_index(i,j,k)),&
               &        get_data_point_y(nicam, get_polygon_index(i,j,k))
        end do
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
    use jcf_misc, only : jcf_avail_fid
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j, k
    integer :: font_size, font_angle, font_no, grid_index
    real(kind=8) :: lon, lat
    character(len=256) :: data_str

    integer :: fid

    font_size = 6
    font_angle = 0
    font_no = 1

    fid = jcf_avail_fid()
    open(unit = fid, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do k = 1, lall
      do j = 1, NY
        do i = 1, NX
          lon = get_data_point_x(nicam, get_polygon_index(i,j,k))
          lat = get_data_point_y(nicam, get_polygon_index(i,j,k))
          grid_index = 100*i+j
          write(data_str, "(2f10.5,3I2,A,I3)") lon, lat+0.25, font_size, font_angle, font_no, " BC", k
          write(fid, *) trim(data_str)
          write(data_str, "(2f10.5,3I2,A,I5)") lon, lat, font_size, font_angle, font_no, " TC ", grid_index
          write(fid, *) trim(data_str)
        end do
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
    use jcf_misc, only : jcf_avail_fid
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j, k, p
    character(len=128) :: grid_file_name

    integer :: fid

    do k = 1, lall
      write(grid_file_name, "(A,I4.4)") trim(file_name),k-1
      open(unit = fid, FILE = trim(grid_file_name), access = "SEQUENTIAL", ERR = 1000)
      do j = 1, NY
        do i = 1, NX
          if ((i==1).and.(j==1)) then
            do p = 1, 5
              write(fid, *) get_point_x(nicam, get_polygon_index(i,j,k), p), get_point_y(nicam, get_polygon_index(i,j,k), p)
            end do
            write(fid, *) get_point_x(nicam, get_polygon_index(i,j,k), 1), get_point_y(nicam, get_polygon_index(i,j,k), 1)
          else
            do p = 1, 6
              write(fid, *) get_point_x(nicam, get_polygon_index(i,j,k), p), get_point_y(nicam, get_polygon_index(i,j,k), p)
            end do
            write(fid, *) get_point_x(nicam, get_polygon_index(i,j,k), 1), get_point_y(nicam, get_polygon_index(i,j,k), 1)
          end if
          write(fid, *) ">"
        end do
      end do
      if (k == 1) then ! add north pole graph
        do p = 1, 5
          write(fid, *) get_point_x(nicam, TOP_POLE, p), get_point_y(nicam, TOP_POLE, p)
        end do
        write(fid, *) get_point_x(nicam, TOP_POLE, 1), get_point_y(nicam, TOP_POLE, 1)
        write(fid, *) ">"
      end if

      if (k == 6) then ! add south pole graph
        do p = 1, 5
          write(fid, *) get_point_x(nicam, BOT_POLE, p), get_point_y(nicam, BOT_POLE, p)
        end do
        write(fid, *) get_point_x(nicam, BOT_POLE, 1), get_point_y(nicam, BOT_POLE, 1)
        write(fid, *) ">"
      end if

      close(fid)
    end do

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_volume_grid


end module nicam_grid
