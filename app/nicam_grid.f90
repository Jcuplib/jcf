!------------------------------------------------------------------------------------
!> @author
!! Arakawa T.
!
!> @brief
!! read coco grid and output GMT file
!! 
!! NICAM mesh definition:
!!
!!         np6          np5
!!                 vp6
!!                    
!!           vp1        vp5
!!                     
!!      np1        dp        np4
!!                     
!!          vp2        vp4
!!             
!!                 vp3
!!          np2         np3
!!
!!   
!------------------------------------------------------------------------------------

module nicam_grid
  use jcf_mesh_base, only : mesh_type
  implicit none
  private

  public :: nicam
  public :: NORTH_POLE
  public :: SOUTH_POLE

  public :: init_grid ! subroutine (file_name)
  public :: get_grid_size ! subroutine (nnx, nny, nnl)
  public :: read_grid
  public :: get_data_lon ! real(kind=8) function (i, j, k)
  public :: get_data_lat ! real(kind=8) function (i, j, k)
  public :: get_data ! real(kind=8) function (i, j, k)
  public :: get_mask ! logical function (i, j)
  public :: get_num_of_volume_point ! integer function (i, j, k)
  public :: get_volume_grid_x ! real(kind=8) function (i, j, k, l, p)
  public :: get_volume_grid_y ! real(kind=8) function (i, j, k, l, p)
  public :: get_landmask ! integer funcion(i, j, l)
!!$  public :: write_gmt_data
  public :: get_polygon_index
  public :: index2ijl

  integer, parameter :: ADM_TI = 1
  integer, parameter :: ADM_TJ = 2
  integer, parameter :: ADM_KNONE = 1

  integer, parameter :: GRD_XDIR=1
  integer, parameter :: GRD_YDIR=2
  integer, parameter :: GRD_ZDIR=3

  integer, parameter :: FILE_UNIT = 666
  integer :: NX
  integer :: NY
  integer, parameter :: NUM_OF_DOMAIN = 10 ! number of outer domain
  integer, parameter,   private       :: GDUMMY=1
  integer :: ADM_gall
  integer :: ADM_gall_1d
  integer :: ADM_lall
  integer :: ADM_lin  ! number of inner domain, ADM_lall/NUM_OF_DOMAIN
  integer :: ADM_lin_1d ! sqrt(ADM_lin)
  integer :: NORTH_POLE
  integer :: SOUTH_POLE
  integer, parameter :: NUM_OF_VOLUME_POINT = 6

  real(kind=8), allocatable :: GRD_x(:,:,:,:) 
  real(kind=8), allocatable :: GRD_xt(:,:,:,:,:)
  real(kind=8), allocatable :: landmask(:,:,:)

  integer, allocatable :: region_info(:,:,:) ! region_info(Target_rgn:Target_direc, My_Direc, My_Rgn_ID)

  type grid_type
    real(kind=8) :: lon
    real(kind=8) :: lat
    real(kind=8) :: landmask
  end type grid_type

  type(grid_type), pointer :: data_grid(:,:,:)
  type(grid_type), pointer :: volume_grid(:,:,:,:)


  type(mesh_type) :: nicam

  integer, parameter :: STR_LEN = 128

  character(len=STR_LEN) :: hgrid_fname
  character(len=STR_LEN) :: rgnmngfname
  character(len=STR_LEN) :: landmask_fname

contains

  !===================================================================================
  !> @breaf initialize nicam grid
  !!
  !!     ADM_gall = (NX+2)*(NY+2)  
  !!
  !> @param[in] nnx number of i direction grid points 
  !> @param[in] nny number of j direction grid points
  !> @param[in] nnl number of ADM_lall 
  subroutine init_grid(cnf_file)
    use jcf_sphere_lib, only : init_sphere_lib
    use jcf_mesh_base, only : init_mesh, init_polygon
    use mod_logfile, only : LOG_FID
    implicit none

    character(len=*), intent(IN) :: cnf_file
    integer :: glevel = 5
    integer :: rlevel = 0

    namelist / nicamgrd /&
         & glevel, rlevel,&
         & rgnmngfname, hgrid_fname, landmask_fname

    integer :: ierr
    integer :: dout, din, i, j, l

    integer,parameter :: ctl_fid=111

    open(CTL_FID,             &
         file=cnf_file,   &
         form='formatted',    &
         status='old',        &
         iostat=ierr)
    if(ierr/=0) then
      write(*,*) 'Cannot open PARAMETER file!'
      stop
    end if
    rewind(ctl_fid)
    read(ctl_fid, nml = nicamgrd,end=999 )
999 close(ctl_fid)
    write(*,nml=nicamgrd)

    NX  = 2**(glevel-rlevel)
    NY  = NX
    ADM_lin   = (2**rlevel)**2
    ADM_lall  = (2**rlevel)**2 * NUM_OF_DOMAIN
    ADM_gall_1d = GDUMMY+NX+GDUMMY
    ADM_gall = ADM_gall_1d*ADM_gall_1d
    ADM_lin_1d = 2**rlevel

    write(LOG_FID, *) "Msg : Sub[init_grid]/Mod[nicam_grid]"
    write(LOG_FID, *) "----- nicam grid initialize "
    write(LOG_FID, '(A,I5)') " ----- glelve : ", glevel
    write(LOG_FID, '(A,I5)') " ----- rlelve : ", rlevel
    write(LOG_FID, '(A,I5)') " ----- NX, NY : ", NX
    write(LOG_FID, '(A,I5)') " ----- ADM_lall : ", ADM_lall
    write(LOG_FID, '(A,I5)') " ----- ADM_gall : ", ADM_gall
    write(LOG_FID, '(A,I5)') " ----- ADM_lin_1d : ", ADM_lin_1d
    write(LOG_FID, *)

    allocate(GRD_x(ADM_gall, 1, ADM_lall, 3))
    allocate(GRD_xt(ADM_gall, 1, ADM_lall, 2, 3))
    allocate(landmask(ADM_gall, 1, ADM_lall))

    allocate(region_info(2,4,ADM_lall))

    allocate(data_grid(0:NX+1, 0:NY+1, ADM_lall))
    allocate(volume_grid(0:NX+1, 0:NY+1, ADM_lall, 2))

    call init_sphere_lib()
    call init_mesh(nicam, NX*NY*ADM_lall+2, NX*NY*ADM_lall+2, 2*NX*NY*ADM_lall)

    NORTH_POLE = NX*NY*ADM_lall+1
    SOUTH_POLE = NX*NY*ADM_lall+2

    do dout = 1, NUM_OF_DOMAIN
      do din = 1, ADM_lin
        l = ADM_lin*(dout-1) + din
        do j = 1, NY
          do i = 1, NX
            if ((din==1).and.(i==1).and.(j==1)) then
              call init_polygon(nicam, get_polygon_index(i,j,l), 5) ! hexagonal initialization
            else
              call init_polygon(nicam, get_polygon_index(i,j,l), 6) ! hexagonal initialization
            end if
          end do
        end do
      end do
    end do

    call init_polygon(nicam, NORTH_POLE, 5)
    call init_polygon(nicam, SOUTH_POLE, 5)

  end subroutine init_grid

  !===================================================================================

  subroutine get_grid_size(nnx, nny, nnl)
    implicit none
    integer, intent(OUT) :: nnx, nny, nnl

    nnx = NX
    nny = NY
    nnl = ADM_lall

  end subroutine get_grid_size

  !===================================================================================

  subroutine read_grid() 
    implicit none

    call input_hgrid(hgrid_fname)
    call input_landmask(landmask_fname)
    call read_rgn_info(rgnmngfname) 

    call cal_data_grid()
    call set_nicam_grid()

  end subroutine read_grid

  !===================================================================================

  subroutine input_hgrid(basename)
    implicit none

    character(LEN=*), intent(in) :: basename
    integer :: l, d, k, n
    character(128) :: fname
    integer :: ierr

    character(len=5) :: digit

    integer,parameter :: fid = 55

    do l=1,ADM_lall

      write(digit,'(I5.5)') l-1
      fname=trim(basename)//'.rgn'//digit

      open(fid,file=trim(fname),status='old',form='unformatted')

      read(fid) ADM_gall_1d

      read(fid) GRD_x(:,ADM_KNONE,l,GRD_XDIR)
      read(fid) GRD_x(:,ADM_KNONE,l,GRD_YDIR)
      read(fid) GRD_x(:,ADM_KNONE,l,GRD_ZDIR)

      read(fid) GRD_xt(:,ADM_KNONE,l,ADM_TI:ADM_TJ,GRD_XDIR)
      read(fid) GRD_xt(:,ADM_KNONE,l,ADM_TI:ADM_TJ,GRD_YDIR)
      read(fid) GRD_xt(:,ADM_KNONE,l,ADM_TI:ADM_TJ,GRD_ZDIR)

      close(fid)

    end do

  end subroutine input_hgrid

  !===================================================================================

  subroutine input_landmask(basename)
    implicit none

    character(LEN=*), intent(in) :: basename
    integer :: l, d, k, n
    character(128) :: fname
    integer :: ierr
    integer :: ADM_gall_1d

    character(len=5) :: digit

    integer,parameter :: fid = 66

    do l=1,ADM_lall

      write(digit,'(I5.5)') l-1
      fname=trim(basename)//'.rgn'//digit

      open(fid,file=trim(fname),status='old',form='unformatted',access = "direct", recl=ADM_gall*8)

      read(fid, rec=1) landmask(:, ADM_KNONE, l)

      close(fid)

    end do


  end subroutine input_landmask

  !===================================================================================

  subroutine read_rgn_info(file_name)
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: l
    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn
    integer :: rgnid, sw(2), nw(2), ne(2), se(2)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se
    integer :: ierr

    integer,parameter :: rgn_fid = 398

    open(rgn_fid,             &
         file=file_name,   &
         form='formatted',    &
         status='old')

    rewind(rgn_fid)
    read(rgn_fid,nml=rgn_info)
    do l = 1, num_of_rgn
      read(rgn_fid, nml=rgn_link_info)
      region_info(1, 1, rgnid) = sw(1)
      region_info(2, 1, rgnid) = sw(2)
      region_info(1, 2, rgnid) = nw(1)
      region_info(2, 2, rgnid) = nw(2)
      region_info(1, 3, rgnid) = ne(1)
      region_info(2, 3, rgnid) = ne(2)
      region_info(1, 4, rgnid) = se(1)
      region_info(2, 4, rgnid) = se(2)
    end do

    close(rgn_fid)

    return

  end subroutine read_rgn_info

  !===================================================================================

  real(kind=8) function get_data_lon(i,j,k)
    implicit none
    integer, intent(IN) :: i, j, k

    get_data_lon = data_grid(i,j,k)%lon

  end function get_data_lon

  !===================================================================================

  real(kind=8) function get_data_lat(i,j,k)
    implicit none
    integer, intent(IN) :: i, j, k

    get_data_lat = data_grid(i,j,k)%lat

  end function get_data_lat

  !===================================================================================

  subroutine cal_data_grid()
    use jcf_sphere_lib, only : xyz2latlon
    implicit none
    real(kind=8) :: x, y, z
    real(kind=8) :: lat, lon
    integer :: i, j, k, ij

    do k = 1, ADM_lall
      do j = 0, NY+1
        do i = 0, NX+1
          ij = (NX+2)*j+i+1

          x = GRD_x(ij,1,k,1)
          y = GRD_x(ij,1,k,2)
          z = GRD_x(ij,1,k,3)

          call xyz2latlon(x,y,z,lat,lon)

          data_grid(i,j,k)%lon = mod(lon,360.d0)
          data_grid(i,j,k)%lat = lat

          data_grid(i,j,k)%landmask = landmask(ij,1,k)

          x = GRD_xt(ij,1,k,1,1)
          y = GRD_xt(ij,1,k,1,2)
          z = GRD_xt(ij,1,k,1,3)

          call xyz2latlon(x,y,z,lat,lon)

          volume_grid(i,j,k,1)%lon = mod(lon,360.d0)
          volume_grid(i,j,k,1)%lat = lat

          x = GRD_xt(ij,1,k,2,1)
          y = GRD_xt(ij,1,k,2,2)
          z = GRD_xt(ij,1,k,2,3)

          call xyz2latlon(x,y,z,lat,lon)

          volume_grid(i,j,k,2)%lon = mod(lon,360.d0)
          volume_grid(i,j,k,2)%lat = lat

        end do
      end do
    end do

  end subroutine cal_data_grid

  !===================================================================================

  integer function get_data_index(i, j, l)
    implicit none
    integer, intent(IN) :: i, j, l

    get_data_index = i+NX*(j-1)+NX*NY*(l-1)

  end function get_data_index

  !===================================================================================

  integer function get_volume_index(i, j, l, k)
    implicit none
    integer, intent(IN) :: i, j, l, k

    get_volume_index = k+2*(i-1)+2*NX*(j-1)+2*NX*NY*(l-1)

  end function get_volume_index

  !===================================================================================

  integer function get_polygon_index(i,j,l)
    implicit none
    integer, intent(IN) :: i, j, l

    get_polygon_index = i+NX*(j-1)+NX*NY*(l-1)

  end function get_polygon_index

  !===================================================================================

  subroutine index2ijl(idx, i, j, l)
    implicit none
    integer, intent(IN) :: idx
    integer, intent(OUT) :: i, j, l

    i = mod(idx-1, NX)+1
    l = int((idx-1)/(NX*NX))+1
    j = int((idx-NX*NX*(l-1)-1)/NX)+1

  end subroutine index2ijl

  !===================================================================================

  real(kind=8) function get_data(i, j, k)
    use jcf_mesh_base, only : get_mesh_data => get_data
    implicit none
    integer, intent(IN) :: i, j, k

    get_data = get_mesh_data(nicam, get_polygon_index(i, j, k))

  end function get_data

  !===================================================================================

  logical function get_mask(i, j, k)
    use jcf_mesh_base, only : get_mesh_mask => get_mask
    implicit none
    integer, intent(IN) :: i, j, k

    get_mask = get_mesh_mask(nicam, get_polygon_index(i,j,k))

  end function get_mask

  !===================================================================================

  integer function get_num_of_volume_point(i, j, k)
    use jcf_mesh_base, only : get_num_of_point
    integer, intent(IN) :: i, j, k

    get_num_of_volume_point = get_num_of_point(nicam, get_polygon_index(i, j, k))

  end function get_num_of_volume_point

  !===================================================================================

  real(kind=8) function get_data_grid_x(i,j,l)
    use jcf_mesh_base, only : get_data_point_x
    implicit none
    integer, intent(IN) :: i, j, l

    get_data_grid_x = get_data_point_x(nicam, get_polygon_index(i,j,l))

  end function get_data_grid_x

  !===================================================================================

  real(kind=8) function get_data_grid_y(i,j,l)
    use jcf_mesh_base, only : get_data_point_y
    implicit none
    integer, intent(IN) :: i, j, l

    get_data_grid_y = get_data_point_y(nicam, get_polygon_index(i,j,l))

  end function get_data_grid_y

  !===================================================================================

  real(kind=8) function get_volume_grid_x(i, j, l, p)
    use jcf_mesh_base, only : get_point_x
    implicit none
    integer, intent(IN) :: i, j, l, p

    get_volume_grid_x = get_point_x(nicam, get_polygon_index(i, j, l), p)

  end function get_volume_grid_x

  !===================================================================================

  real(kind=8) function get_volume_grid_y(i, j, l, p)
    use jcf_mesh_base, only : get_point_y
    implicit none
    integer, intent(IN) :: i, j, l, p

    get_volume_grid_y = get_point_y(nicam, get_polygon_index(i, j, l), p)

  end function get_volume_grid_y

  !===================================================================================

  integer function get_landmask(i, j, l)
    implicit none
    integer, intent(IN) :: i, j, l

    get_landmask = int(data_grid(i,j,l)%landmask)

  end function get_landmask

  !===================================================================================

  subroutine set_nicam_grid() 
    use jcf_mesh_base, only : set_data_point_location, set_data_point, &
         set_volume_point_location, set_volume_point, set_next_polygon
    implicit none
    integer :: sw_rgn, nw_rgn, ne_rgn, se_rgn ! region info
    integer :: sw_drc, nw_drc, ne_drc, se_drc ! direction info
    integer :: next_rgn
    integer :: polygon_index
    integer :: i, j, l, dout, din

    do l = 1, ADM_lall
      do j = 1, NY
        do i = 1, NX
          call set_data_point_location(nicam, get_data_index(i,j,l), data_grid(i,j,l)%lon, data_grid(i,j,l)%lat)
          call set_volume_point_location(nicam, get_volume_index(i,j,l,1), volume_grid(i,j,l,1)%lon, volume_grid(i,j,l,1)%lat)
          call set_volume_point_location(nicam, get_volume_index(i,j,l,2), volume_grid(i,j,l,2)%lon, volume_grid(i,j,l,2)%lat)
        end do
      end do
    end do

    call set_data_point_location(nicam, NORTH_POLE, 0.d0, 90.d0)
    call set_data_point_location(nicam, SOUTH_POLE, 0.d0, -90.d0)

    do l =1, ADM_lall
      do j = 1, NY
        do i = 1, NX
          call set_data_point(nicam, get_polygon_index(i,j,l), get_data_index(i,j,l), .true.)
        end do
      end do
    end do

    call set_data_point(nicam, NORTH_POLE, NORTH_POLE, .true.)
    call set_data_point(nicam, SOUTH_POLE, SOUTH_POLE, .true.)

    do dout = 1, NUM_OF_DOMAIN
      do din = 1, ADM_lin

        l = ADM_lin*(dout-1)+din

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
          write(0,*) "set_nicam_grid swdrc error"
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
          write(0,*) "set_nicam_grid nw_drc error"
          stop
        end select

        i = 1 ; j = 1
        if (din == 1) then
          if (l <= ADM_lall/2) then ! north
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

        if ((din>=2).and.(din<=ADM_lin_1d)) then ! sw side of outer domain
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
            write(0,*) "sw side of outer domain, direc error"
            stop
          end select
          cycle
        end if

        if (mod(din-1, ADM_lin_1d) == 0) then ! nw side of outer domain
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
            write(0,*) "nw side of outer domain, direc error"
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

    call set_volume_point(nicam, NORTH_POLE, get_volume_index(1,NY,ADM_lin*1-ADM_lin_1d+1,2))
    call set_volume_point(nicam, NORTH_POLE, get_volume_index(1,NY,ADM_lin*2-ADM_lin_1d+1,2))
    call set_volume_point(nicam, NORTH_POLE, get_volume_index(1,NY,ADM_lin*3-ADM_lin_1d+1,2))
    call set_volume_point(nicam, NORTH_POLE, get_volume_index(1,NY,ADM_lin*4-ADM_lin_1d+1,2))
    call set_volume_point(nicam, NORTH_POLE, get_volume_index(1,NY,ADM_lin*5-ADM_lin_1d+1,2))

    call set_volume_point(nicam, SOUTH_POLE, get_volume_index(NX,1,ADM_lin*5+ADM_lin_1d,1))
    call set_volume_point(nicam, SOUTH_POLE, get_volume_index(NX,1,ADM_lin*6+ADM_lin_1d,1))
    call set_volume_point(nicam, SOUTH_POLE, get_volume_index(NX,1,ADM_lin*7+ADM_lin_1d,1))
    call set_volume_point(nicam, SOUTH_POLE, get_volume_index(NX,1,ADM_lin*8+ADM_lin_1d,1))
    call set_volume_point(nicam, SOUTH_POLE, get_volume_index(NX,1,ADM_lin*9+ADM_lin_1d,1))

    ! set next polygon
    call set_next()

  end subroutine set_nicam_grid

  !===================================================================================

  subroutine set_next()
    use jcf_mesh_base, only : set_next_polygon
    use jcf_mesh_base, only : get_polygon_ptr, polygon_type

    implicit none
    integer :: sw_rgn, nw_rgn, ne_rgn, se_rgn ! region info
    integer :: sw_drc, nw_drc, ne_drc, se_drc ! direction info
    integer :: polygon_index
    integer :: i,j,l, dout, din

    type(polygon_type), pointer :: pol

    do dout = 1, NUM_OF_DOMAIN
      do din = 1, ADM_lin

        l = ADM_lin*(dout-1)+din

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
    call set_next_polygon(nicam, NORTH_POLE, get_polygon_index(1,NY,ADM_lin*2-ADM_lin_1d+1))
    call set_next_polygon(nicam, NORTH_POLE, get_polygon_index(1,NY,ADM_lin*3-ADM_lin_1d+1))
    call set_next_polygon(nicam, NORTH_POLE, get_polygon_index(1,NY,ADM_lin*4-ADM_lin_1d+1))
    call set_next_polygon(nicam, NORTH_POLE, get_polygon_index(1,NY,ADM_lin*5-ADM_lin_1d+1))
    call set_next_polygon(nicam, NORTH_POLE, get_polygon_index(1,NY,ADM_lin*1-ADM_lin_1d+1)) ! 2013/05/17 mod

    call set_next_polygon(nicam, SOUTH_POLE, get_polygon_index(NX,1,ADM_lin*5+ADM_lin_1d))
    call set_next_polygon(nicam, SOUTH_POLE, get_polygon_index(NX,1,ADM_lin*6+ADM_lin_1d))
    call set_next_polygon(nicam, SOUTH_POLE, get_polygon_index(NX,1,ADM_lin*7+ADM_lin_1d))
    call set_next_polygon(nicam, SOUTH_POLE, get_polygon_index(NX,1,ADM_lin*8+ADM_lin_1d))
    call set_next_polygon(nicam, SOUTH_POLE, get_polygon_index(NX,1,ADM_lin*9+ADM_lin_1d))

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
      write(0,*) "set_sw_side_polygon direc error"
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
      write(0,*) "set_nw_side_polygon direc error"
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
      write(0,*) "set_ne_side_polygon direc error"
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
      write(0,*) "set_se_side_polygon direc error"
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

    l = ADM_lin*(dout-1)+din
    if (din == 1) then
      if (dout <= NUM_OF_DOMAIN/2) then ! north
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

    if ((din>=2).and.(din<=ADM_lin_1d)) then ! sw side of outer domain
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
        write(0,*) "sw side of outer domain, direc error"
        stop
      end select
      return
    end if


    if (mod(din-1, ADM_lin_1d) == 0) then ! nw side of outer domain
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
        write(0,*) "nw side of outer domain, direc error"
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

    l = ADM_lin*(dout-1)+din

    if (din == ADM_lin) then     
      if (dout <= NUM_OF_DOMAIN/2) then ! north
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

    if ((din >= ADM_lin-ADM_lin_1d+1).and.(din <= ADM_lin)) then ! ne side of outer domain
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
        write(0,*) "nw side of outer domain, direc error"
        stop
      end select
      return
    end if

    if (mod(din, ADM_lin_1d) == 0) then ! se side of outer domain
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
        write(0,*) "nw side of outer domain, direc error"
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

    l = ADM_lin*(dout-1)+din

    if (din == ADM_lin_1d) then ! south corner
      if ( l <= NUM_OF_DOMAIN/2) then ! north
        goto 1000
      else ! south
        polygon_index = get_polygon_index(NX,1,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,sw_rgn))
        call set_next_polygon(nicam, polygon_index, SOUTH_POLE)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,1   ,se_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX  ,2   ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(NX-1,1   ,l))
      end if
      return
    end if


    if (din < ADM_lin_1d) then ! sw side
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
        write(0,*) "nw side of outer domain, direc error"
        stop
      end select
      return
    end if

    if (mod(din, ADM_lin_1d) == 0) then ! se side
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
        write(0,*) "nw side of outer domain, direc error"
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

    l = ADM_lin*(dout-1)+din

    if (din == ADM_lin-ADM_lin_1d+1) then ! north corner
      if (dout <= NUM_OF_DOMAIN/2) then ! north
        polygon_index = get_polygon_index(1,NY,l)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,nw_rgn))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY-1,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(2   ,NY  ,l))
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,ne_rgn))
        call set_next_polygon(nicam, polygon_index, NORTH_POLE)
        call set_next_polygon(nicam, polygon_index, get_polygon_index(1   ,NY  ,nw_rgn))
      else ! south
        goto 1000
      end if
      return
    end if

    if (mod(din-1, ADM_lin_1d) == 0) then ! nw side
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
        write(0,*) "set_nicam_grid l error"
        stop
      end select
      return
    end if

    if (din> ADM_lin-ADM_lin_1d+1) then ! ne side
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
        write(0,*) "set_nicam_grid l error"
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

  subroutine write_gmt_data(file_name)
    implicit none
    character(len=*), intent(IN) :: file_name

    call write_data_grid(trim(file_name)//"data_grid.dat")
    call write_data_grid_index(trim(file_name)//"data_grid_index.dat")
    call write_volume_grid(trim(file_name)//"volume_grid.dat")

  end subroutine write_gmt_data

  !===================================================================================

  subroutine write_data_grid(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j, k

    open(unit = FILE_UNIT, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do k = 1, ADM_lall
      do j = 1, NY
        do i = 1, NX
          write(FILE_UNIT, *) get_data_point_x(nicam, get_polygon_index(i,j,k)), &
               get_data_point_y(nicam, get_polygon_index(i,j,k))
        end do
      end do
    end do

    close(FILE_UNIT)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_data_grid

  !===================================================================================

  subroutine write_data_grid_index(file_name)
    use jcf_mesh_base, only : get_data_point_x, get_data_point_y
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j, k
    integer :: font_size, font_angle, font_no, grid_index
    real(kind=8) :: lon, lat
    character(len=256) :: data_str

    font_size = 6
    font_angle = 0
    font_no = 1

    open(unit = FILE_UNIT, FILE = trim(file_name), access = "SEQUENTIAL", ERR = 1000)

    do k = 1, ADM_lall
      do j = 1, NY
        do i = 1, NX
          lon = get_data_point_x(nicam, get_polygon_index(i,j,k))
          lat = get_data_point_y(nicam, get_polygon_index(i,j,k))
          grid_index = 100*i+j
          write(data_str, "(2f10.5,3I2,A,I3)") lon, lat+0.25, font_size, font_angle, font_no, " BC", k
          write(FILE_UNIT, *) trim(data_str)
          write(data_str, "(2f10.5,3I2,A,I5)") lon, lat, font_size, font_angle, font_no, " TC ", grid_index
          write(FILE_UNIT, *) trim(data_str)
        end do
      end do
    end do

    close(FILE_UNIT)

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_data_grid_index

  !===================================================================================

  subroutine write_volume_grid(file_name)
    use jcf_mesh_base, only : get_point_x, get_point_y
    implicit none
    character(len=*), intent(IN) :: file_name
    integer :: i, j, k, p
    character(len=128) :: grid_file_name


    do k = 1, ADM_lall
      write(grid_file_name, "(A,I4.4)") trim(file_name),k-1
      open(unit = FILE_UNIT, FILE = trim(grid_file_name), access = "SEQUENTIAL", ERR = 1000)
      do j = 1, NY
        do i = 1, NX
          if ((i==1).and.(j==1)) then
            do p = 1, 5
              write(FILE_UNIT, *) get_point_x(nicam, get_polygon_index(i,j,k), p), get_point_y(nicam, get_polygon_index(i,j,k), p)
            end do
            write(FILE_UNIT, *) get_point_x(nicam, get_polygon_index(i,j,k), 1), get_point_y(nicam, get_polygon_index(i,j,k), 1)
          else
            do p = 1, 6
              write(FILE_UNIT, *) get_point_x(nicam, get_polygon_index(i,j,k), p), get_point_y(nicam, get_polygon_index(i,j,k), p)
            end do
            write(FILE_UNIT, *) get_point_x(nicam, get_polygon_index(i,j,k), 1), get_point_y(nicam, get_polygon_index(i,j,k), 1)
          end if
          write(FILE_UNIT, *) ">"
        end do
      end do
      if (k == 1) then ! add north pole graph
        do p = 1, 5
          write(FILE_UNIT, *) get_point_x(nicam, NORTH_POLE, p), get_point_y(nicam, NORTH_POLE, p)
        end do
        write(FILE_UNIT, *) get_point_x(nicam, NORTH_POLE, 1), get_point_y(nicam, NORTH_POLE, 1)
        write(FILE_UNIT, *) ">"
      end if

      if (k == 6) then ! add south pole graph
        do p = 1, 5
          write(FILE_UNIT, *) get_point_x(nicam, SOUTH_POLE, p), get_point_y(nicam, SOUTH_POLE, p)
        end do
        write(FILE_UNIT, *) get_point_x(nicam, SOUTH_POLE, 1), get_point_y(nicam, SOUTH_POLE, 1)
        write(FILE_UNIT, *) ">"
      end if

      close(FILE_UNIT)
    end do

    return

1000 continue
    write(0,*) "file "//trim(file_name)//" open error"
    stop

  end subroutine write_volume_grid

  !===================================================================================

!!$subroutine graph_grid()
!!$  use gmt_lib, only : init_graph, set_graph_file, set_contour_color_scale, graph_global_map, &
!!$                      graph_polygon_line, graph_polygon_index, finalize_graph, graph_coastline
!!$  implicit none
!!$  real(kind=8) :: x(6), y(6), data
!!$  real(kind=4) :: min_lon, max_lon, min_lat, max_lat
!!$  integer :: i, j, k, p
!!$  
!!$
!!$  call init_graph("nicam_grid")
!!$
!!$  call set_graph_file("nicam_grid.eps")
!!$
!!$  call set_contour_color_scale(0.d0, 10.d0, 1.d0)
!!$
!!$  !min_lon = 0.0 ; max_lon = 360.0
!!$  !min_lat = -88 ; max_lat = 88 
!!$
!!$  !call graph_global_map(0.04, min_lon, max_lon, min_lat, max_lat)
!!$
!!$  min_lon = 270.0 ; max_lon = 360.0
!!$  min_lat = -88.0 ; max_lat = 0.0
!!$  call graph_global_map(0.18, min_lon, max_lon, min_lat, max_lat)
!!$
!!$  do k = 1, ADM_lall
!!$    do j = 1, NY
!!$      do i = 1, NX
!!$        do p = 1, get_num_of_volume_point(i, j, k)
!!$          x(p) = get_volume_grid_x(i,j,k,p) ; y(p) = get_volume_grid_y(i,j,k,p)
!!$        end do
!!$        data = 5.d0
!!$        if ((min_lon <= maxval(x)).and.(minval(x) <= max_lon).and.(min_lat <= maxval(y)).and.(minval(y) <= max_lat)) then
!!$          call graph_polygon_line(get_num_of_volume_point(i,j,k), x, y)
!!$        end if
!!$  
!!$      end do
!!$    end do
!!$  end do
!!$
!!$  do k = 1, ADM_lall
!!$    do j = 1, NY
!!$      do i = 1, NX
!!$        do p = 1, get_num_of_volume_point(i, j, k)
!!$          x(p) = get_volume_grid_x(i,j,k,p) ; y(p) = get_volume_grid_y(i,j,k,p)
!!$        end do
!!$        if ((min_lon <= maxval(x)).and.(minval(x) <= max_lon).and.(min_lat <= maxval(y)).and.(minval(y) <= max_lat)) then
!!$          x(1) = get_data_grid_x(i,j,k) ; y(1) = get_data_grid_y(i,j,k)
!!$          call grapH_polygon_index(x(1), y(1), 1, k, i*100+j)
!!$        end if
!!$      end do
!!$    end do
!!$  end do
!!$
!!$  call graph_coastline()
!!$  
!!$  call finalize_graph()
!!$
!!$end subroutine graph_grid
!!$
  !===================================================================================

!!$subroutine graph_index()
!!$  use gmt_lib, only : init_graph, set_graph_file, set_contour_color_scale, graph_global_map, &
!!$                      graph_polygon_data, graph_polygon_index, finalize_graph, graph_coastline
!!$  implicit none
!!$  real(kind=8) :: x(6), y(6), data
!!$  real(kind=4) :: min_lon, max_lon, min_lat, max_lat
!!$  integer :: i, j, k, p
!!$  
!!$
!!$  call init_graph("nicam_grid")
!!$
!!$  call set_graph_file("nicam_landmask.eps")
!!$
!!$  call set_contour_color_scale(0.d0, 10.d0, 1.d0)
!!$
!!$  min_lon = 0.0 ; max_lon = 360.0
!!$  min_lat = -88 ; max_lat = 88 
!!$
!!$  call graph_global_map(0.04, min_lon, max_lon, min_lat, max_lat)
!!$
!!$  !min_lon = 270.0 ; max_lon = 360.0
!!$  !min_lat = -88.0 ; max_lat = 0.0
!!$  !call graph_global_map(0.18, min_lon, max_lon, min_lat, max_lat)
!!$
!!$  do k = 1, ADM_lall
!!$    do j = 1, NY
!!$      do i = 1, NX
!!$        do p = 1, get_num_of_volume_point(i, j, k)
!!$          x(p) = get_volume_grid_x(i,j,k,p) ; y(p) = get_volume_grid_y(i,j,k,p)
!!$        end do
!!$        data = data_grid(i,j,k)%landmask
!!$        if ((min_lon <= maxval(x)).and.(minval(x) <= max_lon).and.(min_lat <= maxval(y)).and.(minval(y) <= max_lat)) then
!!$          call graph_polygon_data(get_num_of_volume_point(i,j,k), x, y, data)
!!$        end if
!!$  
!!$      end do
!!$    end do
!!$  end do
!!$
!!$  !do k = 1, ADM_lall
!!$  !  do j = 1, NY
!!$  !    do i = 1, NX
!!$  !      do p = 1, get_num_of_volume_point(i, j, k)
!!$  !        x(p) = get_volume_grid_x(i,j,k,p) ; y(p) = get_volume_grid_y(i,j,k,p)
!!$  !      end do
!!$  !      if ((min_lon <= maxval(x)).and.(minval(x) <= max_lon).and.(min_lat <= maxval(y)).and.(minval(y) <= max_lat)) then
!!$  !        x(1) = get_data_grid_x(i,j,k) ; y(1) = get_data_grid_y(i,j,k)
!!$  !        call grapH_polygon_index(x(1), y(1), 1, k, i*100+j)
!!$  !      end if
!!$  !    end do
!!$  !  end do
!!$  !end do
!!$
!!$  call graph_coastline()
!!$  
!!$  call finalize_graph()
!!$
!!$end subroutine graph_index


end module nicam_grid
