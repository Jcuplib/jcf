!> Module for handling mapping table file.
!!
!! \par Composition of mapping table file.
!!
!! Mapping table file consists of header and correspondance records.
!!
!! Header consists of following five records;
!! 1. character(len=mtf_clen) :: crecv ! name and description of receiver compo.
!! 2. integer                 :: nrecv ! num  of polygons of receiver compo.
!! 3. character(len=mtf_clen) :: csend ! name and description of sender compo.
!! 4. integer                 :: nsend ! num of polygons of sender compo.
!! 5. integer(2)              :: nrec, ncf ! num of correspondants records,
!!                              and num of coeffs in one correspondance record.
!! .
!!
!! The first __word__ (before first whitespace) of `crecv`/`csend` is
!! considered as the name of receiver/sender name, the rest is an
!! description.  For example:
!! 
!! ~~~~~~~~~~~~~~~~~~~~~~~~~
!!  cname = crecv(1:index(crecv,' ')-1) !! name
!!  cdesc = crecv(index(crecv,' ')+1:) !! description
!!  if ( cname .ne. yourcomponame) then
!!     !! you have wrong data
!!  end if
!! ~~~~~~~~~~~~~~~~~~~~~~~~~
!! or,
!! ~~~~~~~~~~~~~~~~~~~~~~~~~
!!  if ( index(crecv, yourcomponame) > 0 ) then
!!    !! do somthing for the right component, ...
!!  end if
!! ~~~~~~~~~~~~~~~~~~~~~~~~~
!! 
!! These header information and other attributes for maptab file are
!! stored as an instance of maptabfile_t.
!!
!! Each of correspondance data record are as follows in one record (line):
!! - integer :: ridx ! index of receiver polygon.
!! - integer :: sidx ! index of sender polygon.
!! - real(kind=8),allocatable :: coef(:) ! array of coefficients
!! .
!! 
!! One correspondence recored should be stored as one element of
!! maptab_t array.
!!
!! Mapping table file may be a text (formatted) file or a binary
!! (unformatted) file.  In both format, contents described above are
!! same.
!!
!! All of public routines whose name start with `jcf_mtf_` are act as
!! like an member function of maptabfile_t, so take mapfiletab_t instance
!! `self` as the first argument.
!!
!! Similarly, public routines whose name start with `jcf_mtab_` are
!! for the maptab_t.
!!
!!
!! \par Simple usage: for writing.
!! 
!! For writing mapping table file with this module, you have to:
!! 1. Prepare maptabfile_t instance (`mtf` for example) by calling
!!    jcf_mtf_init() and jcf_mtf_set().  You have to set __ALL__ of
!!    `name`, `binary`, `crecv`, `nrecv`, `csend`, `nsend`, `nrec` and
!!    `ncf`.
!!    Then open mapping table file by calling jcf_mtf_open() with
!!    `act="Write"`.
!!
!! 2. Prepare correspondance data as an array of maptab_t (`mtab(:)`
!!    for example).  You should allocate mtab(:) by calling
!!    jcf_mtab_alloc() with the same value of `nrec` and `ncf` above.
!!    If you want to know necessary record size, you should use
!!    count_correspondance() in jcf_mesh_base.
!!    Then set each element of `mtab` by calling jcf_mtab_set() in do loop.
!!
!! 3. Write out header info first with jcf_mtf_write_head(), then
!!    write out whole correspondance data by jcf_mtf_write_whole_records().
!!
!! 4. Do not forget to call jcf_mtf_close() when everything done.   
!! .
!!
!! Of cource you may not use `mtab` in 2. above, and use
!! jcf_mtf_write_one_record() within a do loop in 3..
!!
!!
!! \par Simple usage: for reading.
!! 
!! For reading mapping table file with this module, you have to:
!! 1. Prepare maptabfile_t instance (say `mtf`) as above for writing,
!!    but you have only to set `name` in calling jcf_mtf_set().
!!    Then open mapping table file by calling jcf_mtf_open() with
!!    `act="Read"`.
!!
!! 2. Before preparing an array of maptab_t, you should read out
!!    header info by jcf_mtf_read_head(), and check if this file is
!!    the right one you want, by comparing `csend`/`crecv` with your
!!    component name.
!!
!! 3. Allocate an array of maptab_t by calling jcf_mtab_alloc() with
!!    the size obtained from header info above.
!! 
!! 4. Read out correspondance data by calling jcf_mtf_read_whole_records().
!!
!! 5. Again, do not forget to call jcf_mtf_close() when everything done.
!! .
!!
!! See also test_mtf_read_write.f90 in test/ directory.
!!
!! \warning
!! Currently, num of coefs in one correspondence record MUST NOT be
!! more than 5 for text format.
!!
!! \note
!! Routines in this module does NO message output except quite a few
!! warining or when some error happens.  In such case, output is to
!! STDERR(unit=0), and `call exit(1)` in any error.
!!
!!
!! \author
!! Takahiro INOUE <tinoue@rist.jp>
!!
module jcf_maptabfile
  use jcf_misc, only:&
       & avail_fid => jcf_avail_fid
  implicit none 

  public :: maptabfile_t
  public :: maptab_t

  public :: jcf_mtf_init
  public :: jcf_mtf_dump
  public :: jcf_mtf_set
  public :: jcf_mtf_open
  public :: jcf_mtf_close
  public :: jcf_mtf_read_head
  public :: jcf_mtf_read_whole_records
  public :: jcf_mtf_read_one_record
  public :: jcf_mtf_write_head
  public :: jcf_mtf_write_whole_records
  public :: jcf_mtf_write_one_record

  public :: jcf_mtab_alloc
  public :: jcf_mtab_dealloc
  public :: jcf_mtab_set

  !! \todo should be public ??
  integer,parameter,public :: dp_k = 8

  integer,parameter,public :: mtf_fnlen = 1024 !< file name length
  integer,parameter,public :: MTF_CLEN = 80    !< crecv/csend length

  !>  class for maptabfile.
  !!
  !! \note All members are __private__, except istat.
  !!
  !! \note Initial values for each member are defined in jcf_mtf_init().
  type maptabfile_t
    private
    integer,public :: istat     !< non-zero if something happens.
    integer :: lun              !< logical unit number for maptabfile.
    character(len=mtf_fnlen) :: name !< maptab filename
    logical :: binary                !< .false. if text format
    logical :: opened                !< .true. after opened
    character(len=mtf_clen) :: crecv !< name/desc of receiver
    integer :: nrecv                 !< num of receiver polygons
    character(len=mtf_clen) :: csend !< name/desc of sender 
    integer :: nsend                 !< num of sender polygons
    integer :: nrec                  !< num of correspondence
    integer :: ncf                   !< num of coefs in one correspondence
    integer :: irec                  !! file position in maptab file
  end type maptabfile_t


  !> Class mapping table correspondance record.
  !!
  !! Should be used as an array of this type, size is `nrec` of maptabfile_t.
  !!
  !! \note size of coef MUST be the same with `ncf` of maptabfile_t.
  type maptab_t
    integer      :: ridx                !< index of receiver polygon
    integer      :: sidx                !< index of sender polygon
    real(kind=8),allocatable :: coef(:) !< coefficients
  end type maptab_t


  integer,   parameter,public :: mtf_imiss = -999    !< integer missing value.
  real(dp_k),parameter,public :: mtf_fmiss = -999.d0 !< double precision missing value.

  !                                               12345678901234
  character(len=*),parameter,private :: rec_form='(2I10,5E24.15)'
  character(len=*),parameter,private :: head_form(5) = (/&
       &  '(1A80)','(1I16)','(1A80)','(1I16)','(2I16)' /)


contains
  !!========================================================================
  !> Initialize maptabfile_t.
  !!
  subroutine jcf_mtf_init(self)
    type(maptabfile_t),intent(out) :: self !< maptabfile_t instance

    self=maptabfile_t(  &
         & istat  = 0,  &
         & lun    = 0,  &
         & name   = 'XXXX',&
         & binary = .false.,&
         & opened = .false.,&
         & crecv  = '', &
         & nrecv  = 0,  &
         & csend  = '', &
         & nsend  = 0,  &
         & nrec   = -1, &
         & ncf    =  1, &
         & irec   = -1  )

    return
  end subroutine jcf_mtf_init


  !!========================================================================
  !> Dump maptabfile_t
  !!
  !! Information will be out to `unit` if specified, otherwise to STDOUT(=6).
  !!
  !! `head` is prepended if specified.
  !! 
  subroutine jcf_mtf_dump(self,unit,head)
    type(maptabfile_t),intent(in) :: self !< maptabfile_t instance
    integer,           intent(in),optional :: unit !< logical unit num for output.
    character(len=*),  intent(in),optional :: head !< heading for dump info.
    integer :: lun = 6

    if ( present(unit) ) lun=unit

    if ( present(head) )  then
      write(lun,'(1X,A)')   head
    else
      write(lun,'(1X,A)')   '##### dump maptabfile_t #####'
    end if
    write(lun,'(A,A)')  '       name: ',trim(self%name)
    write(lun,'(A,L5)') '     binary: ',self%binary
    write(lun,'(A,L5)') '     opened: ',self%opened
    write(lun,'(A,A)')  '      crecv: ',trim(self%crecv)
    write(lun,'(A,I10)')'      nrecv: ',self%nrecv
    write(lun,'(A,A)')  '      csend: ',trim(self%csend)
    write(lun,'(A,I10)')'      nsend: ',self%nsend
    write(lun,'(A,I10)')'       nrec: ',self%nrec
    write(lun,'(A,I10)')'        ncf: ',self%ncf
!!$    write(lun,'(A,I10)')'      irec: ',self%irec
!!$    write(lun,'(A,I10)')'     istat: ',self%istat

    return
  end subroutine jcf_mtf_dump


  !!========================================================================
  !> Set some attributes of maptabfile_t in advance.
  !!
  !! Each argument  corresponds to a member with the same name.
  !! 
  subroutine jcf_mtf_set(self, &
       & name   , &
       & binary ,&
       & crecv  ,&
       & nrecv  ,&
       & csend  ,&
       & nsend  ,&
       & nrec   ,&
       & ncf )

    type(maptabfile_t),intent(inout) :: self        !< maptabfile_t instance
    character(len=*) ,intent(in),optional :: name   !< name of maptab file
    logical          ,intent(in),optional :: binary !< .true. if binary file
    character(len=*) ,intent(in),optional :: crecv  !< name/desc of receiver
    integer          ,intent(in),optional :: nrecv  !< num of polygons of receiver 
    character(len=*) ,intent(in),optional :: csend  !< name/desc of sender
    integer          ,intent(in),optional :: nsend  !< num of polygons of sender
    integer          ,intent(in),optional :: nrec   !< num of correspondance records.
    integer          ,intent(in),optional :: ncf    !< num of coefs in one correspondance.


    if ( present(name) ) self%name=name
    if ( present(binary) ) self%binary=binary

    if ( present(crecv) ) self%crecv = crecv
    if ( present(nrecv) ) self%nrecv = nrecv
    if ( present(csend) ) self%csend = csend
    if ( present(nsend) ) self%nsend = nsend
    if ( present(nrec) ) self%nrec = nrec
    if ( present(ncf) ) self%ncf = ncf
       

    return

  end subroutine jcf_mtf_set

  !!========================================================================
  !> Open mappingtable file with setting of maptabfile_t.
  !!
  !! If `fname` specified, open it and override self\%name.
  !!
  !! `act` can be arbitary length, but only first letter matters.
  !!
  !! Currently, whether specified file is binary or not is checked internally,
  !! not depending on preset self\%binary, so you don't need to preset it.
  !!
  !! If something happens, self\%istat .ne. 0
  !!
  subroutine jcf_mtf_open(self, act, fname)
    type(maptabfile_t),intent(inout) :: self  !< maptabfile.
    character(len=*),  intent(in)    :: act   !< 'R' or 'W'
    character(len=*),  intent(in), optional :: fname !< name of maptablefile.

    integer :: lun

    logical :: lex=.false.
    logical :: lbin=.false.

    character(len=mtf_fnlen) :: name
    character(len=16)  :: cform !! 'formatted' or 'unformatted'
    character(len=16)  :: cstat !! 'old' or 'unknown' for 'R'ead or 'Write'

    character(len=80)  :: emsg =''

    integer :: ios=0

    if ( present( fname ) ) then
      self%name=fname
    end if

    if ( self%name .eq. 'XXXX') then
      write(0,*)'ERR:jcf_mtf_open:fname NOT specified.'
      call exit(1)
    end if


    select case ( act(1:1) )
    case( 'R', 'r' )
      inquire(file=self%name,exist=lex)

      if ( .not. lex ) then
        write(0,*)'ERR:jcf_mtf_open: file NOT exist: '//trim(self%name)
        self%istat=-1
        call exit(1)
      end if
      cstat='old'
      call check_if_binary_file(lbin,self%name)
    case( 'W', 'w' )
      cstat='unknown'
      lbin = self%binary
    case default
      write(0,*)'ERR:jcf_mtf_open: Invarid act:'//trim(act)
      call exit(1)
    end select


    if (lbin) then
      cform = 'unformatted'
    else
      cform = 'formatted'
    end if

    lun = avail_fid()
    open(unit=lun,file=self%name,form=cform,status=cstat,iostat=ios,iomsg=emsg)
    self%istat=ios
    if ( ios .ne. 0 ) then
      write(0,*)'ERR:jcf_mtf_open:Open error: '//trim(emsg)
      return
    end if

    rewind(lun)

    self%lun=lun
    self%binary=lbin
    self%opened=.true.
    self%istat=0
    self%irec=-1

    !! Is this preferable ??
!!$    call jcf_mtf_dump(self,unit=0,head='dbg:jcf_mtf_open:Opened:')
    return

  end subroutine jcf_mtf_open



  !!========================================================================
  !> Close mappingtable file.
  !!
  !! \warning `self` is re-initialized after this.
  subroutine jcf_mtf_close( self )
    type(maptabfile_t),intent(inout) :: self

    close(self%lun)

    call jcf_mtf_init( self )

    return
  end subroutine jcf_mtf_close


  !!========================================================================
  !> Read out header of maptabfile and set maptabfile_t members.
  !!
  !! All arguments other than `self` are optional, you can specify only you need.
  !!
  !! \note
  !! After calling this routine, file point is set the head of first record.
  !! 
  subroutine jcf_mtf_read_head(self,crecv,nrecv,csend,nsend,nrec,ncf)
    type(maptabfile_t),intent(inout) :: self  !< maptab file
    character(len=MTF_CLEN),intent(out),optional :: crecv !< name and desc of receiver 
    integer,                intent(out),optional :: nrecv !< num of receiver polygon
    character(len=MTF_CLEN),intent(out),optional :: csend !< name and desc of sender
    integer,                intent(out),optional :: nsend !< num of sender polygon
    integer,                intent(out),optional :: nrec  !< num of record
    integer,                intent(out),optional :: ncf   !< num of coefs in one correspondance.

    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_read_head: File NOT Opened Yet!'
      call exit(1)
    end if

    if ( self%binary ) then
      rewind(self%lun)
      read(self%lun) self%crecv
      read(self%lun) self%nrecv
      read(self%lun) self%csend
      read(self%lun) self%nsend 
      read(self%lun) self%nrec, self%ncf
    else
      rewind(self%lun)
      read(self%lun,head_form(1)) self%crecv
      read(self%lun,head_form(2)) self%nrecv
      read(self%lun,head_form(3)) self%csend
      read(self%lun,head_form(4)) self%nsend 
      read(self%lun,head_form(5)) self%nrec, self%ncf  
    end if

    self%istat = 0
    self%irec  = 0

    if ( present( crecv ) ) crecv = self%crecv
    if ( present( nrecv ) ) nrecv = self%nrecv
    if ( present( csend ) ) csend = self%csend
    if ( present( nsend ) ) nsend = self%nsend
    if ( present( nrec  ) ) nrec  = self%nrec
    if ( present( ncf   ) ) ncf   = self%ncf

    return
  end subroutine jcf_mtf_read_head


  !!========================================================================
  !> Read whole records and return as an array of maptab_t.
  !!
  !! \note Must allocate mtab(:) BEFORE calling me.
  !! 
  subroutine jcf_mtf_read_whole_records(self, mtab)
    type(maptabfile_t),intent(inout) :: self    !< maptabfile_t instance
    type(maptab_t),    intent(inout) :: mtab(:) !< all correspondance info.

    integer :: nrec
    integer :: n

    if ( .not. self%opened ) then
      write(0,*)'ERR:jcf_mtf_read_whole_records: File NOT Opened Yet!'
      call exit(1)
    end if

    if ( self%nrec .ne. size(mtab) ) then
      write(0,*)'WRN:jcf_mtf_read_whole_records:maptab size mismatch:'
      write(0,*)'    header:',self%nrec
      write(0,*)'      data:',size(mtab)
    end if

    nrec = min(self%nrec,size(mtab))

    if ( self%binary ) then
      do n = 1, nrec
        read(self%lun)  mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(:)
      end do
    else
      do n = 1, nrec
        read(self%lun,rec_form) mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(:)
      end do
    end if

    self%istat=0
    self%irec=nrec
    return

  end subroutine jcf_mtf_read_whole_records

  !!========================================================================
  !> Read and return one (next) record.
  !!
  !! \note coef(:) must be allocated in advance.
  !! 
  subroutine jcf_mtf_read_one_record(self,ridx,sidx,coef)
    type(maptabfile_t),intent(inout) :: self !< maptabfile_t instance
    integer,intent(out) :: sidx              !< sender index
    integer,intent(out) :: ridx              !< receiver index
    real(dp_k),intent(out):: coef(:)         !< coefficients

    real(dp_k) :: cf(2)
    integer :: ios=0
    character(len=80) :: emsg

    cf(:) = mtf_fmiss
    if ( self%binary ) then
      read(self%lun,iostat=ios) ridx, sidx, coef(:)
    else
      read(self%lun,rec_form,iostat=ios) ridx, sidx, coef(:)
    end if

    self%istat = ios
    self%irec = self%irec+1

    return
  end subroutine jcf_mtf_read_one_record


  !!========================================================================
  !> Write out header of maptabfile.
  !!
  !! All arguments other than `self` overrides a corresponding member
  !! of `self` before write.
  !!
  !! \note
  !! After calling this routine, file point is set the head of first record.
  !! 
  subroutine jcf_mtf_write_head(self, crecv, nrecv, csend, nsend, nrec, ncf )
    type(maptabfile_t),intent(inout) :: self       !< maptab file
    character(len=*), intent(in),optional :: crecv !< name and desc of receiver 
    integer,          intent(in),optional :: nrecv !< num of receiver polygon
    character(len=*), intent(in),optional :: csend !< name and desc of sender
    integer,          intent(in),optional :: nsend !< num of sender polygon
    integer,          intent(in),optional :: nrec  !< num of record
    integer,          intent(in),optional :: ncf   !< num of coefs in one correspondance.

    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_write_head: File NOT Opened Yet!'
      call exit(1)
    end if

    if ( present(crecv) ) self%crecv = crecv
    if ( present(nrecv) ) self%nrecv = nrecv
    if ( present(csend) ) self%csend = csend
    if ( present(nsend) ) self%nsend = nsend
    if ( present(nrec ) ) self%nrec  = nrec
    if ( present(ncf  ) ) self%ncf   = ncf

    if ( self%binary ) then
      rewind(self%lun)
      write(self%lun) self%crecv
      write(self%lun) self%nrecv
      write(self%lun) self%csend
      write(self%lun) self%nsend 
      write(self%lun) self%nrec, self%ncf
    else
      rewind(self%lun)
      write(self%lun,head_form(1)) self%crecv
      write(self%lun,head_form(2)) self%nrecv
      write(self%lun,head_form(3)) self%csend
      write(self%lun,head_form(4)) self%nsend 
      write(self%lun,head_form(5)) self%nrec, self%ncf  
    end if

    self%istat=0
    self%irec=0

    return
  end subroutine jcf_mtf_write_head


  !!========================================================================
  !> Write whole record as an array of maptab_t,
  !!
  subroutine jcf_mtf_write_whole_records(self, mtab)
    type(maptabfile_t),intent(inout) :: self    !< maptabfile_t instance
    type(maptab_t),    intent(inout) :: mtab(:) !< all correspondance info.


    integer :: nrec

    integer :: n

    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_write_whole_records: File NOT Opened Yet!'
      call exit(1)
    end if

    nrec = size(mtab)

    if ( self%binary ) then
      do n = 1, nrec
        write(self%lun) mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(:)
      end do
    else
      do n = 1, nrec
        write(self%lun,rec_form) mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(:)
      end do
    end if

    self%nrec=nrec
    self%irec=nrec

    return

  end subroutine jcf_mtf_write_whole_records



  !!========================================================================
  !> Write one record.
  subroutine jcf_mtf_write_one_record(self,ridx,sidx,coef)
    type(maptabfile_t),intent(inout) :: self !< maptabfile_t
    integer,intent(in) :: sidx               !< sender index
    integer,intent(in) :: ridx               !< receiver index
    real(dp_k),intent(in):: coef(:)           !< coefficients


    integer :: ios=0
    character(len=80) :: emsg

    if ( self%binary ) then
      write(self%lun) ridx, sidx, coef(:)
    else
      write(self%lun,rec_form) ridx, sidx, coef(:)
    end if

    self%istat = 0
    self%irec = self%irec+1

    return
  end subroutine jcf_mtf_write_one_record


  !!========================================================================
  !> Allocate array ob maptab_t.
  !!
  !! `mt` MUST NOT ALLOCATED in advance.
  !! 
  subroutine jcf_mtab_alloc(mt,nrec,ncf)
    type(maptab_t),intent(out),allocatable :: mt(:) !< array of maptab_t.
    integer       ,intent(in)  :: nrec              !< num of correspondence records.
    integer       ,intent(in)  :: ncf               !< num of coefs in one correspondence.

    integer :: n
    allocate( mt(nrec) )

    do n=1, nrec
      allocate( mt(n)%coef(ncf) )
    end do

  end subroutine jcf_mtab_alloc


  !!========================================================================
  !> Deallocate array ob maptab_t
  subroutine jcf_mtab_dealloc(mt)
    type(maptab_t),intent(inout),allocatable :: mt(:) !< array of maptab_t

    integer :: nrec  !< num of correspondence records.
    integer :: n

    nrec = size(mt)

    do n=1, nrec
      deallocate( mt(n)%coef )
    end do

    deallocate( mt )

  end subroutine jcf_mtab_dealloc



  !!========================================================================
  !> Set maptab_t.
  !!
  !! \note
  !! Do not specify intent(out) for `mt`. If so, mt%coef seems to be
  !! deallocated on entry in gfortran ver.4.8.2.
  !!
  subroutine jcf_mtab_set(mt,ridx,sidx,coef)
    type(maptab_t),intent(inout) :: mt      !< one element of maptab_t
    integer       ,intent(in)    :: ridx    !< index of receiver polygon
    integer       ,intent(in)    :: sidx    !< index of sender polygon
    real(kind=8)  ,intent(in)    :: coef(:) !< coefficients.

    integer :: n
    integer :: nn

    if ( size(mt%coef).ne.size(coef)) then
      write(0,*)'WRN:jcf_mtab_set:size mismatch:'
      write(0,*)'    argument:',size(coef)
      write(0,*)'      struct:',size(mt%coef)
    end if

    nn = min(size(mt%coef),size(coef))

    mt%ridx = ridx
    mt%sidx = sidx
    mt%coef(1:nn) = coef(1:nn)
    return

  end subroutine jcf_mtab_set





!!$========================================================================
!!$ Private routine(s)
!!$========================================================================


  !!========================================================================
  !> check if this file is binary or not.
  subroutine check_if_binary_file(res, name)
    character(len=*),intent(in) :: name
    logical,intent(out) :: res

    integer :: lun = 467
    character(len=MTF_CLEN) :: rdesc !! name and description of receiver
    integer :: nn

    integer :: ios,ios1,ios2,ios3
    character(len=80) :: emsg=''

    logical :: lbin,ltxt

    integer :: irec, sidx, ridx
    real(dp_k) :: cf
    integer :: nop !< num_of_operation

    !write(0,*)'dbg:check_if_binary_file:name=',trim(name)

    !! check if text(formatted) file ??
    ios=0
    emsg=''
    open(unit=lun,file=name,form='formatted',status='old',iostat=ios,iomsg=emsg)
!!$    write(0,*)'dbg:check_if_binary_file:open txt,ios,emsg:',ios,trim(emsg)
    
    ios=0
    emsg=''
    read(lun,*,iostat=ios,iomsg=emsg) ! skip desc
    read(lun,*,iostat=ios,iomsg=emsg) nn
!!$    write(0,*)'dbg:check_if_binary_file:read txt:ios,emsg:',ios,trim(emsg)
!!$    write(0,*)'dbg:check_if_binary_file:read txt:nrecv:',nn

    close(lun)

    ltxt = ( ios .eq. 0 )

!!$    write(0,*)'dbg:check_if_binary_file:ltxt:',ltxt

    !! binary ??
    ios=0
    emsg=''
    open(unit=lun,file=name,form='unformatted',status='old',iostat=ios,iomsg=emsg)
!!$    write(0,*)'dbg:check_if_binary_file:open bin:ios,emsg',ios,trim(emsg)

    ios=0
    emsg=''
    read(lun,iostat=ios,iomsg=emsg) !! skip
    read(lun,iostat=ios,iomsg=emsg) nn
!!$    write(0,*)'dbg:check_if_binary_file:read bin:ios,emsg:',ios,trim(emsg)
!!$    write(0,*)'dbg:check_if_binary_file:read bin:nrecv:',nn

    close(lun)

    lbin = ( ios==0 )
!!$    write(0,*)'dbg:check_if_binary_file:lbin:',lbin


    if ( lbin .eqv. ltxt ) then !! something wrong
      write(0,*)'ERR:check_if_binary_file:lbin and ltxt is equivalent',lbin,ltxt
      call exit(1)
    end if

    res = lbin

  end subroutine check_if_binary_file

end module jcf_maptabfile
