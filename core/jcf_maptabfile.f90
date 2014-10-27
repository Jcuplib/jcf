!> module for handling mapping table file.
!!
!! Mapping table file consists of header and data records.
!!
!! Header consists of following records;
!! - character(len=80) :: name and description of receiver compo.
!! - integer :: num of polygons of receiver compo.
!! - character(len=80) :: name and description of sender compo.
!! - integer :: num of polygons of sender compo.
!! - integer :: num of coeffisient records.
!! .
!!
!! Each of data record are as follows in one record:
!! - integer :: sequential id
!! - integer :: index of receiver polygon.
!! - integer :: index of sender polygon.
!! - real(kind=8) :: coefficient 1
!! - read(kind=8) :: coefficient 2
!!
!! Mapping table file may text (formatted) file or binary
!! (unformatted) file.  In each format, contents described above are
!! same.
!!
!! All of public routines are act as like an member function of
!! maptabfile_t, so take mapfiletab_t `self` as first argument.
!! 
!!
!!

module jcf_maptabfile
  use jcf_misc, only:&
       & avail_fid => jcf_avail_fid
  implicit none 

  private !! is default
  public :: maptabfile_t
  public :: maptab_t

  public :: jcf_mtf_init
  public :: jcf_mtf_set
  public :: jcf_mtf_open
  public :: jcf_mtf_close
  public :: jcf_mtf_read_head
  public :: jcf_mtf_read_whole_records
  public :: jcf_mtf_read_one_record
  public :: jcf_mtf_write_head
  public :: jcf_mtf_write_whole_records

  integer,parameter :: dp_k = 8 !! define somewhere.

  integer,parameter,public :: mtf_fnlen = 1024
  integer,parameter,public :: MTF_CLEN = 80

  type maptabfile_t
    private
    integer,public :: istat = 0        !< non-zero if something happens.
    integer :: lun
    character(len=mtf_fnlen) :: name = 'XXXX'
    logical :: binary = .false. !< .F. if text format
    logical :: opened = .false. !< .T. after opened
  end type maptabfile_t


  !> mapping table record.
  !!
  !! should be used as an array of this type.
  type maptab_t
!!$    integer      :: id       !< record id
    integer      :: ridx     !< index of receiver polygon
    integer      :: sidx     !< index of sender polygon
    real(kind=8) :: coef(2)  !< coefficient
  end type maptab_t


  integer,parameter,public :: mtf_imiss =  -999
  real(dp_k),parameter,public :: mtf_fmiss = -999.d0

  character(len=*),parameter :: rec_form='(2I10,2E20.10)'
contains
  !!========================================================================
  !> Initialize maptabfile_t
  subroutine jcf_mtf_init(self)
    type(maptabfile_t),intent(out) :: self

    self=maptabfile_t( &
         & lun  = 0, &
         & name = 'XXXX',&
         & binary=.false.,&
         & opened=.false.,&
         & istat=0 )

    return
  end subroutine jcf_mtf_init


  !!========================================================================
  !> Dump maptabfile_t
  subroutine jcf_mtf_dump(self,unit,head)
    type(maptabfile_t),intent(in) :: self
    integer,intent(in),optional :: unit
    character(len=*),optional :: head
    integer :: lun = 6

    if ( present(unit) ) lun=unit

    if ( present(head) )  then
      write(lun,'(A)')   head
    else
      write(lun,'(A)')   '##### dump maptabfile_t #####'
    end if
    write(lun,'(A,A)') '      name:',trim(self%name)
    write(lun,'(A,L5)')'    binary:',self%binary
    write(lun,'(A,L5)')'    opened:',self%opened
    write(lun,'(A,I5)')'     istat:',self%istat

    return
  end subroutine jcf_mtf_dump


  !!========================================================================
  !> Set some attributes of maptabfile_t in advance.
  !!
  subroutine jcf_mtf_set(self, &
       & name   , &
       & binary   )
    type(maptabfile_t),intent(inout) :: self
    character(len=*),intent(in),optional :: name
    logical,intent(in),optional :: binary

    if ( present(name) ) self%name=name
    if ( present(binary) ) self%binary=binary
       
    return

  end subroutine jcf_mtf_set

  !!========================================================================
  !> Open mappingtable file specified by fname.
  !!
  !! If something happens, self%istat .ne. 0
  !!
  subroutine jcf_mtf_open(self, act, fname)
    type(maptabfile_t),intent(inout) :: self  !< maptabfile.
    character(len=*),  intent(in)    :: act   !< 'R' or 'W'
    character(len=*),  intent(in), optional    :: fname !< name of maptablefile.

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
!!$      write(0,*)'dbg:jcf_mtf_open:lbin:',lbin
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
!!$    write(0,*)'dbg:jcf_mtf_open:lun:',lun
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

    call jcf_mtf_dump(self,unit=0,head='jcf_mtf_open:Opened:')
    return

  end subroutine jcf_mtf_open



  !!========================================================================
  !> close mappingtable file.
  subroutine jcf_mtf_close( self )
    type(maptabfile_t),intent(inout) :: self

    close(self%lun)

    self%name='XXXX'
    self%lun=0
    self%opened=.false.
    self%istat=0

    return
  end subroutine jcf_mtf_close


  !!========================================================================
  !> read header.
  !!
  !! \note
  !! After calling this routine, file point is set the head of first record.
  !! 
  subroutine jcf_mtf_read_head(self, crecv, nrecv, csend, nsend, nrec )
    type(maptabfile_t),intent(inout) :: self  !< maptab file
    character(len=MTF_CLEN),intent(out)  :: crecv !< name and desc of receiver 
    integer,            intent(out)  :: nrecv !< num of receiver polygon
    character(len=MTF_CLEN),intent(out)  :: csend !< name and desc of sender
    integer,            intent(out)  :: nsend !< num of sender polygon
    integer,            intent(out)  :: nrec  !< num of record


    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_read_head: File NOT Opened Yet!'
      call exit(1)
    end if

    if ( self%binary ) then
      rewind(self%lun)
      read(self%lun) crecv
      read(self%lun) nrecv
      read(self%lun) csend
      read(self%lun) nsend 
      read(self%lun) nrec  
    else
      rewind(self%lun)
      read(self%lun,'(A80)') crecv
      read(self%lun,'(I16)'  ) nrecv
      read(self%lun,'(A80)') csend
      read(self%lun,'(I16)'  ) nsend 
      read(self%lun,'(I16)'  ) nrec  
    end if

    self%istat=0
    return
  end subroutine jcf_mtf_read_head



  !!========================================================================
  !> Read whole record and return array of maptab_t,
  !!
  !! Must allocate mtab(:) BEFORE calling me.
  !! 
  subroutine jcf_mtf_read_whole_records(self, mtab, nrec)
    type(maptabfile_t),intent(inout) :: self
    type(maptab_t),intent(out) :: mtab(:)
    integer, intent(in) :: nrec

    integer :: n

    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_read_whole_records: File NOT Opened Yet!'
      call exit(1)
    end if

    if ( self%binary ) then
      do n = 1, nrec
        read(self%lun)  mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(1:2)
      end do
    else
      do n = 1, nrec
        read(self%lun,rec_form) mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(1:2)
      end do
    end if

    return

  end subroutine jcf_mtf_read_whole_records

  !!========================================================================
  !> Read and return one record.
  subroutine jcf_mtf_read_one_record(self,ridx,sidx,coeff1,coeff2)
!!$F2003    class(maptabfile_t),intent(inout) :: self
    type(maptabfile_t),intent(inout) :: self
    integer,intent(out) :: sidx !< sender index
    integer,intent(out) :: ridx !< receiver index
    real(dp_k),intent(out):: coeff1 !< coefficient
    real(dp_k),intent(out),optional :: coeff2 !< second coefficient, if available.

    real(dp_k) :: cf(2)
    integer :: ios=0
    character(len=80) :: emsg

    cf(:) = mtf_fmiss
    if ( self%binary ) then
      read(self%lun,iostat=ios,iomsg=emsg) ridx, sidx, cf(1:2)
    else
      read(self%lun,rec_form,iostat=ios,iomsg=emsg) ridx, sidx, cf(1:2)
    end if

    if ( ios .ne. 0 ) then
      write(0,*)'ERR:jcf_mtf_read_one_record:Read ERRor: '//trim(emsg)
      call exit(1)
    end if

    coeff1 = cf(1)
    if (present(coeff2)) coeff2 = cf(2)


    return
  end subroutine jcf_mtf_read_one_record


  !!========================================================================
  !> write header.
  !!
  !! \note
  !! After calling this routine, file point is set the head of first record.
  !! 
  subroutine jcf_mtf_write_head(self, crecv, nrecv, csend, nsend, nrec )
    type(maptabfile_t),intent(inout) :: self  !< maptab file
    character(len=*),intent(in)   :: crecv !< name and desc of receiver 
    integer,            intent(in)   :: nrecv !< num of receiver polygon
    character(len=*),intent(in)   :: csend !< name and desc of sender
    integer,            intent(in)   :: nsend !< num of sender polygon
    integer,            intent(in)   :: nrec  !< num of record

    character(len=MTF_CLEN) :: cdesc

    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_write_head: File NOT Opened Yet!'
      call exit(1)
    end if

!!$    write(0,*)'dbg:jcf_mtf_write_head:self%lun,self%binary:',self%lun,self%binary
    if ( self%binary ) then
      rewind(self%lun)
      cdesc=trim(crecv)
      write(self%lun) cdesc
      write(self%lun) nrecv
      cdesc=trim(csend)
      write(self%lun) cdesc
      write(self%lun) nsend 
      write(self%lun) nrec  
    else
      rewind(self%lun)
      cdesc=trim(crecv)
      write(self%lun,'(A80)') cdesc
      write(self%lun,'(I16)') nrecv
      cdesc=trim(csend)
      write(self%lun,'(A80)') cdesc
      write(self%lun,'(I16)') nsend 
      write(self%lun,'(I16)') nrec  
    end if

    self%istat=0
    return
  end subroutine jcf_mtf_write_head


  !!========================================================================
  !> Write whole record as an array of maptab_t,
  !!
  subroutine jcf_mtf_write_whole_records(self, mtab, nrec)
    type(maptabfile_t),intent(inout) :: self
    type(maptab_t),intent(out) :: mtab(:)
    integer, intent(in) :: nrec

    integer :: n

    if ( .not. self%opened ) then
      write(0,*)'ERR: jcf_mtf_write_whole_records: File NOT Opened Yet!'
      call exit(1)
    end if

!!$    write(0,*)'dbg:jcf_mtf_write_whole_records:self%lun,self%binary:',self%lun,self%binary
    if ( self%binary ) then
      do n = 1, nrec
        write(self%lun) mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(1:2)
      end do
    else
      do n = 1, nrec
        write(self%lun,rec_form) mtab(n)%ridx, mtab(n)%sidx, mtab(n)%coef(1:2)
      end do
    end if

    return

  end subroutine jcf_mtf_write_whole_records


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
