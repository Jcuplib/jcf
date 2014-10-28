program test_mtf_read_write
  use jcf_maptabfile
  use jcf_maptabfile,only:&
       & clen=>mtf_clen, &
       & fmiss=>mtf_fmiss
  implicit none 

  character(len=*),parameter :: crecv='COCO Tri-polar 360x256 for test'
  character(len=*),parameter :: csend='NICAM GL05RL01 for test'

  integer,parameter :: nrecv = 1
  integer,parameter :: nsend = 1
  integer,parameter :: nrec = 5

  type(maptabfile_t) :: mtf
  character(len=32) :: fbase='maptabfile_test'
  character(len=32) :: fname

  type(maptab_t),allocatable :: mtab(:)

  integer :: nr,ns,nn
  character(len=CLEN) :: cr, cs

  integer :: ri, si
  real(kind=8) :: cf1,cf2


  integer :: n


  allocate(mtab(nrec))
  do n =1, nrec
    call jcf_mtab_set(mtab(n),ridx=1,sidx=1,coef=(/1.d0,fmiss/))
  end do

  write(*,*)'##### Write mapping table #####'

  fname=trim(fbase)//'.txt'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf, binary=.false., name=fname)
  call jcf_mtf_open(mtf,act='Write')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_write_head(mtf, &
       & crecv=crecv ,&
       & nrecv=nrecv ,&
       & csend=csend ,&
       & nsend=nsend ,&
       & nrec=nrec    )
  write(*,*)'Status after write_head:',mtf%istat

  call jcf_mtf_write_whole_records(mtf, mtab)

  call jcf_mtf_close(mtf)
  deallocate(mtab)


  write(*,*)
  write(*,*)'##### Read whole mapping table #####'

  call jcf_mtf_set(mtf, binary=.false., name=fname)
  call jcf_mtf_open(mtf,act='Read')
  call jcf_mtf_read_head(mtf,&
       & crecv=cr  ,&
       & nrecv=nr  ,&
       & csend=cs  ,&
       & nsend=ns  ,&
       & nrec=nn    )
  write(*,*)'Read Header:'
  write(*,*)'    crecv:',cr
  write(*,*)'    nrecv:',nr
  write(*,*)'    csend:',cs
  write(*,*)'    nsend:',ns
  write(*,*)'    nrec :',nn

  allocate(mtab(nn))

  call jcf_mtf_read_whole_records(mtf,mtab)

  write(*,*)'Read records:'
  do n=1,nn
    write(*,*)n, mtab(n)%ridx,mtab(n)%sidx,mtab(n)%coef(:)
  end do

  call jcf_mtf_close(mtf)
  
  write(*,*)
  write(*,*)'##### Read mapping table one by one #####'
  call jcf_mtf_set(mtf, binary=.false., name=fname)
  call jcf_mtf_open(mtf,act='Read')
  call jcf_mtf_read_head(mtf,&
       & crecv=cr  ,&
       & nrecv=nr  ,&
       & csend=cs  ,&
       & nsend=ns  ,&
       & nrec=nn    )

  write(*,*)'Read records:'
  do n = 1, nn
    call jcf_mtf_read_one_record(mtf, ri, si, cf1, cf2)
    write(*,*)n,ri,si,cf1,cf2
  end do

  call jcf_mtf_close(mtf)


end program test_mtf_read_write
