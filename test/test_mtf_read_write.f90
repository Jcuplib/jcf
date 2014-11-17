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
  integer,parameter :: ncf = 2

  type(maptabfile_t) :: mtf
  character(len=32) :: fbase='maptabfile_test'
  character(len=32) :: fname

  type(maptab_t),allocatable :: mtab(:)

  integer :: nr,ns,nn,nc
  character(len=CLEN) :: cr, cs

  integer :: ri, si
  real(kind=8),allocatable :: cf(:)


  integer :: n


  call jcf_mtab_alloc(mtab,nrec,ncf)
  do n =1, nrec
    call jcf_mtab_set(mtab(n),ridx=1,sidx=1,coef=(/0.1d0,0.2d0/))
  end do

  write(*,*)'##### Write mapping table(binary) #####'
  fname=trim(fbase)//'.dat'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf    ,&
       & name   = fname   ,&
       & binary = .true. ,&
       & crecv  = crecv   ,&
       & nrecv  = nrecv   ,&
       & csend  = csend   ,&
       & nsend  = nsend   ,&
       & nrec   = nrec    ,&
       & ncf    = ncf     )

  call jcf_mtf_open(mtf,act='Write')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_write_head(mtf)
  write(*,*)'Status after write_head:',mtf%istat
  call jcf_mtf_write_whole_records(mtf, mtab)
  write(*,*)'Status after write_whole_records:',mtf%istat
  call jcf_mtf_close(mtf)


  write(*,*)'##### Write mapping table(text) #####'
  fname=trim(fbase)//'.txt'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf    ,&
       & name   = fname   ,&
       & binary = .false. ,&
       & crecv  = crecv   ,&
       & nrecv  = nrecv   ,&
       & csend  = csend   ,&
       & nsend  = nsend   ,&
       & nrec   = nrec    ,&
       & ncf    = ncf     )

  call jcf_mtf_open(mtf,act='Write')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_write_head(mtf)
  write(*,*)'Status after write_head:',mtf%istat
  call jcf_mtf_write_whole_records(mtf, mtab)
  write(*,*)'Status after write_whole_records:',mtf%istat
  call jcf_mtf_close(mtf)

  call jcf_mtab_dealloc(mtab)


  write(*,*)
  write(*,*)'##### Read whole mapping table(binary) #####'

  fname=trim(fbase)//'.txt'
  call jcf_mtf_set(mtf,name=fname)
  call jcf_mtf_open(mtf,act='Read')
  call jcf_mtf_read_head(mtf,&
       & crecv=cr  ,&
       & nrecv=nr  ,&
       & csend=cs  ,&
       & nsend=ns  ,&
       & nrec=nn   ,&
       & ncf = nc )
  write(*,*)'Read Header:'
  write(*,*)'    crecv:','"'//cr//'"'
  write(*,*)'    nrecv:',nr
  write(*,*)'    csend:','"'//cs//'"'
  write(*,*)'    nsend:',ns
  write(*,*)'    nrec :',nn
  write(*,*)
  write(*,*)'    receiver compo:','"'//cr(1:index(cr,' ')-1)//'"'
  write(*,*)'    receiver desc.:','"'//cr(index(cr,' ')+1:)//'"'
  write(*,*)'    sender compo  :','"'//cs(1:index(cs,' ')-1)//'"'
  write(*,*)'    sender desc.  :','"'//cs(index(cs,' ')+1:)//'"'

  if ( index(crecv, 'COCO') > 0 ) then
    write(*,*)'COCO is a receiver !!'
  else
    write(*,*)'COCO is NOT a receiver!!'
  end if

  call jcf_mtab_alloc(mtab,nn,nc)
  call jcf_mtf_read_whole_records(mtf,mtab)
  write(*,*)'Read records:'
  do n=1,nn
    write(*,*)n, mtab(n)%ridx,mtab(n)%sidx,mtab(n)%coef(:)
  end do
  call jcf_mtf_close(mtf)
  
  write(*,*)
  write(*,*)'##### Read mapping table(binary) one by one #####'
  fname=trim(fbase)//'.txt'
  call jcf_mtf_set(mtf, name=fname)
  call jcf_mtf_open(mtf,act='Read')
  call jcf_mtf_read_head(mtf,&
       & crecv = cr  ,&
       & nrecv = nr  ,&
       & csend = cs  ,&
       & nsend = ns  ,&
       & nrec  = nn  ,&
       & ncf   = nc  )
  write(*,*)'Read records:'
  allocate( cf(nc) )
  do n = 1, nn
    call jcf_mtf_read_one_record(mtf, ri, si, cf)
    write(*,*)n,ri,si,cf
  end do

  call jcf_mtf_close(mtf)


end program test_mtf_read_write
