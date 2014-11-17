program test_mtf_open
  use jcf_maptabfile
  implicit none 

  type(maptabfile_t) :: mtf

  character(len=32) :: fbase='maptabfile_test'
  character(len=32) :: fname

  fname=trim(fbase)//'.dat'
  write(*,*)'### Test 1: Open Binary for Write:'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf, binary=.true., name=fname)
  call jcf_mtf_open(mtf,act='Write')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_write_head(mtf, &
       & crecv='Receiver',&
       & nrecv=1         ,&
       & csend='Sender'  ,&
       & nsend=1         ,&
       & nrec=1           )
  write(*,*)'Status after write_head:',mtf%istat
  call jcf_mtf_close(mtf)
  write(*,*)'Status after close:',mtf%istat

  write(*,*)
  write(*,*)'### Test 2: Open Binary for Read:'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf, binary=.true., name=fname)
  call jcf_mtf_open(mtf,act='Read')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_close(mtf)
  write(*,*)'Status after close:',mtf%istat

  write(*,*)
  write(*,*)'### Test 3: Open Text for Write:'
  fname=trim(fbase)//'.txt'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf, binary=.false., name=fname)
  call jcf_mtf_open(mtf,act='Write')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_write_head(mtf, &
       & crecv='Receiver',&
       & nrecv=1         ,&
       & csend='Sender'  ,&
       & nsend=1         ,&
       & nrec=1           )
  write(*,*)'Status after write_head:',mtf%istat
  call jcf_mtf_close(mtf)
  write(*,*)'Status after close:',mtf%istat

  write(*,*)
  write(*,*)'### Test 4: Open Text for Read:'
  call jcf_mtf_init(mtf)
  call jcf_mtf_set(mtf, binary=.false., name=fname)
  call jcf_mtf_open(mtf,act='Read')
  write(*,*)'Status after open:',mtf%istat
  call jcf_mtf_close(mtf)
  write(*,*)'Status after close:',mtf%istat

  call exit(0)

end program test_mtf_open
  
