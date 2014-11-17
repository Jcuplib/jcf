#!/bin/sh

TESTS="\
cal_great_circle_area \
cal_great_circle_center_rect \
get_length \
is_cross_line \
is_in_latlon \
is_inner \
is_on_line \
is_same_line \
is_same_point \
latlon2xyz \
xyz2latlon \
mtf_open_close \
mtf_read_write \
"


for t in $TESTS; do \
  exe=./test_$t;
  ifile=input_$t.txt;
  ofile=output_$t.txt;
  echo -n "Do test for $t" ... ;
  if  [ ! -x $exe ] ; then\
    echo "Cannot execute $exe." && exit 1
  else
    $exe < $ifile > $ofile && echo "done."
  fi
  echo ""
done
