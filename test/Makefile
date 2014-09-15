### Makefile for test.

include ../Mkinclude

INCFLAGS = -I../core
LDFLAGS = -L../core
LDLIBS = -ljcf

LIB=../core/libjcf.a

EXE = \
 test_xyz2latlon \
 test_latlon2xyz \
 test_cal_great_circle_area \
 test_get_length \
 test_is_same_point \
 test_is_on_line \
 test_is_same_line

TRASH= *~ *.o *.lst

all: exe
exe: $(EXE)
lib: $(LIB)
	make -C ../core all install

$(EXE): $(LIB)

clean:
	$(RM) $(TRASH) $(EXE)

distclean: clean
	$(RM) $(EXE)


%: %.f90
	$(LINK.F)  $^ $(LDLIBS) -o $@


