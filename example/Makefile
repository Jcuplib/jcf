### Makefile for examples.

include ../Mkinclude

EXE	= mkmappingtable
OBJS	= \
 $(EXE).o \
 nicamio_grid.o \
 mod_coef.o \
 nicam_grid.o \
 io_grid.o \
 mod_logfile.o

### Borrow app/ module files.
SRCS_BORROW = nicam_grid.f90 io_grid.f90 mod_logfile.f90

OUTPUTS	= mapping_table_nicam_to_io.txt mk_mappingtable.log


LDLIBS 	= $(LIBDIR)/libjcf.a

all: exe

exe: $(EXE)

$(EXE): $(OBJS)
	$(LINK.f90) $(OBJS) $(LDLIBS) -o $@

clean:
	$(RM) $(TRASH)

realclean: clean
	$(RM) $(EXE) $(OUTPUTS)

distclean: realclean
	$(RM) $(SRCS_BORROW)

$(SRCS_BORROW):
	ln -fs ../app/$@ .


### dependency
mkmappingtable.o: nicamio_grid.o mod_logfile.o
nicamio_grid.o: nicam_grid.o io_grid.o mod_coef.o mod_logfile.o
mod_coef.o:
nicam_grid.o: mod_logfile.o
io_grid.o:
mod_logfile.o:




.PHONY : all exe clean realclean distclean

LINK.f90 = $(FC) $(FFLAGS) $(LDFLAGS) $(TARGET_ARCH)
