include ../sysdep/Mkinclude.$(JCF_SYS)

OBJS	= \
 nicam_grid.o \
 io_grid.o \
 mod_logfile.o
SRCS	= $(OBJS:.o=.f90)


INCDIR	= ../include
LIB	= ../lib/libjcf.a

all: objs

objs: $(OBJS)

clean: 
	$(RM) $(TRASH)

realclean: clean

distclean: realclean

nicam_grid.o : mod_logfile.o
