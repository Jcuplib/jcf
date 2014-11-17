include ../Mkinclude

OBJS	= \
 nicam_grid.o \
 io_grid.o \
 coco_grid.o
SRCS	= $(OBJS:.o=.f90)


INCDIR	= ../include
LIB	= ../lib/libjcf.a

all: objs

objs: $(OBJS)

clean: 
	$(RM) $(TRASH)

realclean: clean

distclean: realclean
