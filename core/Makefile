include ../Mkinclude

LIB	= libjcf.a

OBJS	= \
 jcf_misc.o \
 jcf_sphere_lib.o \
 jcf_spherical_lib.o \
 jcf_mesh_base.o \
 jcf_coef_base.o \
 jcf_maptabfile.o 

MODS	= $(OBJS:.o=.mod)


all: lib

lib: $(LIB)

clean:
	$(RM) $(TRASH)

realclean: clean
	$(RM) $(LIB)

distclean: realclean


install: $(LIB) $(INCDIR) $(LIBDIR)
	$(INSTALL) $(LIB) $(LIBDIR)
	$(INSTALL) $(MODS) $(INCDIR)

$(LIB): $(OBJS)
	$(AR) $@ $^ && $(RANLIB) $@

$(BINDIR) $(INCDIR) $(LIBDIR):
	mkdir -p $@


jcf_coef_base.o: jcf_mesh_base.o jcf_sphere_lib.o
jcf_mesh_base.o: jcf_sphere_lib.o
jcf_sphere_lib.o : jcf_spherical_lib.o
jcf_maptabfile.o : jcf_misc.o
