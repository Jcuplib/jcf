# -*- makefile -*-  for JCF

### JCF_SYS should be defined as an environmental variable.
# JCF_SYS	= GNU
# JCF_SYS	= Intel
# JCF_SYS	= PGI
# JCF_SYS	= FX10

include ./sysdep/Mkinclude.$(JCF_SYS)

all lib install:
	$(MAKE) -C core $@

clean realclean distclean:
	$(MAKE) -C core $@

distclean: clean-bin

clean-bin:
	$(RM) ./bin/*
	$(RM) ./lib/*
	$(RM) ./include/*

