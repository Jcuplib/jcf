# -*- makefile -*-  for JCF

include ./Mkinclude

all lib install:
	$(MAKE) -C core $@

clean realclean distclean:
	$(MAKE) -C core $@
	$(MAKE) -C app $@
	$(MAKE) -C example $@

distclean: clean-bin

clean-bin:
	$(RM) ./bin/*
	$(RM) ./lib/*
	$(RM) ./include/*

