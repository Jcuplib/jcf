# -*- makefile -*-  for JCF

include ./Mkinclude

all: install
install:
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

core:
	$(MAKE) -C $@
app example: install
	$(MAKE) -C $@

.PHONY: core app example
