# -*- makefile -*-  for JCF

include ./Mkinclude

all: install
install:
	$(MAKE) -C core $@

clean realclean distclean:
	$(MAKE) -C core $@
	$(MAKE) -C app $@
	$(MAKE) -C example $@
	$(MAKE) -C test $@

distclean: clean-bin

clean-bin:
	$(RM) ./bin/*
	$(RM) ./lib/*
	$(RM) ./include/*

core:
	$(MAKE) -C $@
app example: install
	$(MAKE) -C $@
test: core
	$(MAKE) -C $@

.PHONY: core app example test
