# 
#

.PHONY: clean realclean distclean

all:
	cd AMD; make; cd -
	cd BTF; make; cd -
	cd COLAMD; make; cd -
	cd KLU; make; cd -
	cd UMFPACK; make; cd -
	make -f makefile.local

clean:
	cd AMD; make clean; cd -
	cd BTF; make clean; cd -
	cd COLAMD; make clean; cd -
	cd KLU; make clean; cd -
	cd UMFPACK; make clean; cd -
	make clean -f makefile.local

realclean:
	cd AMD; make purge; cd -
	cd BTF; make purge; cd -
	cd COLAMD; make purge; cd -
	cd KLU; make purge; cd -
	cd UMFPACK; make purge; cd -
	make realclean -f makefile.local

distclean:
	cd AMD; make purge; cd -
	cd BTF; make purge; cd -
	cd COLAMD; make purge; cd -
	cd KLU; make purge; cd -
	cd UMFPACK; make purge; cd -
	make realclean -f makefile.local
	-rm -rf AMD BTF KLU UMFPACK UFconfig COLAMD
	\rm -rf *.tar.gz

# Local Variables:
# mode: makefile
# End:

