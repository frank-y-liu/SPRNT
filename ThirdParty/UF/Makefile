# 
#

.PHONY: clean realclean distclean

include ../../Make.rules

all:
	cd SuiteSparse/AMD; make -j$(NPROC); cd -
	cd SuiteSparse/BTF; make -j$(NPROC); cd -
	cd SuiteSparse/COLAMD; make -j$(NPROC); cd -
	cd SuiteSparse/KLU; make -j$(NPROC); cd -
	cd SuiteSparse/UMFPACK; make -j$(NPROC); cd -
	make -j$(NPROC) -f makefile.local

clean:
	cd SuiteSparse/AMD; make clean; cd -
	cd SuiteSparse/BTF; make clean; cd -
	cd SuiteSparse/COLAMD; make clean; cd -
	cd SuiteSparse/KLU; make clean; cd -
	cd SuiteSparse/UMFPACK; make clean; cd -
	make clean -f makefile.local

realclean:
	cd SuiteSparse/AMD; make purge; cd -
	cd SuiteSparse/BTF; make purge; cd -
	cd SuiteSparse/COLAMD; make purge; cd -
	cd SuiteSparse/KLU; make purge; cd -
	cd SuiteSparse/UMFPACK; make purge; cd -
	make realclean -f makefile.local

distclean:
	cd SuiteSparse/AMD; make purge; cd -
	cd SuiteSparse/BTF; make purge; cd -
	cd SuiteSparse/COLAMD; make purge; cd -
	cd SuiteSparse/KLU; make purge; cd -
	cd SuiteSparse/UMFPACK; make purge; cd -
	make realclean -f makefile.local
	-rm -rf SuiteSparse
	\rm -rf *.tar.gz
	-rm -rf makefile.local

# Local Variables:
# mode: makefile
# End:

