# 
#

.PHONY: clean realclean distclean

all:
	cd SuiteSparse/AMD; make; cd -
	cd SuiteSparse/BTF; make; cd -
	cd SuiteSparse/COLAMD; make; cd -
	cd SuiteSparse/KLU; make; cd -
	cd SuiteSparse/UMFPACK; make; cd -
	make -f makefile.local

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

# Local Variables:
# mode: makefile
# End:

