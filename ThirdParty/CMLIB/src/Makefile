
# ================================
# Makefile for generic Unix system
# ================================

#  Set the macro FORTRAN  to the compiler you would like to use
#  Set the macro OPTS     to the compiler options you would like to use
#  Set the macro NOPT     for problem routines - currently set to -O
#  Set the macro RANLIB   to the processor that writes tables-of-contents 
#                         on archive files.
#  Set the macro MACHCON  to the machine type (see the file 
#                         ../../utils/r1mach.template)
#  Set the macro INPT     to the default tests to run

include	make.inc

# --------------------------------------------------------------------------

SRCDIR = ../../src
SCRIPTDIR = ../../utils

SUBLIB = \
   adapt   amoslib blas    bocls   bspline bvsup   cdriv  cluster \
   cpqr79  cpzero  dbsplin ddassl  ddriv   depac   dnl2sn dqrlss \
   dtensbs eispack fc      fcnpak  fftpkg  fishpak fnlib  licepak \
   lindrv  linpakc linpakd linpaks lotps   mxentrp nl2sn  odrpack \
   pchipd  pchips  q1da    quadpkd quadpks rv      sdassl sdriv \
   sglss   slrpack slvblk  snls1e  splp    sqrlss  ssort  sudsods \
   tensbs  twodq   uncmin  vfft    vhs3    vsfft   xblas  xerror \
   zeroin
# --------------------------------------------------------------------------

all: libcm testall

libcm: $(SUBLIB) constants $(SPECIALS)
	@echo " " ; echo "Updating library file ..." ; echo " "
	csh $(SCRIPTDIR)/Update $(LIBRARY) '$(RANLIB)'

# --------------------------------------------------------------------------

# ... special replacement sources for this machine

constants : src/machcon/i1mach.f

src/machcon/i1mach.f:
	if [ ! -d src/machcon ];then mkdir -p src/machcon ;fi
	csh $(SCRIPTDIR)/genmachcon $(MACHCON)
	mv i1mach.f r1mach.f d1mach.f src/machcon
	$(SCRIPTDIR)/createmk


# --------------------------------------------------------------------------

#  ... the following targets run the cmlib quickchecks
#  the output is written to both stdout and to a file in
#  ARCH/LIBDIR 
#

testall :
	@echo " "
	@echo "Testing all of cmlib (this make take a while) ..."
	@echo " "
	csh $(SCRIPTDIR)/Testall "$(FORTRAN)" "$(NOPTS)" "" "$(INPT)" "$(LIBRARY)"

#
#  runs a test for a given sublibrary 
#
test :
	@csh $(SCRIPTDIR)/Test "$(PGM)" "$(FORTRAN)" "$(NOPTS)" "" "$(INPT)" "$(LIBRARY)"

# --------------------------------------------------------------------------

#
# clean removes all .o files for all architectures
# and recreates the Makefiles in the src directory for cmlib
# It does not remove libcm.a, the test result output files
# or the special source versions for an architecture
#
clean :
	rm -f src/[a-e]*/*.o
	rm -f src/[f-l]*/*.o
	rm -f src/[m-z]*/*.o
	rm -f src/*/Makefile
	rm -f LOCK job.cmd make.inc job.cmd.tmp
	utils/crmkall
	rm -f [A,H,I,L,S]*/lib*/make.inc
	rm -f [A,H,I,L,S]*/lib*/*/*/*.o
	rm -f [A,H,I,L,S]*/lib*/*.o
	rm -f [A,H,I,L,S]*/lib*/core
#
# veryclean removes all versions of cmlib for all architectures
# this includes all special source versions for an architecture
#
veryclean :
	rm -f src/*/Makefile
	rm -f WorkingPlatform job.cmd make.inc LOCK
	utils/crmkall
	rm -r -f [A,H,I,L,S]*

$(SUBLIB):
	( cd $(SRCDIR)/$@ ; $(MAKE) )
