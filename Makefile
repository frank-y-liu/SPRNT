#
#  Copyright (C) 2011, 2014 International Business Machines
#
#  Author:  Frank Liu, IBM
#
#  All Rights Reserved. This program and the accompanying materials
#  are made available under the terms of the Eclipse Public License v1.0
#  which accompanies this distribution, and is available at
#  http://www.eclipse.org/legal/epl-v10.html
#

# 
# available methods:
#
# make dep 		: build 3rd party libraries
# make     		: build library and bin
# make install 		: install, note that default prefix is .
# make uninstall
# make test 		: test bin 
# make clean 		: clear .o files 
# make realclean	: clean .o, .a and bin files
# make distclean	: remove (almost) everything
#

.PHONY: clean realclean distclean dep install uninstall test

ROOT=.
include $(ROOT)/Make.rules

prefix=.

all:
	cd $(SRC_HOME); make -j $(NPROC)
	cd $(SPT_HOME); make -j $(NPROC)

test:
	cd $(TEST_HOME); make 

install:
	install -m 0644 $(SRC_HOME)/base.h ${prefix}/include
	install -m 0644 $(SRC_HOME)/blaspr.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/dspline.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/flexvec.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/matrix.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/sptcomplex.h ${prefix}/include 
	install -m 0644  $(SRC_HOME)/spttimemeas.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/ngraph.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/node.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/options.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/sources.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/subcatch.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/tdarray.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/waveforms.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/xsection.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/xsmodel.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/xstrans.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/solver_def.h ${prefix}/include
	install -m 0644  $(SRC_HOME)/solver_interface.h ${prefix}/include
	install -m 0644 $(SRC_HOME)/libsprnt.a ${prefix}/lib
ifeq (${ARCH},Darwin)
	install -m 0644 $(UF_HOME)/libsolvers.dylib ${prefix}/lib
else
	install -m 0644 $(UF_HOME)/libsolvers.so ${prefix}/lib
endif
	install -m 0755 $(SPT_HOME)/sprnt ${prefix}/bin
	install -m 0755  $(UTL_HOME)/run_spt.sh ${prefix}/bin
	install -m 0755  $(UTL_HOME)/chk_output.sh ${prefix}/bin
	install -m 0755  $(UTL_HOME)/extract_time.sh ${prefix}/bin
	install -m 0755  $(UTL_HOME)/extract_node.sh ${prefix}/bin

uninstall:
	cd ${prefix}/bin; \rm -rf *; cd -
	cd ${prefix}/include; \rm -rf *.h; cd -
	cd ${prefix}/lib; \rm -rf *.a; \rm -rf *.so; \rm -rf *.dylib; cd -

dep:
ifeq (${ARCH},Darwin)
	cd $(BLAS_HOME); ./fetch_blas.sh; make; 
endif
	cd $(SPLINE_HOME); ./fetch_cmlib.sh; make;
	cd $(UF_HOME); ./fetch_ufsparse.sh; make;

clean:
	cd $(SRC_HOME); make clean
	cd $(SPT_HOME); make clean
	cd $(TEST_HOME); make clean

testclean:
	cd $(TEST_HOME); make clean

realclean:
	cd $(SRC_HOME); make realclean;
	cd $(SPT_HOME); make realclean;
ifeq (${ARCH}, Darwin)
	cd $(BLAS_HOME); make realclean;
endif
	cd $(SPLINE_HOME); make realclean;
	cd $(UF_HOME); make realclean;
	cd $(TEST_HOME); make realclean;

distclean:
	cd $(SRC_HOME); make realclean;
	cd $(SPT_HOME); make realclean;
ifeq (${ARCH}, Darwin)
	cd $(BLAS_HOME); make distclean;
endif
	cd $(SPLINE_HOME); make distclean;
	cd $(UF_HOME); make distclean;
	cd $(TEST_HOME); make realclean;
	-rm -rf Make.rules
	-rm -rf Makefile
	-rm -rf UFconfig.mk
	-rm -rf uf_makefile.local
	-rm -rf config.status
	-rm -rf config.log
	-rm -rf aclocal.m4

# Local Variables:
# mode: makefile
# End:

