include ../../make.inc
#
#

all: \
bndacc.o bndsol.o bsplvd.o bsplvn.o bvalue.o cv.o diff.o efc.o\
efcmn.o fc.o fcmn.o h12.o hfti.o interv.o lpdp.o lsei.o lsi.o\
wnlit.o wnlsm.o wnnls.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<

