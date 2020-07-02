include ../../make.inc
#
#

all: \
cfod.o de2.o deabm.o debdf.o defehl.o derkf.o derkfs.o hstart.o\
intrp.o intyd.o lsod.o pjac.o slvs.o step2.o stod.o vnorm.o vnwrms.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<

