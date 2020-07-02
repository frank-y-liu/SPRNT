include ../../make.inc
#
#

all: \
caxpy.o ccopy.o cdotc.o cdotu.o crotg.o cscal.o csscal.o cswap.o\
dasum.o daxpy.o dcopy.o ddot.o dnrm2.o drot.o drotg.o drotm.o\
drotmg.o dscal.o dsdot.o dswap.o icamax.o idamax.o isamax.o sasum.o\
saxpy.o scasum.o scnrm2.o scopy.o sdot.o sdsdot.o snrm2.o srot.o\
srotg.o srotm.o srotmg.o sscal.o sswap.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<

