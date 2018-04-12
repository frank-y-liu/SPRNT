include ../../make.inc
#
#

all: \
bakslv.o chlhsn.o choldc.o d1fcn.o d2fcn.o dfault.o dogdrv.o dogstp.o\
forslv.o fstocd.o fstofd.o grdchk.o heschk.o hookdr.o hookst.o\
hsnint.o lltslv.o lnsrch.o mvmltl.o mvmlts.o mvmltu.o optchk.o\
optdrv.o optif0.o optif9.o optstp.o qraux1.o qraux2.o qrupdt.o\
result.o sclmul.o secfac.o secunf.o sndofd.o tregup.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<

