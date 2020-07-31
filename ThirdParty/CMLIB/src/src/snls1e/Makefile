include ../../make.inc
#
#

all: \
chkder.o dogleg.o enorm.o fdjac1.o fdjac3.o lmpar.o qform.o qrfac.o\
qrsolv.o r1mpyq.o r1updt.o rwupdt.o scov.o snls1.o snls1e.o snsq.o\
snsqe.o sos.o soseqs.o sossol.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<

