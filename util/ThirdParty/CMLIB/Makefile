## make file for the dbsplin spline package

.PHONY:   clean realclean distclean

AR = ar
ARFLAGS = ruv

FC = gfortran
FFLAGS	= -O3

RANLIB = ranlib

TARGET = libdbsplin.a

ALL = TARGET

SRCS=	$(wildcard src/dbsplin/*.f) $(wildcard src/dtensbs/*.f) \
	$(wildcard src/machcon/*.f) $(wildcard src/xerror/*.f)

OBJS=	$(patsubst %.f,%.o, $(wildcard src/dbsplin/*.f)) \
	$(patsubst %.f,%.o, $(wildcard src/dtensbs/*.f)) \
	$(patsubst %.f,%.o, $(wildcard src/machcon/*.f)) \
	$(patsubst %.f,%.o, $(wildcard src/xerror/*.f)) 

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	-rm -rf $(OBJS)

realclean:
	-rm -rf $(OBJS) $(TARGET) 

distclean:
	-rm -rf $(OBJS) $(TARGET) 
	-rm -rf src
	-rm -rf dbsplin_archive
	-rm -rf dtensbs_archive
	-rm -rf dbsplin_depend

# DO NOT DELETE
