## make file for the dbsplin spline package

.PHONY:   clean realclean distclean

include ../../Make.rules

ARFLAGS = ruv
FFLAGS	= $(OPT_FLAG)

TARGET = libblas.a

ALL = TARGET

SRCS=	$(wildcard src/BLAS-3.6.0/*.f)

OBJS=	$(patsubst %.f,%.o, $(wildcard src/BLAS-3.6.0/*.f))

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	-rm -rf $(OBJS)

realclean:
	-rm -rf $(OBJS) $(TARGET) 

distclean:
	-rm -rf $(OBJS) $(TARGET) 
	-rm -rf src/*

# DO NOT DELETE
