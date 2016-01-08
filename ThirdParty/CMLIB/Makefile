## make file for the dbsplin spline package

.PHONY:   clean realclean distclean

include ../../Make.rules

ARFLAGS = ruv
FFLAGS	= $(OPT_FLAG)

TARGET = libdbsplin.a

ALL = TARGET

SRCS=	$(wildcard src/src/dbsplin/*.f) $(wildcard src/src/dtensbs/*.f) \
	$(wildcard src/src/machcon/*.f) $(wildcard src/src/xerror/*.f)

OBJS=	$(patsubst %.f,%.o, $(wildcard src/src/dbsplin/*.f)) \
	$(patsubst %.f,%.o, $(wildcard src/src/dtensbs/*.f)) \
	$(patsubst %.f,%.o, $(wildcard src/src/machcon/*.f)) \
	$(patsubst %.f,%.o, $(wildcard src/src/xerror/*.f)) 

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
