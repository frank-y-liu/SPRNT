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

## makefile for API

.PHONY: clean realclean

ROOT=..
include ${ROOT}/Make.rules

# locations to find the libraries, hard coded here for the time being
SRCLIB=$(LIBHOME)/lib
SRCINC=$(LIBHOME)/include

# which BLAS to use, this is default from the system
# we can also install a tuned ATLAS 
BLAS = -lblas

# optimzation flags
OPT_FLAG = -O3

# bits
BTS = -m64  ## of course!

# compiler flags
CFLAGS = -Wall -DLINUX ${OPT_FLAG} -I${SRCINC} ${BITS}
CXXFLAGS = -Wall -std=c++0x -DLINUX ${OPT_FLAG} -I${SRCINC} ${BITS}

# linker flag
LINK_FLAG=-L${SRCLIB} -lsprnt \
	-L${SPLINE_HOME} -ldbsplin \
	${BLAS} ${LIBFORTRAN} \
	${BITS} -lm -ldl


##################################################################
##
## should not have to change too much beyond this point
##

## dependency, we simply add all local .C and .c files
HDRS = $(wildcard *.h)
OBJS = $(patsubst %.C,%.o, $(wildcard *.C)) $(patsubst %.c,%.o, $(wildcard *.c))
LIBS = ${SPTLIBS}/libsprnt.a

EXEC = api_demo

all: ${EXEC}

$(OBJS): $(HDRS)
$(EXEC): $(OBJS)
	$(CXX) $^ -o ${EXEC} ${OPT_FLAG} ${LINK_FLAG}

## other options
clean:
	rm -rf $(OBJS)

realclean:
	rm -rf $(OBJS) $(EXEC) 

distclean:
	rm -rf $(OBJS) $(EXEC) 

## end

