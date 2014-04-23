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

## makefile for SPRNT 

.PHONY: clean realclean

# where to find Make.rules, where the customerization is done
ROOT=..
include ${ROOT}/Make.rules

##################################################################
##
## should not have to change too much beyond this point
##

CFLAGS=-Wall -DLINUX ${OPT_FLAG} -I${SRC_HOME} ${BITS}
CXXFLAGS=-Wall -std=c++0x -DLINUX ${OPT_FLAG} -I${SRC_HOME} ${BITS}

LINK_FLAGS=-L${SRC_HOME} -lsprnt \
	-L${SPLINE_HOME} -ldbsplin \
	${BLAS} ${LIBFORTRAN} \
	${BITS} -lm -ldl

## 
HDRS = $(wildcard *.h)
OBJS = $(patsubst %.C,%.o, $(wildcard *.C)) $(patsubst %.c,%.o, $(wildcard *.c))
LIBS = ${SRC_HOME}/libsprnt.a ${SPLINE_HOME}/libdbsplin.a

EXEC = sprnt

# dependencies
all: $(EXEC)

$(OBJS):$(HDRS)
$(EXEC):${OBJS} ${LIBS}
	$(CXX) ${OBJS} -o ${EXEC} ${OPT_FLAG} ${LINK_FLAGS}

## other options
clean:
	rm -rf $(OBJS)

realclean:
	rm -rf $(OBJS) $(EXEC) 

## end

