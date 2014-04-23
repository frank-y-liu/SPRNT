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

## $Id$

## makefile for SPRNT library

.PHONY: clean realclean

## where to find Make.rules, in which we set up a bunch of things
ROOT=..
include ${ROOT}/Make.rules

##################################################################
##
## should not have to change anything beyond this point
##

CFLAGS=-Wall -DLINUX ${OPT_FLAG} \
	${BITS} ${MYDEBUG} 

CXXFLAGS=-Wall -DLINUX ${OPT_FLAG} \
	${BITS} ${MYDEBUG} 

ARFLAGS=ruv

HDRS=   $(wildcard *.h)
OBJS=   $(patsubst %.C,%.o, $(wildcard *.C)) $(patsubst %.c,%.o, $(wildcard *.c))

## 
TARGET = libsprnt.a

${TARGET}: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@


# dependencies, blanket
$(OBJS):$(HDRS)

## other options
clean:
	rm -rf $(OBJS)

realclean:
	rm -rf $(OBJS) $(TARGET)

## end







