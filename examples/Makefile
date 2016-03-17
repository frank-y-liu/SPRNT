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

## makefile for testing 

.PHONY: test

ROOT=..
include ${ROOT}/Make.rules

# locations to find the executable 
CMD=LD_LIBRARY_PATH=../ThirdParty/UF ${SPT_HOME}/sprnt
DIFF=diff -q
UNZIP=gzip -q -d

CASES= $(wildcard *.spt)
OUTPUT= $(addsuffix .output.dat, $(wildcard *.spt))

# RULES
%.spt.output.dat : %.spt
	$(CMD) $?
	-$(UNZIP) $?.output.dat.gz

all: test

test : $(OUTPUT)
	$(foreach f,$^, $(DIFF) $(f) golden/$(f); )

clean:
	-rm -rf *.output.dat
	-rm -rf *.output.dat.gz

realclean:
	-rm -rf *.output.dat
	-rm -rf *.output.dat.gz
## end

