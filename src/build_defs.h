/*********************************************************************
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
*********************************************************************/

//
// nothing to see here, move on
//

#ifndef _VERSION_H
#define _VERSION_H

#include <stdio.h>

#define VERSION_MAJOR 1
#define VERSION_MINOR 3
#define VERSION_REVISION 2

#ifdef __GNUG__

#define BUILD_DATE __DATE__
#define BUILD_TIME __TIME__
#define COMPILER_VER __VERSION__

#else

#define BUILD_DATE "unknown"
#define BUILD_TIME "unknown"
#define COMPILER_VER "unknown"

#endif

void PrintHeader(FILE *F) {
  fprintf(F,"  \n");
  fprintf(F,"   SPRNT: Simulation Program for River Networks\n");
  fprintf(F,"      version %d.%d.%d (%s) : %s  %s\n", 
	  VERSION_MAJOR, VERSION_MINOR, VERSION_REVISION, GIT_HASH, BUILD_DATE, BUILD_TIME);
  fprintf(F,"      compiled with %s\n", COMPILER_VER);
  fprintf(F,"  \n");
}

#undef VERSION_MAJOR
#undef VERSION_MINOR
#undef VERSION_REVISION

#undef BUILD_DATE
#undef BUILD_TIME
#undef COMPILER_VER

#endif

// Local Variables:
// mode: c++
// End:
