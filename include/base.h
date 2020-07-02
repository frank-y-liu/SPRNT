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

/********************************************************************
  History
  04/20/11    FYL   created
********************************************************************/

#ifndef	_BASE_H
#define	_BASE_H

#include <assert.h>
#include <stdlib.h>
#include <float.h>

/* needed for communication with drss */
enum {OK=0, WARNING=1, ERROR=-1};

/* Internal behavior control */
enum {MAXFIXSZ=5000, MAXDENSESIZE=255, MAXGR=1};

/* X section types */
enum XsecType { RECT=0, TRAP, SPLINE, INTRINSIC, XS_UNKNOWN};

/* communicate with the solver */
enum SolverType { S_UMF, S_KLU, S_UNKNOWN };

/* type of the mass or momentum equations */
enum MMType { A_SRC, A_DEP, MAS_STV, Q_SRC, Q_DEP, MOM_STV, M_UNKNOWN };

/* type of the independent sources */
enum SourceType { SR_PWL, SR_UNKNOWN };

/* default size of the source length */
enum { DRSS_A_DFT=3, DRSS_DEP_DFT=5, DRSS_PWL_DFT=25, DRSS_SRC_DFT=100, DRSS_NN_DFT=500 };

#define TUNIT 1.0  // 1 second is the minimal time when processing Q sources
/*
  Common Macros, inline function type, not macros
*/
template<class T> inline T MIN(T x, T y)       { return(x < y ? x : y); };
template<class T> inline T MIN(T x, T y, T z)  { return MIN(MIN(x,y),z); };
template<class T> inline T MAX(T x, T y)       { return(x < y ? y : x); };
template<class T> inline T MAX(T x, T y, T z)  { return MAX(MAX(x,y),z); };
template<class T> inline void SWAP(T& x, T& y) { T v=x; x=y; y=v; };
template<class T> inline T ABS(T x)            { return (x>=0) ? x : -x; };
template<class T> inline int INT_SIGN(T x)     { return (x>=0) ? 1 : -1; };
template<class T> inline double DBL_SIGN(T x)  { return (x>=0) ? 1.0 : -1.0; };

/* constants have been moved to options */

/* flags to control what to print */
#define PRT_FR 0x20 //Added by JustinYu
#define PRT_XY 0x10
#define PRT_Q  0x08
#define PRT_A  0x04
#define PRT_D  0x02
#define PRT_Z  0x01

/* deal with file compresion */
#ifdef HAVE_LIBZ

#include <zlib.h>
#define sptFile gzFile
#define sptOpenToRead(a) gzopen((a),"r")
#define sptGets(buf,len, fp) gzgets(fp, buf, len)
#define sptRewind(a) gzrewind(a)
#define sptOpenToWrite(a) gzopen((a),"w")
#define sptFprintf(a, ...) gzprintf(a, ##__VA_ARGS__)
#define sptClose(a) gzclose(a)

#else

#define sptFile FILE*
#define sptOpenToRead(a) fopen((a),"r")
#define sptGets(buf,len, fp) fgets(buf, len, fp)
#define sptRewind(a) rewind(a)
#define sptOpenToWrite(a) fopen((a),"w")
#define sptFprintf(a, ...) fprintf(a,__VA_ARGS__)
#define sptClose(a) fclose(a)

#endif

#endif // _BASE_H

// Local Variables:
// mode: c++
// End:
