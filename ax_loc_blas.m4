# ============================================================================
#
# SYNOPSIS
#
#   AX_LOC_BLAS
#
# DESCRIPTION
#
#   Check for BLAS library path, not comprehensive
#
#   site-specification locations can be specified by BLAS_LIB
#

# Frank Liu, Sept 2015,
# 

# LICENSE
#
#   Copyright (c) 2015 Frank Liu <frankliu@us.ibm.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
#

AC_DEFUN([AX_LOC_BLAS],[
ax_loc_blas_ok=no
LOC_BLAS_ARCHIVE=""
# try customer flag first
if test $ax_loc_blas_ok = no; then
   if test "x$BLAS_LIB" != x; then
       saved_LIBS="$LIBS";
       LIBS="$BLAS_LIB $LIBS"
       AC_MSG_CHECKING([check user defined flags specified in BLAS_LIB ...])
       AC_TRY_LINK_FUNC(dgemm_, [ax_loc_blas_ok=yes
                                 LOC_BLAS_LIBS=$BLAS_LIB
                                 LOC_BLAS_ARCHIVE="-Wl,--whole-archive $BLAS_LIB -Wl,--no-whole-archive"
                                 ],[])
       AC_MSG_RESULT($ax_loc_blas_ok)
       LIBS="$saved_LIBS"
   fi   
fi

# try system default
if test $ax_loc_blas_ok = no; then
   AC_CHECK_LIB(blas, dgemm_, [ax_loc_blas_ok=yes
                              LOC_BLAS_LIBS="-lblas"],[
   AC_MSG_ERROR([unable to find any useable blas]);
			      ])
fi
AC_SUBST([LOC_BLAS_LIBS])
AC_SUBST([LOC_BLAS_ARCHIVE])
]) dnl done
