# A few checks to make sure we have required components
# 
# The configure file will not be removed and should be fairly portable
# To create new configure file, do the following:
#    aclocal -I . 
#    autoconf
#
# use -I flag to make sure that the local m4 will be used.

AC_INIT(sprint_config, version_0.9)
AC_MSG_NOTICE([Checking required libraries for SPRNT])

AC_PREFIX_DEFAULT([.])

# need compilers
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77

# of course
AC_CHECK_TOOL(MAKEB, make, [
    AC_MSG_ERROR([need make to continue...])
])

# need library utilies
AC_CHECK_TOOL(RANLIB, ranlib, [
    AC_MSG_ERROR([need ranlib to continue ...])
])

AC_CHECK_TOOL(AR, ar, [
    AC_MSG_ERROR([need ar to continue ...])
])

# need curl to fetch dependencies
AC_CHECK_PROG(CURL_CHECK,curl,yes)
if test x"$CURL_CHECK" != x"yes" ; then
    AC_MSG_ERROR([need curl to continue ...])
fi

# need dirname for shell scripts
AC_CHECK_PROG(DIRNAME_CHECK,dirname,yes)
if test x"$DIRNAME_CHECK" != x"yes" ; then
    AC_MSG_ERROR([need dirname to continue ...])
fi

# need readlink for shell scripts
AC_CHECK_PROG(REALPATH_CHECK,readlink,yes)
if test x"$REALPATH_CHECK" != x"yes" ; then
    AC_MSG_ERROR([need realpath to continue ...])
fi

# math library should be here, just in case 
AC_CHECK_LIB(m,cbrt,[],[
   AC_MSG_ERROR([something is wrong with the math lib -m ...]);
])

# check zlib, define a local flag and handle it internally
AX_CHECK_ZLIB([
AC_MSG_NOTICE([Found zlib. I/O compression will be enabled])
AC_SUBST(LOC_HAVE_LIBZ,[-DHAVE_LIBZ])
],[
AC_MSG_NOTICE([Unable to find zlib.])
])

# check libgfortran, should be installed with gfortran
AC_CHECK_LIB(gfortran,_gfortran_fget,[],[
   AC_MSG_ERROR([unable to locate libgfortran]);
])

## check blas, determines blas flag
AX_LOC_BLAS

# check if c++0x is available
#  if yes, use <unordered_map> for string->integer mapping
#  if not, use <map> (about twice as slow as <unordered_map>)
AX_LOC_CXX_COMPILE_STDCXX_0X

AC_CONFIG_FILES([
Make.rules
Makefile
UFconfig.mk
uf_makefile.local
])

AC_OUTPUT

AC_MSG_NOTICE([==============================================================])
AC_MSG_NOTICE([Required libraries seem to be in place ...... ])
AC_MSG_NOTICE([Build the thirdparty libraries before])
AC_MSG_NOTICE([building SPRNT. Steps: ])
AC_MSG_NOTICE([ make dep])
AC_MSG_NOTICE([ make ])
AC_MSG_NOTICE([ make test (optional)])
AC_MSG_NOTICE([ make install])
AC_MSG_NOTICE([Options to clean])
AC_MSG_NOTICE([ make uninstall - as the name says])
AC_MSG_NOTICE([ make clean     - clean the libraries only ])
AC_MSG_NOTICE([ make realclean - linbraries + binary])
AC_MSG_NOTICE([ make distclean - remove downloadd 3rd party components and Make files])
AC_MSG_NOTICE([==============================================================])

## we are done
