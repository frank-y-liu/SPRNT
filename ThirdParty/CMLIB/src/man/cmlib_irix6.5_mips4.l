.TH CMLIB 3F "June 2000"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
\fB/usr/local/lib64/libcm.a\fP 
.br
\fB/usr/local/lib64/libcm_i8r8d16.a\fP 
.br
\fB/usr/local/lib64/libcm_f90.a\fP 
.br
\fB/usr/local/lib64/libcm_f90_default64.a\fP 
.br
\fB/usr/local/lib32/libcm.a\fP 
.br
\fB/usr/local/lib32/libcm_i8r8d16.a\fP 
.br
\fB/usr/local/lib32/libcm_f90.a\fP 
.br
\fB/usr/local/lib32/libcm_f90_i8r8d16.a\fP 
.LP
\fBcmdoc\fP [ \fB\-s\fP [ \fIsublibrary\fP ] ] .\|.\|. [ \fIsubprogram ] .\|.\|.
.br
\fBcmsrc\fP \fIsubprogram\fP .\|.\|.
.br
\fBcmsrcx\fP \fIsubprogram\fP .\|.\|.
.SH DESCRIPTION
.LP
CMLIB is the NIST Core Math Library, a collection of Fortran-callable
subprograms for solving a wide range of mathematical problems. 
For IRIX 6 systems eight implementations are currently supported, 
.LP
.ta 3i
-lcm			f77 -64
.br
-lcm_i8r8d16		f77 -64 -i8 -r8 -d16
.br
-lcm_f90		f90 -64
.br
-lcm_f90_default64	f90 -default64
.br
-lcm			f77 -n32 -mips4
.br
-lcm_i8r8d16		f77 -n32 -i8 -r8 -d16 -mips4
.br
-lcm_f90		f90 -n32 -mips4
.br
-lcm_f90_i8r8d16	f90 -n32 -i8 -r8 -d16 -mips4
.LP
In order to correctly link to VFNLIB, user code must be compiled
with the compiler and options indicated in the second column of
the table above.  For example, each of the following is legal.
.sp
.br
f77 -64 prog.f -lcm
.sp
f77 -64 -i8 -r8 -d16 prog.f -lcm_i8r8d16
.sp
f90 -64 prog.f -lcm_f90
.sp
f90 -default64 prog.f -lcm_f90_default64
.sp
f77 -n32 -mips4 prog.f -lcm
.sp
f77 -n32 -mips4 -i8 -r8 -d16 prog.f -lcm_i8r8d16
.sp
f90 -n32 -mips4 prog.f -lcm_f90
.sp
f90 -n32 -mips4 -i8 -r8 -d16 prog.f -lcm_f90_i8r8d16
.LP
Documentation for  CMLIB comes in three forms:  
(a) an overall summary,
(b) sublibrary summaries, and 
(c) detailed documentation on individual subprograms.  
Sublibraries are subsets of the library which may be
considered as separate entities, such as the EISPACK or FISHPAK
packages.  The overall summary contains a list of sublibraries and
their purposes.  The command \fBcmdoc\fP is used to retrieve this
documentation; it finds documentation on the items selected and
deposits it on standard output.  For example,
.LP
.ta 2i
\fBcmdoc \-s\fP			gets the overall summary
.br
\fBcmdoc \-s eispack\fP		gets the EISPACK summary 
.br
\fBcmdoc sgefa sgesl\fP		gets details on SGEFA and SGESL.
.LP
To discover which CMLIB routine is useful for your application, refer
to the NIST Guide to Available Mathematical Software (GAMS).  GAMS
can be accessed by pointing your Web browser at http://math.nist.gov/.
.LP
Two additional commands are available to retrieve the Fortran source
for CMLIB subprograms.  \fBcmsrc\fP retrieves the Fortran source for
each subprogram whose name is listed on the command line.  \fBcmsrcx\fP
retrieves the same source as \fBcmsrc\fP but also includes the source
for all external references within the library which are required to
use the named subprograms.
.SH FILES
.ta 3i
In \fB/usr/local/lib64\fP ...
.br
\fBlibcm.a\fP			SGI IRIX 6 R10xxx library
.br
\fBlibcm_i8r8d16.a\fP		SGI IRIX 6 R10xxx library
.br
\fBlibcm_f90.a\fP		SGI IRIX 6 R10xxx library
.br
\fBlibcm_f90_default64.a\fP	SGI IRIX 6 R10xxx library
.sp
In \fB/usr/local/lib32\fP ...
.br
\fBlibcm.a\fP			SGI IRIX 6 R10xxx library
.br
\fBlibcm_i8r8d16.a\fP		SGI IRIX 6 R10xxx library
.br
\fBlibcm_f90.a\fP		SGI IRIX 6 R10xxx library
.br
\fBlibcm_f90_i8r8d16.a\fP	SGI IRIX 6 R10xxx library
.sp
\fB/usr/local/apps/cmlib/doc\fP	documentation directory
.br
\fB/usr/local/apps/cmlib/src\fP	source directory
.br
\fB/usr/local/bin/cmdoc\fP	documentation retrieval script
.br
\fB/usr/local/bin/cmsrc\fP	source retrieval script
.br
\fB/usr/local/bin/cmsrcx\fP	source+externals retrieval 
.SH SEE\ ALSO
\fBf77\fP(1), \fBf90\fP(1)
