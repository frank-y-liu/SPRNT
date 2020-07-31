.TH CMLIB 3F "June 2000"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
\fB/usr/local/lib64/libcm.a\fP 
.br
\fB/usr/local/lib64/libcm_i8r8d16.a\fP 
.sp
.br
\fB/usr/local/lib32/libcm.a\fP 
.br
\fB/usr/local/lib32/libcm_i8r8d16.a\fP 
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
.LP
Five implementations are currently supported for SGI systems
.sp
Library		compiled with	
.sp
\fBlibcm.a\fP          \fBf77 -64 -mips4 -O2 -static \fP
.br
\fBlibcm_i8r9d16.a\fP  \fBf77 -64 -mips4 -i8 -r8 -d16 -O2 -static\fP
.br
\fBlibcm.a\fP          \fBf77 -n32 -mips3 -O2 -static\fP
.br
\fBlibcm_i8r9d16.a\fP  \fBf77 -n32 -mips3 -i8 -r8 -d16 -O2 -static\fP
.LP
To link to the library \fBlibcm.a\fP, use the \fBf77\fP or \fBf90\fP command.
You must link to a library that is compatible with the options that you used
to compile the calling program.  The possible combinations are
For example,
.sp
.br
 f77 -64 -mips4 -o myprog myprog.f -lcm
.br
 f77 -64 -mips4 -i8 -r8 -d16 -o myprog myprog.f -lcm_i8r8d16
.sp
.br
 f77 -n32 -mips3 -o myprog myprog.f -lcm
.br
 f77 -n32 -mips3 -i8 -r8 -d16 -o prog prog.f -lcm_i8r8d16
.sp
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
Two additional commands are available to retrieve the Fortran source
for CMLIB subprograms.  \fBcmsrc\fP retrieves the Fortran source for
each subprogram whose name is listed on the command line.  \fBcmsrcx\fP
retrieves the same source as \fBcmsrc\fP but also includes the source
for all external references within the library which are required to
use the named subprograms.
.LP
A problem-oriented index to the routines in CMLIB is available in the
NIST Guide to Available Mathematical Software (GAMS).  GAMS can be
accessed by poiting your Web browser to \fBhttp://math.nist.gov/gams/\fP.
.SH FILES
\fB/usr/local/lib64/libcm.a\fP         IRIX mips4 f77 library
.br
\fB/usr/local/lib64/libcm_i8r8d16.a\fP IRIX mips4 f77 library
.br
\fB/usr/local/lib32/libcm.a\fP         IRIX n32 mips3 f77 library
.br
\fB/usr/local/lib32/libcm_i8r8d16.a\fP IRIX n32 mips3 f77 library
.br
\fB/usr/local/apps/cmlib/doc\fP        documentation directory
.br
\fB/usr/local/apps/cmlib/src\fP        source directory
.br
\fB/usr/local/bin/cmdoc\fP       documentation retrieval script
.br
\fB/usr/local/bin/cmsrc\fP       source retrieval script
.br
\fB/usr/local/bin/cmsrcx\fP      source+externals retrieval script
.SH SEE\ ALSO
\fBf77\fP(1)
