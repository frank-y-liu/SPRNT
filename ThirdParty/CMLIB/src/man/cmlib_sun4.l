.TH CMLIB 3F "August 1999"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
\fB/usr/local/lib/libcm.a\fP 
.br
\fB/usr/local/lib/libcm8.a\fP
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
Two versions of CMLIB are available:
.LP
.B libcm.a
.in +5
This is the standard version of the library.
.in -5
.sp
.B libcm8.a
.in +5
This version is compiled with \fB\-r8\fP option.  The \fB\-r8\fP option 
requests automatic precision doubling.  It causes single 
precision variables to occupy 8 bytes, and double 
precision  variables to occupy 16 bytes.
.in -5
.LP
To link to the library \fBlibcm.a\fP, use the option \fB\-lcm\fP
on the \fBf77\fP command.
To link to the library \fBlibcm8.a\fP, use the options \fB\-r8 \-lcm8\fP
instead.
For example,
.sp
.B f77 -o myprog myprog.f -lcm
.br
.B f77 -r8 -o myprog myprog.f -lcm8
.sp
.LP
Note that programs compiled with \fB\-r8\fP will \fInot\fP work
correctly when linked with \fBlibcm.a\fP.  Similarly, programs linked
to \fBlibcm8.a\fP \fImust\fP be compiled with \fB\-r8\fP.
Unfortunately, no reasonable error message is generated when a
mismatch occurs.
Also, order is important when listing object and library files
on these statements.
Files containing references to subprograms in libraries should preceed
the library file name (or the \fB\-l\fP) on the command line.
.LP
Both versions of CMLIB are compiled with \fBf77\fP optimization option \fB-O3\fP.
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
.SH FILES
.ta 2i
\fB/usr/local/lib/libcm.a\fP       Sun 4 library
.br
\fB/usr/local/lib/libcm8.a\fP      Sun 4 \fB\-r8\fP library
.br
\fB/usr/local/apps/cmlib/doc\fP    documentation directory
.br
\fB/usr/local/apps/cmlib/src\fP    source directory
.sp
.ta 2i
\fBIn /usr/local/bin\fP ...
.br
\fBcmdoc\fP      documentation retrieval script
.br
\fBcmsrc\fP      source retrieval script
.br
\fBcmsrcx\fP     source+externals retrieval script
.SH SEE\ ALSO
\fBf77\fP(1)
.SH BUGS
.LP
Odrpack may be unsuitable for use on workstations in single precision,
unless the \fB\-r8\fP version is used.
.LP
The special function routines in fnlib cannot deliver results to full machine
precision when compiled with \fB\-r8\fP.  These codes will deliver at most 15
correct digits in single precision, and 31 correct digits in double precision.
