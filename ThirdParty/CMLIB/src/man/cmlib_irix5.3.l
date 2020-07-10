.TH CMLIB 3F "September 1999"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
\fB/usr/local/lib/libcm.a\fP 
.br
\fB/usr/local/lib/libcm_mips2.a\fP
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
Two versions of CMLIB are implemented for IRIX 5 systems:
.LP
.B libcm.a
.in +5
This is the standard version of the library.
.in -5
.sp
.B libcm_mips2.a
.in +5
This version is compiled with \fB\-mips2\fP option and will only run
on IRIX 5 systems with R4xxx CPUs.  The use of this option can 
generate code that runs 10-50% faster.
.in -5
.LP
To link to the library \fBlibcm.a\fP, use the option \fB\-lcm\fP
on the \fBf77\fP command.
To link to the library \fBlibcm_mips2.a\fP, use the options \fB\-mips2 \-lcm_mips2\fP
instead.
For example,
.sp
.B f77 -o myprog myprog.f -lcm
.br
.B f77 -mips2 -o myprog myprog.f -lcm_mips2
.sp
.LP
Note that programs compiled with \fB\-mips2\fP will \fInot\fP work
correctly when linked with \fBlibcm.a\fP.  
.LP
Both versions of CMLIB are also compiled with \fBf77\fP options \fB-O2 -static\fP.
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
.br
\fB/usr/local/lib/libcm.a\fP       SGI IRIX 5 R3xxx library
.br
\fB/usr/local/lib/libcm_mips2.a\fP SGI IRIX 5 R4xxx library
.sp
.br
\fB/usr/local/bin/cmdoc\fP         documentation retrieval script
.br
\fB/usr/local/bin/cmsrc\fP         source retrieval script
.br
\fB/usr/local/bin/cmsrcx\fP        source+externals retrieval
.sp
\fB/usr/local/apps/cmlib/doc\fP    documentation directory
.br
\fB/usr/local/apps/cmlib/src\fP    source directory
.SH SEE\ ALSO
\fBf77\fP(1)
