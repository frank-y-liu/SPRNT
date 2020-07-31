.TH CMLIB 3F "June 2000"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
\fB/usr/local/lib/libcm.a\fP 
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
To link to the g77 version of CMLIB, use the options
\fB\-L/usr/local/lib -lcm\fP on the \fBg77\fP command.

For example,
.sp
.B g77 -o myprog myprog.f -L/usr/local/lib -lcm
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
\fB/usr/local/lib/libcm.a\fP       g77 version of library
.br
\fB/usr/local/apps/cmlib/doc\fP    documentation directory
.br
\fB/usr/local/apps/cmlib/src\fP    source directory
.sp
.br
\fB/usr/local/bin/cmdoc\fP   documentation retrieval script
.br
\fB/usr/local/bin/cmsrc\fP   source retrieval script
.br
\fB/usr/local/bin/cmsrcx\fP  source+externals retrieval
.SH SEE\ ALSO
\fBg77\fP(1)
