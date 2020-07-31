.TH CMLIB 3F "20 Sep 99"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
/usr/local/lib/libcm.a 

\fBcmdoc\fP [ \fB\-s\fP [ sublibrary ] ] .\|.\|. [ subprogram ] .\|.\|.
.br
\fBcmsrc\fP subprogram .\|.\|.
.br
\fBcmsrcx\fP subprogram .\|.\|.
.SH DESCRIPTION

.I CMLIB
is the NIST Core Math Library, a collection of Fortran-callable
subprograms for solving a wide range of mathematical problems. 

To link to CMLIB, use the option \fB\-l\fIcm\fP on the
.I f77
or 
.I ld
statement.
On the HP9000, CMLIB was compiled with f77 options -O -K.

Documentation for 
.I CMLIB
comes in three forms:
(a) an overall summary, 
(b) sublibrary summaries, and 
(c) detailed documentation on individual subprograms.
Sublibraries are subsets of the library which may be considered
as separate entities, such as the EISPACK or FISHPAK packages.
The overall summary contains a list of sublibraries and their purposes.
The command
.I cmdoc
is used to retrieve this documentation; it finds documentation on the
items selected and deposits it on standard output.
For example,

.I cmdoc \-s
gets the overall summary
.br
.I cmdoc \-s eispack
gets the EISPACK summary 
.br
.I cmdoc sgefa sgesl
gets details on SGEFA and SGESL.

Two additional commands are available to retrieve the Fortran source for 
.I CMLIB
subprograms.
.I cmsrc 
retrieves the Fortran source for each subprogram whose name is listed
on the command line.
.I cmsrcx
retrieves the same source as 
.I cmsrc
but also includes the source for all external references within the
library which are required to use the named subprograms.

Another way of discovering the appropriate CMLIB routine for a given problem
is to use the NIST Guide to Available Mathematical Software (GAMS).
The gams command will connect you to an interactive
problem-oriented catalog of more than 6500 software modules available
for use at NIST.
.SH FILES
.ta 2i
/usr/local/lib/libcm.a     library file
.br
/usr/local/apps/cmlib/src  source directory
.br
/usr/local/apps/cmlib/doc  doc directory
.br
/usr/local/bin/cmdoc       documentation retrieval script
.br
/usr/local/bin/cmsrc       source retrieval script
.br
/usr/local/bin/cmsrcx      source+externals retrieval
.SH SEE\ ALSO
f77(1), gams(1)
