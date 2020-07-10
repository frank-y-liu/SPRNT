.TH CMLIB 3F "June 2000"
.SH NAME
CMLIB \- NIST Core Math Library (Version 3, April 1988)
.SH SYNOPSIS
.br
/usr/local/lib/libcm.a 
.br
/usr/local/lib/libcm_autodbl.a 
.br
/usr/local/lib/libcm_f90.a 
.br
/usr/local/lib/libcm_f90_autodbl.a 
.LP
\fBcmdoc\fP [ \fB\-s\fP [ sublibrary ] ] .\|.\|. [ subprogram ] .\|.\|.
.br
\fBcmsrc\fP subprogram .\|.\|.
.br
\fBcmsrcx\fP subprogram .\|.\|.
.SH DESCRIPTION
.LP
.I CMLIB
is the NIST Core Math Library, a collection of Fortran-callable
subprograms for solving a wide range of mathematical problems. 
Four implementations are maintained on the RS/6000: 
a standard Fortran 77 and Fortran 90 implementation, and two compiled 
using automatic precision doubling option (-qautodbl=dblpad).
All are compiled with the -O optimization flag.
.LP 
One selects a version of the library by including appropriate options 
on the compilation command.
.sp 1
.in +4
.ti -4
(1) To link to the Fortran 77 standard version use
 
  xlf -o prog prog.f -lcm
 
.ti -4
(2) To link to the Fortran 77 version compiled with automatic precision doubling, use
 
  xlf -qautodbl=dblpad -o prog prog.f -lcm_autodbl
 
.ti -4
(3) To link to the Fortran 90 standard version use
 
  xlf90 -o prog prog.f -lcm_f90
 
.ti -4
(4) To link to the Fortran 90 version compiled with automatic precision doubling, use
 
  xlf90 -qautodbl=dblpad -o prog prog.f -lcm_f90_autodbl
 
.LP
It is very important to use the version of the library that is 
compatible with the xlf options that you are using.  An incorrect
combination will lead to incorrect results, segmentation faults or bus errors.  
.LP
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
.LP
.I cmdoc \-s
gets the overall summary
.br
.I cmdoc \-s eispack
gets the EISPACK summary 
.br
.I cmdoc sgefa sgesl
gets details on SGEFA and SGESL.
.LP
Two additional commands are available to retrieve the Fortran source for 
.I CMLIB
subprograms.
.I cmsrc 
retrieves the Fortran source for each subprogram whose name is listed
on the command line.
.I cmsrcx
retrieves the same source as
.I cmsrc
but also includes the source for all external references within
the library which are required to use the named subprograms.
.SH FILES
.ta 3i
In /usr/local/lib ...
.br
libcm.a		library archive
.br
libcm_autodbl.a	library archive (autodbl)
.br
libcm_f90.a			library archive
.br
libcm_f90_autodbl.a		library archive (autodbl)
.sp
/usr/local/apps/cmlib/doc	documentation directory
.br
/usr/local/apps/cmlib/src	source directory
.br
/usr/local/bin/cmdoc	documentation retrieval script
.br
/usr/local/bin/cmsrc	source retrieval script
.br
/usr/local/bin/cmsrcx	source+externals retrieval
.SH SEE\ ALSO
xlf(1)
.SH BUGS
.LP
Odrpack may be unsuitable for use on workstations in single precision.
