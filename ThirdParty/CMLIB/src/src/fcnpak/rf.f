      REAL FUNCTION rf(x,y,z,ier)
c***BEGIN PROLOGUE  rf
c***DATE WRITTEN   790801   (yymmdd)
c***REVISION DATE  830622   (yymmdd)
c***CATEGORY NO.  b5
c***KEYWORDS  elliptic integral,incomplete,complete,
c             integral of the first kind.
c             duplication theorem,taylor series,
c***AUTHOR  carlson b.c. ames laboratory-doe,iowa state university
c           notis e.m.   ames laboratory-doe,iowa state university
c                        ames, iowa  50011
c           pexton r.l.  lawrence livermore national laboratory
c                        livermore, california  94550
c***PURPOSE  the routine calculates an approximation result to
c           rf(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -1/2
c                      (1/2)(t+x)    (t+y)    (t+z)    dt,
c           where x, y, and z are nonnegative and at most one of them
c           is zero.  if one of them  is zero, the integral  is complete
c***DESCRIPTION
c
c ......................................................................
c
c   1.     rf
c          evaluate an incomplete (or complete) elliptic integral
c          of the first kind
c          standard fortran function routine
c          real version
c          the routine calculates an approximation result to
c          rf(x,y,z) = integral from zero to infinity of
c
c                               -1/2     -1/2     -1/2
c                     (1/2)(t+x)    (t+y)    (t+z)    dt,
c
c          where x, y, and z are nonnegative and at most one of them
c          is zero.  if one of them  is zero, the integral  is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c
c   2.     calling sequence
c          rf( x, y, z, ier )
c
c          parameters on entry     values assigned by the calling subrou
c
c          x      - real,nonnegative variable
c
c          y      - real,nonnegative variable
c
c          z      - real,nonnegative variable
c
c
c
c          on return     values assigned by the rf routine
c
c          rf     - real approximation to the integral
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine. it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c          x, y, z are unaltered.
c
c
c   3.    error messages
c
c
c         value of ier assigned by the rf subroutine
c
c                  value assigned         error message printed
c                  ier = 1                amin1(x,y,z) .lt. 0.0e0
c                      = 2                amin1(x+y,x+z,y+z) .lt. lolim
c                      = 3                amax1(x,y,z) .gt. uplim
c
c
c
c   4.     control parameters
c
c                  values of lolim,uplim,and errtol are set by the
c                  subroutine.
c
c          lolim and uplim determine the valid range of x, y and z
c
c          lolim  - lower limit of valid arguments
c
c                   not less than 5 * (machine minimum).
c
c          uplim  - upper limit of valid arguments
c
c                   not greater than (machine maximum) / 5.
c
c
c                     acceptable values for:   lolim      uplim
c                     ibm 360/370 series   :   3.0e-78     1.0e+75
c                     cdc 6000/7000 series :   1.0e-292    1.0e+321
c                     univac 1100 series   :   1.0e-37     1.0e+37
c                     cray                 :   2.3e-2466   1.09e+2465
c                     vax 11 series        :   1.5e-38     3.0e+37
c
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the subroutine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c
c
c          errtol - relative error due to truncation is less than
c                   errtol ** 6 / (4 * (1-errtol)  .
c
c
c
c              the accuracy of the computed approximation to the integra
c              can be controlled by choosing the value of errtol.
c              truncation of a taylor series after terms of fifth order
c              introduces an error less than the amount shown in the
c              second column of the following table for each value of
c              errtol in the first column.  in addition to the truncatio
c              error there will be round-off error, but in practice the
c              total error from both sources is usually less than the
c              amount given in the table.
c
c
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    3.0e-19
c                           3.0e-3    2.0e-16
c                           1.0e-2    3.0e-13
c                           3.0e-2    2.0e-10
c                           1.0e-1    3.0e-7
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c
c***LONG DESCRIPTION
c
c              rf special comments
c
c
c
c          check by addition theorem: rf(x,x+z,x+w) + rf(y,y+z,y+w)
c          = rf(0,z,w), where x,y,z,w are positive and x * y = z * w.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral rf(x,y,z).
c
c
c          on output:
c
c
c          x, y, z are unaltered.
c
c
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c          expense of robustness.
c
c
c
c    ###################################################################
c    ###################################################################
c
c
c                  special functions via rf
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c
c                  legendre form of elliptic integral of 1st kind
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c                                            2         2   2
c                  F(phi,k) = sin(phi) rf(cos (phi),1-k sin (phi),1)
c
c
c                                 2
c                  K(k) = rf(0,1-k ,1)
c
c        In these formulae k is known as the modulus.  Elliptic
c        integrals are often tabulated as functions of k*k or arcsin(k).
c
c
c
c
c                  bulirsch form of elliptic integral of 1st kind
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c                                         2 2    2
c                  el1(x,kc) = x rf(1,1+kc x ,1+x )
c
c
c
c
c
c
c
c                  lemniscate constant a
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c                       1      4 -1/2
c                  a = int (1-s )    ds = rf(0,1,2) = rf(0,2,1)
c                       0
c
c
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c
c
c          subroutines or functions needed
c              - xerror
c              - r1mach
c              - fortran abs, amax1,amin1, sqrt
c
c***REFERENCES
c
c
c          carlson,b.c.and notis,e.m.
c          algorithms for incomplete elliptic integrals
c          acm transactions on mathematical software,vol.7,no.3,
c          sept,1981,pages 398-403
c
c
c          carlson,b.c.
c          computing elliptic integrals by duplication
c          numer.math.33,(1979),1-16
c
c
c          carlson,b.c.
c          elliptic integrals of the first kind
c          siam j.math.anal.8(1977),231-242
c
c
c
c***ROUTINES CALLED  xerror,r1mach
c
c
c***END PROLOGUE  rf
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c                  rf routine begins
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c
c
      integer ier
c
c
      real lolim, uplim, epslon, errtol
c
      real c1, c2, c3, e2, e3, lamda
c
      real mu, s, x, xn, xndev
c
c
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
c
c
c
c***FIRST EXECUTABLE STATEMENT  rf
c
c
      errtol=(4.0*r1mach(3))**(1.0/6.0)
c
c
      lolim = 5.0e0 * (r1mach(1) )
c
      uplim = (r1mach(2) ) / 5.0e0
c
      c1=1.0e0/24.0e0
c
      c2=3.0e0/44.0e0
c
      c3=1.0e0/14.0e0
c
c
c
c
c         call error handler if necessary.
c
c
    5 rf=0.0e0
      if( amin1(x,y,z).lt.0.0e0) then
      ier=1
      call xerror (36h*** rf   error amin1(x,y,z).lt.0.0e0,36,1,1 )
      return
      endif
      if (amin1(x+y,x+z,y+z).lt.lolim) then
      ier=2
      call xerror (42h*** rf   error amin1(x+y,x+z,y+z).lt.lolim,42,1,2)
      return
      endif
      if (amax1(x,y,z).gt.uplim) then
      ier=3
      call xerror (36h*** rf   error amax1(x,y,z).gt.uplim,36,3,1 )
      return
      endif
c
c
c
c
   20 ier = 0
      xn = x
      yn = y
      zn = z
c
c
c
   30 mu = (xn+yn+zn)/3.0e0
      xndev = 2.0e0 - (mu+xn)/mu
      yndev = 2.0e0 - (mu+yn)/mu
      zndev = 2.0e0 - (mu+zn)/mu
      epslon = amax1( abs(xndev), abs(yndev), abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      go to 30
c
c
c
c
   40 e2 = xndev*yndev - zndev*zndev
      e3 = xndev*yndev*zndev
      s = 1.0e0 + (c1*e2-0.10e0-c2*e3)*e2 + c3*e3
      rf = s/ sqrt(mu)
c
   50 return
c
c
c
c
      END
