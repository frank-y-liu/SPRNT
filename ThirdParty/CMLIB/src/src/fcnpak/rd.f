      REAL FUNCTION rd(x,y,z,ier)
c***BEGIN PROLOGUE  rd
c***DATE WRITTEN   790801   (yymmdd)
c***REVISION DATE  830622   (yymmdd)
c***CATEGORY NO.  b5
c***KEYWORDS  elliptic integral,incomplete,complete,
c             integral of the second kind.
c             duplication theorem,taylor series,
c***AUTHOR  carlson b.c. ames laboratory-doe,iowa state university
c           notis e.m.   ames laboratory-doe,iowa state university
c                        ames, iowa  50011
c           pexton r.l.  lawrence livermore national laboratory
c                        livermore, california  94550
c***PURPOSE  the routine calculates an approximation result to
c            rd(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -3/2
c                      (3/2)(t+x)    (t+y)    (t+z)    dt,
c            where x and y are nonnegative, x + y is positive, and z is
c            positive.  if x or y is zero, the integral is complete.
c***DESCRIPTION
c
c ......................................................................
c
c   1.     rd
c          evaluate an incomplete (or complete) elliptic integral
c          of the second kind
c          standard fortran function routine
c          real version
c          the routine calculates an approximation result to
c          rd(x,y,z) = integral from zero to infinity of
c                              -1/2     -1/2     -3/2
c                    (3/2)(t+x)    (t+y)    (t+z)    dt,
c          where x and y are nonnegative, x + y is positive, and z is
c          positive.  if x or y is zero, the integral is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c   2.     calling sequence
c          rd( x, y, z, ier )
c
c          parameters on entry     values assigned by the calling subrou
c
c          x      - real,nonnegative variable
c
c          y      - real,nonnegative variable
c
c                   x + y is positive
c
c          z      - real,positive variable
c
c
c
c          on return     values assigned by the rd subroutine
c
c          rd     - real approximation to the integral
c
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine. it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c
c          x, y, z are unaltered.
c
c   3.    error messages
c
c
c         value of ier assigned by the rd subroutine
c
c                  value assigned         error message printed
c                  ier = 1                amin1(x,y) .lt. 0.0e0
c                      = 2                amin1(x + y, z ) .lt. lolim
c                      = 3                amax1(x,y,z) .gt. uplim
c
c
c   4.     control parameters
c
c                  values of lolim,uplim,and errtol are set by the
c                  subroutine.
c
c          lolim and uplim determine the valid range of x, y, and z
c
c          lolim  - lower limit of valid arguments
c
c                    not less  than 2 / (machine maximum) ** (2/3).
c
c          uplim  - upper limit of valid arguments
c
c                    not greater than (0.1e0 * errtol / machine
c                    minimum) ** (2/3), where errtol is described below.
c                    in the following table it is assumed that errtol wi
c                    never be chosen smaller than 1.0e-5.
c
c
c                    acceptable values for:   lolim      uplim
c                    ibm 360/370 series   :   6.0e-51     1.0e+48
c                    cdc6000/7000 series  :   5.0e-215    2.0e+191
c                    univac1100 series    :   1.0e-25     2.0e+21
c                    cray                 :   3.0e-1644   1.69e+1640
c                    vax 11 series        :   1.0e-25     4.5e+21
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the subroutine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c          errtol    relative error due to truncation is less than
c                    3 * errtol ** 6 / (1-errtol) ** 3/2.
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
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    4.0e-18
c                           3.0e-3    3.0e-15
c                           1.0e-2    4.0e-12
c                           3.0e-2    3.0e-9
c                           1.0e-1    4.0e-6
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c
c
c
c***LONG DESCRIPTION
c
c              rd special comments
c
c
c
c          check: rd(x,y,z) + rd(y,z,x) + rd(z,x,y)
c          = 3 /  sqrt(x * y * z), where x, y, and z are positive.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral rd(x,y,z).
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
c                  special functions via rd and rf
c
c
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                  legendre form of elliptic integral of 2nd kind
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c                                            2         2   2
c                  E(phi,k) = sin(phi) rf(cos (phi),1-k sin (phi),1) -
c
c                     2      3            2         2   2
c                  -(k/3) sin (phi) rd(cos (phi),1-k sin (phi),1)
c
c
c                                 2        2           2
c                  E(k) = rf(0,1-k ,1) - (k/3) rd(0,1-k ,1)
c
c      In these formulae k is known as the modulus.  Elliptic
c      integrals are often tabulated as functions of k*k or arcsin(k).
c
c
c
c                  bulirsch form of elliptic integral of 2nd kind
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                                              2 2    2
c                  el2(x,kc,a,b) = ax rf(1,1+kc x ,1+x ) +
c
c                                              3         2 2    2
c                                 +(1/3)(b-a) x rd(1,1+kc x ,1+x )
c
c
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c                  legendre form of alternative elliptic integral
c                  of 2nd kind
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c                            q     2       2   2  -1/2
c                  D(q,k) = int sin p  (1-k sin p)     dp
c                            0
c
c
c
c                                     3          2     2   2
c                  D(q,k) = (1/3) (sin q)  rd(cos q,1-k sin q,1)
c
c      In these formulae k is known as the modulus.  Elliptic
c      integrals are often tabulated as functions of k*k or arcsin(k).
c
c
c
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c                  lemniscate constant  b
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c                       1    2    4 -1/2
c                  b = int  s (1-s )    ds
c                       0
c
c
c                  b = (1/3) rd (0,2,1)
c
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c                        '
c                  heuman s lambda function
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c
c                  (pi/2) lambda0(a,b) =
c
c                                         2                2
c                       = sin(b) (rf(0,cos (a),1)-(1/3) sin (a) *
c
c                                 2              2         2       2
c                        *rd(0,cos (a),1)) rf(cos (b),1-cos (a) sin (b),
c
c                                 2       3            2
c                       -(1/3) cos (a) sin (b) rf(0,cos (a),1) *
c
c                               2         2       2
c                        *rd(cos (b),1-cos (a) sin (b),1)
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c                  jacobi zeta function
c
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                             2                2       2   2
c                  Z(b,k) = (k/3) sin(b) rf(cos (b),1-k sin (b),1)
c
c
c                                      2            2
c                             *rd(0,1-k ,1)/rf(0,1-k ,1)
c
c                               2       3          2       2   2
c                            -(k /3) sin (b) rd(cos (b),1-k sin (b),1)
c
c      In these formulae k is known as the modulus.  Elliptic
c      integrals are often tabulated as functions of k*k or arcsin(k).
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
c***END PROLOGUE  rd
c
c
c
c
c
c                  rd routine begins
c                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c
c
      integer ier
c
      real lolim, uplim, epslon, errtol
c
      real c1, c2, c3, c4, ea, eb, ec, ed, ef, lamda
c
      real mu, power4, sigma, s1, s2, x, xn, xndev
c
c
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
c
c
c
c***FIRST EXECUTABLE STATEMENT  rd
c
c
      errtol=(r1mach(3)/3.0e0)**(1.0e0/6.0e0)
c
c
      lolim = 2.0e0/(r1mach(2))**(2.0e0/3.0e0)
c
      uplim=r1mach(1)**(1.0e0/3.0e0)
      uplim=(0.10e0*errtol)**(1.0e0/3.0e0)/uplim
      uplim=uplim**2.0e0
c
c
      c1 = 3.0e0/14.0e0
      c2 = 1.0e0/6.0e0
      c3 = 9.0e0/22.0e0
      c4 = 3.0e0/26.0e0
c
c
c
c
c         call error handler if necessary.
c
c
    5 rd=0.0e0
      if( amin1(x,y).lt.0.0e0) then
      ier=1
      call xerror (32h*** rd   error amin1(x,y).lt.0.0,32,1,1 )
      return
      endif
      if (amin1(x+y,z).lt.lolim) then
      ier=2
      call xerror (36h*** rd   error amin1(x+y,z).lt.lolim,36,2,1 )
      return
      endif
      if (amax1(x,y,z).gt.uplim) then
      ier=3
      call xerror (36h*** rd   error amax1(x,y,z).gt.uplim,36,3,1 )
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
      sigma = 0.0e0
      power4 = 1.0e0
c
c
c
   30 mu = (xn+yn+3.0e0*zn)*0.20e0
      xndev = (mu-xn)/mu
      yndev = (mu-yn)/mu
      zndev = (mu-zn)/mu
      epslon = amax1( abs(xndev), abs(yndev), abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      sigma = sigma + power4/(znroot*(zn+lamda))
      power4 = power4*0.250e0
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      go to 30
c
c
c
c
   40 ea = xndev*yndev
      eb = zndev*zndev
      ec = ea - eb
      ed = ea - 6.0e0*eb
      ef = ed + ec + ec
      s1 = ed*(-c1+0.250e0*c3*ed-1.50e0*c4*zndev*ef)
      s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
      rd = 3.0e0*sigma + power4*(1.0e0+s1+s2)/(mu* sqrt(mu))
c
   50 return
c
c
c
c
      END
