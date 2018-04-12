      subroutine ryork (n, x, y, wx, wy, bstart, mxiter, delmax,
     1                  ires, maxd, work, output, res, ier)
c
c-----------------------------------------------------------------------
c
c     Package:  SLRPACK
c
c     Version:  October, 1985
c
c-----------------------------------------------------------------------
c
c     PURPOSE
c     -------
c
c     This user-callable subroutine estimates simple linear regression
c     coefficients when both variables are subject to errors which are
c     not necessarily homogeneous in variance. The algorithm for the
c     estimation is due to York (1966) (see reference below).
c
c
c     DESCRIPTION
c     -----------
c
c     1. The input data are observations (x(i),y(i)), i = 1,...,n, and
c        weights wx(i) and wy(i) associated with observations x(i) and
c        y(i) respectively, for i = 1,...,n.
c
c     2. The function
c
c        S = sum over i (wx(i)*(x(i)-X(i))**2 + wy(i)*(y(i)-Y(i))**2)
c
c        is minimized with respect to the fitted values X(i) and Y(i),
c        i = 1,...,n. The set of points (X(i),Y(i)) i = 1,...,n, are
c        constrained to be collinear.
c
c        This subroutine straightforwardly implements the calculations
c        described in York (1966) (see references below).
c
c     3. The slope estimate is obtained iteratively, with an initial
c        estimate supplied by the user. One suggested initial estimate
c        is the geometric mean of the OLS-y and OLS-x slopes obtained
c        by executing SLRPACK subroutine RGM (see documentation of that
c        subroutine).
c
c     4. At each iteration a cubic equation is solved and the "best" of
c        the up to three possible real-valued solutions is selected as
c        the slope estimate at that iteration.
c
c        There are between one and three solutions of the cubic equation
c        at a given iteration. If there is one real solution (and two co
c        plex solutions), the iterations proceed. If there is more than
c        one real solution, the solution closest to the previous solutio
c        is found and the iterations proceed.
c
c     5. The iterations terminate when 1) the angle between the lines es
c        imated in two successive iterations is less than user-set DELMA
c        or 2) the maximum number of iterations (user-set MXITER) has be
c        reached, whichever occurs first.
c
c     6. The equations for the weighted averages of the x- and y-
c        observations are given by equations (19) in the York reference,
c        and the equations for the two standard deviations are given
c        on pages 1084 and 1085 of the same reference.
c
c     7. Another algorithm for minimizing S, described by Williamson
c        (1968), is implemented in the SLRPACK subroutine RWILL.
c
c
c     INPUT PARAMETERS
c     ----------------
c
c     N       Integer scalar (unchanged on output)
c             Number of observations.
c
c             N must be at least 3.
c
c     X       Real vector dimensioned at least N (unchanged on output)
c             The X-observations.
c
c     Y       Real vector dimensioned at least N (unchanged on output)
c             The Y-observations.
c
c     WX      Real vector dimensioned at least N (unchanged on output)
c             Weights associated with x-observations (commonly the
c             multiplicative inverses of the variances associated with
c             the x-observations).
c
c             WX(i) > 0   for i = 1,...,N
c
c     WY      Real vector dimensioned at least N (unchanged on output)
c             Weights associated with y-observations (commonly the
c             multiplicative inverses of the variances associated with
c             the y-observations).
c
c             WY(i) > 0   for i = 1,...,N
c
c     BSTART  Real scalar (unchanged on output)
c             Initial slope estimate.
c
c     MXITER  Integer scalar (unchanged on output)
c             Maximum number of iterations allowed.
c
c     DELMAX  Real scalar (unchanged on output)
c             Angular distance in radians such that execution terminates
c             if the angular distance between lines estimated in two
c             successive iterations is less than DELMAX.
c
c     IRES    Integer scalar (unchanged on output)
c
c             IRES = 1 if residuals are to be calculated
c             IRES = 0 otherwise
c
c     MAXD    Integer scalar (unchanged on output)
c             Row dimension of RES in the calling program.
c
c             MAXD must be at least N if IRES = 1, otherwise MAXD
c             may be 1.
c
c     WORK    Real vector dimensioned at least 10 + (9*N)
c             Work vector.
c
c     OUTPUT PARAMETERS
c     -----------------
c
c     OUTPUT  Real vector dimensioned at least 7
c
c             OUTPUT(1) = slope estimate at final iteration
c
c             OUTPUT(2) = intercept estimate at final iteration
c
c             OUTPUT(3) = standard deviation of slope estimate
c
c             OUTPUT(4) = standard deviation of intercept estimate
c
c             OUTPUT(5) = weighted average of x observations at final
c                         iteration
c
c             OUTPUT(6) = weighted average of y observations at final
c                         iteration
c
c             OUTPUT(7) = number of iterations
c
c     RES     Real matrix dimensioned at least MAXD by 2
c             The i-th x and y residuals are in RES(i,1) and  RES(i,2)
c             respectively.
c
c     IER     Integer scalar
c
c             IER =  0  if there are no execution time errors or
c                          warnings
c
c             IER =  1  York technique regression line slope estimates
c                          cannot be calculated because the data set is
c                          too small
c
c                       Cannot compute OUTPUT(I) for I=1,...,7
c                          and cannot compute the residuals
c
c             IER =  2  York technique regression line slope estimates
c                          cannot be calculated because all x-observatio
c                          are equal
c
c                       Cannot compute OUTPUT(I) for I=1,2,3,4,7
c                          and cannot compute the residuals
c
c             IER =  3  York technique regression line slope estimates
c                          cannot be calculated because all of the
c                          weights are not positive
c
c                       Cannot compute OUTPUT(I) for I=1,...,7
c                          and cannot compute the residuals
c
c             IER =  4  York technique slope estimates have not converge
c                          in allowed number of iterations
c
c                       Have computed OUTPUT(I) for I=1,...,7 and residu
c                          for the final iteration
c
c     EXAMPLE
c     -------
c
c     Input:
c
c              N = 10
c           IRES =  1
c         MXITER = 20
c         BSTART =  -.54600
c         DELMAX =   .00010
c
c             I        X(I)      WX(I)       Y(I)      WY(I)
c             1        .000   1000.000      5.900      1.000
c             2        .900   1000.000      5.400      1.800
c             3       1.800    500.000      4.400      4.000
c             4       2.600    800.000      4.600      8.000
c             5       3.300    200.000      3.500     20.000
c             6       4.400     80.000      3.700     20.000
c             7       5.200     60.000      2.800     70.000
c             8       6.100     20.000      2.800     70.000
c             9       6.500      1.800      2.400    100.000
c            10       7.400      1.000      1.500    500.000
c
c     Output:
c
c                   IER = 0
c
c                   OUTPUT(1) =  -.481
c                   OUTPUT(2) =  5.480
c                   OUTPUT(3) =   .071
c                   OUTPUT(4) =   .362
c                   OUTPUT(5) =  4.911
c                   OUTPUT(6) =  3.120
c                   OUTPUT(7) =  5.
c
c              Weighted Residuals
c           I         RES(I,1)    RES(I,2)
c           1           .000       -.420
c           2           .000       -.352
c           3           .001        .215
c           4          -.002       -.369
c           5           .019        .385
c           6          -.038       -.316
c           7           .080        .143
c           8          -.234       -.139
c           9          -.084       -.003
c          10           .875        .004
c
c     Note:  Given initial slope bstart = -0.546, after one "iteration"
c            of the York technique, the estimated slope is -0.477, in
c            agreement with York (1966).
c
c     PRECISION
c     ---------
c
c     All calculations are done in single precision.
c
c     LANGUAGE
c     --------
c
c     The routine is written in standard Fortran 77.
c
c     OTHER SUBROUTINES
c     -----------------
c
c     SLRPACK  subroutines  WYORK, RYBAR, RYARNG
c
c     PORT     subroutines  R1MACH
c
c     REFERENCES
c     ----------
c
c     Madansky, A. (1959). The fitting of straight lines when both
c        variables are subject to error. JASA, 54, 173-205.
c
c     Riggs, D. S., Guarnieri, J. A., and Addelman, S. (1978).  Fitting
c        straight lines when both variables are subject to error.  Life
c        Sciences, 22, 1305-1360.
c
c     York, D. (1966). Least-square fitting of a straight line. Canadian
c        Journal of Physics, 44, 1079-1086.
c
c     Williamson, J. H. (1968). Least-squares fitting of a straight line
c        Canadian Journal of Physics, 46, 1845-1847.
c
c     CONTRIBUTORS
c     ------------
c
c     Sally E. Howe
c     Katherine Rensenbrink
c     Amy DelGiorno
c     Gregory Rhoads
c     Scientific Computing Division
c     Center for Applied Mathematics
c     National Bureau of Standards
c     Gaithersburg, MD  20899
c
c     NBS CONTACT
c     -----------
c
c     Sally E. Howe
c     Center for Applied Mathematics
c
c-----------------------------------------------------------------------
c
      real x(*), y(*), wx(*), wy(*), res(maxd,*), output(*), work(*)
      eps = r1mach(4)
c
c------------------------------------------
c     Initialize and check for input errors
c------------------------------------------
c
      do 10 i = 1 , 7
 10      output(i) = 0.0
c
      do 20 i = 1 , maxd
         do 20 j = 1 , 2
 20         res(i,j) = 0.0
c
      if (n .le. 2) then
         ier = 1
         return
      endif
      do 30 i = 1 , n
         if (wx(i) .le. eps .or. wy(i) .le. eps) then
            ier = 3
            return
         endif
 30   continue
c
c------------------------
c     Allocate work space
c------------------------
c
      iang   = 1
      iu     = iang   + 4
      iu2    = iu     + n
      iuv    = iu2    + n
      iv     = iuv    + n
      iv2    = iv     + n
      iw     = iv2    + n
      iw3    = iw     + n
      iwwx   = iw3    + n
      iwwy   = iwwx   + n
      ib     = iwwy   + n
      iphi   = ib     + 3
c
      work(ib+2) = 0
      a          = 0
      sda        = 0
      sdb        = 0
      xbar       = 0
      ybar       = 0
c
c--------------------------------
c     Call the workhorse of ryork
c--------------------------------
c
      call wyork (n, x, y, wx, wy, bstart, mxiter, delmax, ires,
     1            res(1,1), res(1,2), ier, eps, work(iang), work(iu),
     2            work(iu2), work(iuv), work(iv), work(iv2), work(iw),
     3            work(iw3), work(iwwx), work(iwwy), a, work(ib), xbar,
     4            ybar, sdb, sda, iter, work(iphi))
c
c-------------------------------------
c     Place results into output vector
c-------------------------------------
c
      output(1) =  work(ib+2)
      output(2) =  a
      output(3) =  sdb
      output(4) =  sda
      output(5) =  xbar
      output(6) =  ybar
      output(7) =  iter
      return
      end
