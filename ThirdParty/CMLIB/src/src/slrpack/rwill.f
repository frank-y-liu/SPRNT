      subroutine rwill(n, x, y, u, v, bstart, mxiter, delmax,
     1                 output, xres, yres, work, ier)
c
c---------------------------------------------------------------------
c
c   Package:  SLRPACK
c
c   Version:  October, 1985
c
c----------------------------------------------------------------------
c
c     PURPOSE
c     -------
c
c     This user-callable subroutine estimates simple linear regression
c     coefficients when both variables are subject to errors which are
c     not necessarily homogeneous in variance. The algorithm for the
c     estimation is due to Williamson (1968) (see reference below).
c
c
c     DESCRIPTION
c     -----------
c
c     1. The input data are observations (x(i),y(i)), i = 1,...,n, and
c        variances u(i) and v(i) associated with observations x(i) and
c        y(i) respectively, for i = 1,...,n.
c
c     2. The function
c
c        S = sum over i ((x(i)-X(i))**2/u(i) + (y(i)-Y(i))**2/v(i))
c
c        is minimized with respect to the fitted values X(i) and Y(i),
c        i = 1,...,n. The set of points (X(i),Y(i)) i = 1,...,n, are
c        constrained to be collinear.
c
c        This subroutine straightforwardly implements the calculations
c        described in Williamson (1968) (see references).
c
c     3. The slope estimate is obtained iteratively, with an initial
c        estimate supplied by the user. One suggested initial estimate
c        is the geometric mean of the OLS-y and OLS-x slopes obtained
c        by executing SLRPACK subroutine RGM (see documentation of that
c        subroutine).
c
c     4. At each iteration a linear equation is solved and the solution
c        is selected as the slope estimate at that iteration.
c
c     5. The iterations terminate when 1) the angle between the lines es
c        imated in two successive iterations is less than user-set DELMA
c        or 2) the maximum number of iterations (user-set MXITER) has be
c        reached, whichever occurs first.
c
c     6. The equations for the standard deviations and the weighted
c        averages are given in Williamson (1968).
c
c     7. Another algorithm for minimizing S, described by York (1966),
c        is implemented in the SLRPACK subroutine RYORK.
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
c     U       Real vector dimensioned at least N (unchanged on output)
c             U(i) is the variance associated with X(i).
c
c             U(i) > 0.0  for i = 1,...,N.
c
c     V       Real vector dimensioned at least N (unchanged on output)
c             V(i) is the variance associated with Y(i).
c
c             V(i) > 0.0  for i = 1,...,N.
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
c     WORK    Real vector dimensioned at least 5 * N
c             Work vector.
c
c     OUTPUT PARAMETERS
c     -----------------
c
c     OUTPUT  Real vector dimensioned at least 11
c
c             OUTPUT(1) = Slope estimate at final iteration.
c
c             OUTPUT(2) = Intercept estimate at final iteration.
c
c             OUTPUT(3) = Standard deviation of slope estimate.
c
c             OUTPUT(4) = Standard deviation of intercept estimate.
c
c             OUTPUT(5) = Standard deviation of the slope estimate
c                         multiplied by the chi-square test statistic
c                         weighted sum of squares divided by degrees
c                         of freedom.
c
c             OUTPUT(6) = Standard deviation of the intercept estimate
c                         multiplied by the chi-square test statistic
c                         weighted sum of squares divided by degrees
c                         of freedom.
c
c             OUTPUT(7) = Weighted average of x observations at final
c                         iteration.
c
c             OUTPUT(8) = Weighted average of y observations at final
c                         iteration.
c
c             OUTPUT(9) = Root mean square error in the x-direction.
c
c             OUTPUT(10)= Root mean square error in the y-direction.
c
c             OUTPUT(11)= Number of iterations.
c
c     XRES    Real vector dimensioned at least N
c             The x residuals.
c
c     YRES    Real vector dimensioned at least N
c             The y residuals.
c
c     IER     Integer scalar
c
c             IER =  0  if there are no execution time errors or
c                          warnings
c
c             IER =  1  Williamson technique regression line slope
c                          estimates cannot be calculated because the
c                          data set is too small
c
c                       Cannot compute OUTPUT(I) for I=1,...,11
c                          and residuals
c
c             IER =  2  Williamson technique regression line slope
c                          estimates cannot be calculated because all
c                          X-observations are equal
c
c                       Cannot compute OUTPUT(I) for I=1,...,11
c                          and residuals
c
c             IER =  3  Williamson technique regression line slope
c                          estimates cannot be calculated because all
c                          variances are not greater than 0.0
c
c                       Cannot compute OUTPUT(I) for I=1,...,11
c                          and residuals
c
c             IER =  4  Williamson technique regression line slope
c                          estimates cannot be determined because the su
c                          of the weights is zero (equation for W(i) is
c                          on page 1846 of Williamson).
c
c                       Cannot compute OUTPUT(I) for I=1,...,11
c                          and residuals
c
c             IER =  5  Uncertainties in slope and intercept cannot be
c                          calculated since the final slope is 0.
c
c                       Cannot compute OUTPUT(I) for I=3,4,5,6,9,10
c                          and residuals
c
c             IER =  6  Uncertainties in slope and intercept cannot be
c                          determined because certain calculations invol
c                          division by zero
c
c                       Cannot compute OUTPUT(I) for I=3,4,5,6,9,10
c                          and residuals
c
c             IER =  7  Williamson technique slope estimates have not
c                          converged in allowed number of iterations
c
c                       Have computed OUTPUT(I) for I=1,...,11 only for
c                          the last iteration.
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
c             I        X(I)      U(I)        Y(I)      V(I)
c             1       0.000     0.00100     5.900     1.00000
c             2       0.900     0.00100     5.400     0.55556
c             3       1.800     0.00200     4.400     0.25000
c             4       2.600     0.00125     4.600     0.12500
c             5       3.300     0.00500     3.500     0.05000
c             6       4.400     0.01250     3.700     0.05000
c             7       5.200     0.01667     2.800     0.01429
c             8       6.100     0.05000     2.800     0.01429
c             9       6.500     0.55556     2.400     0.01000
c            10       7.400     1.00000     1.500     0.00200
c
c     Output:
c
c                   IER = 0
c
c                   OUTPUT(1) =  -.481
c                   OUTPUT(2) =  5.480
c                   OUTPUT(3) =   .071
c                   OUTPUT(4) =   .362
c                   OUTPUT(5) =   .0576
c                   OUTPUT(6) =   .2919
c                   OUTPUT(7) =  4.911
c                   OUTPUT(8) =  3.120
c                   OUTPUT(9) =   .2888
c                   OUTPUT(10)=   .2776
c                   OUTPUT(11)=  4.
c
c              Weighted Residuals
c
c           I          XRES(I)    YRES(I)
c           1           .000       -.420
c           2           .000       -.352
c           3           .001        .214
c           4          -.002       -.369
c           5           .019        .385
c           6          -.038       -.316
c           7           .080        .143
c           8          -.234       -.139
c           9          -.085       -.003
c          10           .874        .004
c
c     PRECISION
c     ---------
c
c     All calculations are done in single precision.
c
c     OTHER SUBROUTINES
c     -----------------
c
c     PORT  subroutine  R1MACH
c
c     LANGUAGE
c     --------
c
c     The routine is coded in standard Fortran 77.
c
c     REFERENCE
c     ---------
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
c
c     AUTHORS
c     -------
c
c     R. Lindstrom
c     Inorganic Analytical Reasearch Division
c     Center for Analytical Chemistry
c
c     Sally E. Howe
c     Scientific Computing Division
c     Center for Applied Mathematics
c
c     National Bureau of Standards
c     Gaithersburg, MD  20899
c
c     NBS CONTACT
c     -----------
c
c     Sally E. Howe
c     Scientific Computing Division
c
c-----------------------------------------------------------------------
c
c
      dimension x(*), y(*), u(*), v(*), xres(*), yres(*), output(*),
     1          work(*)
c
      eps = r1mach(4)
c
c------------------------------
c   initialize and check output
c------------------------------
c
      do 10 i = 1 , 11
 10      output(i) = 0.0
c
      if (n .le. 2) then
         ier = 1
         return
      endif
c
      do 20 i = 1 , n
         xres(i) = 0.0
 20      yres(i) = 0.0
c
      do 30 i = 1 , n
         if (u(i) .le. eps .or. v(i) .le. eps) then
            ier = 3
            return
         endif
 30   continue
c
c
c----------------------
c   allocate work space
c----------------------
c
      idelx  = 1
      idely  = idelx  + n
      idelz  = idely  + n
      iw     = idelz  + n
      iz     = iw     + n
c
c------------------------------
c   call the workhorse of rwill
c------------------------------
c
      call wwill(n, x, y, u, v, bstart, mxiter, delmax, xres,
     1           yres, ier, eps, iter, a, b, adjsda, adjsdb,
     2           sda, sdb, xbar, ybar, rmsx, rmsy, work(idelx),
     3           work(idely), work(idelz), work(iw), work(iz))
c
c-----------------------------------
c   place results into output vector
c-----------------------------------
c
      output(1) = b
      output(2) = a
      output(3) = sdb
      output(4) = sda
      output(5) = adjsdb
      output(6) = adjsda
      output(7) = xbar
      output(8) = ybar
      output(9) = rmsx
      output(10)= rmsy
      output(11)= iter
c
      return
      end
