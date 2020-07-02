      subroutine rgm (n, x, y, output, ier)
c
c-----------------------------------------------------------------------
c
c   Package:  SLRPACK
c
c   Version:  October, 1985
c
c-----------------------------------------------------------------------
c
c   PURPOSE
c   -------
c
c   This subroutine computes estimates of simple linear regression
c   parameters for geometric mean regression.
c
c   DESCRIPTION
c   -----------
c
c   1.  The input data are observations (x(i), y(i)), i = 1,...,n.
c
c   2.  This technique is appropriate for fitting a straight line
c
c                          Y  =  A  +  B * X
c
c       to observations  (x(i),y(i)),i=1,...,n  where
c
c           x(i) = X(i) + e(i)   and   y(i) = Y(i) + d(i),
c
c           X(i) and Y(i) are unknown population means, and
c
c           e(i) and d(i) are random error terms distributed normally
c           with zero mean and variances s**2 = s(e)**2 = s(d)**2.
c
c   3.  The name geometric mean alludes to the fact that if B(OLS-y)
c       denotes the estimated slope from conventional least squares
c       regression (with the y-observations but not the x-observations
c       subject to error) and if B(OLS-x) denotes the converse, then
c       the geometric mean slope estimate is the geometric mean of
c       B(OLS-y) and B(OLS-x).  This technique for estimating the
c       regression parameters is preferable to OLS-y and to OLS-x
c       for the above type of data because the OLS-y technique is
c       known to give under-estimates of the magnitude of the true
c       slope and over-estimates of the magnitude of the true intercept
c       while OLS-x results have the opposite characteristics (Riggs
c       et al (1978)).
c
c   4.  This subroutine straightforwardly implements the calculations
c       described in the first reference.  The regression line is
c       given by equation (2), the standard deviations of the
c       estimated slope is given by equation (11), and the standard
c       deviation of the estimated intercept is given by equation (14)
c       all in the first reference.  The correlation coefficient is
c       the Pearson product-moment correlation coefficient.
c
c   INPUT PARAMETERS
c   ----------------
c
c   N      Integer scalar (unchanged on output)
c          Number of observations.
c
c          N must be greater than 2.
c
c   X      Real vector dimensioned at least N (unchanged on output)
c          X-observations.
c
c   Y      Real vector dimensioned at least N (unchanged on output)
c          Y-observations.
c
c   OUTPUT PARAMETERS
c   -----------------
c
c   OUTPUT Real vector dimensioned at least 7 -- output
c
c          OUTPUT(1) = Slope of geometric mean regression line
c
c          OUTPUT(2) = y-intercept of geometric mean regression line
c
c          OUTPUT(3) = Standard deviation of the slope
c
c          OUTPUT(4) = Standard deviation of the intercept
c
c          OUTPUT(5) = Average of the x-observations
c
c          OUTPUT(6) = Average of the y-observations
c
c          OUTPUT(7) = Pearson product-moment correlation coefficient
c
c   IER    Integer scalar -- output
c          Execution error indicator.
c
c          IER = 0   No errors
c
c          IER = 1   Cannot compute parameters by geometric mean techniq
c                    because data set is too small
c
c                    Cannot compute OUTPUT(I) for I = 1,...,7
c
c          IER = 2   Geometric mean results cannot be computed because
c                    all x-values are equal
c
c                    Cannot compute OUTPUT(I) for I = 1,2,3,4,7
c
c          IER = 3   Cannot compute standard deviations of slope and int
c                    estimates by geometric mean technique or correlatio
c                    coefficient because all y-values are equal
c
c                    Cannot comput OUTPUT(I) for I = 3, 4, 7
c
c   EXAMPLE
c   -------
c
c   INPUT:
c
c        N = 10
c
c        I          X(I)          Y(I)
c        1          0.0           5.9
c        2          0.9           5.4
c        3          1.8           4.4
c        4          2.6           4.6
c        5          3.3           3.5
c        6          4.4           3.7
c        7          5.2           2.8
c        8          6.1           2.8
c        9          6.5           2.4
c       10          7.4           1.5
c
c   OUTPUT:
c
c         IER = 0
c
c         OUTPUT(1) = -0.5526
c         OUTPUT(2) =  5.8108
c         OUTPUT(3) =  0.0377
c         OUTPUT(4) =  0.2465
c         OUTPUT(5) =  3.8200
c         OUTPUT(6) =  3.7000
c         OUTPUT(7) = -0.9765
c
c   See documentation of SLRPACK subroutines RYORK and RWILL for
c      continuations of this example.
c
c   PRECISION
c   ---------
c
c   All calculations are done in single precision.
c
c   LANGUAGE
c   --------
c
c   The routine is coded in standard Fortran 77.
c
c   OTHER SUBROUTINES USED
c   ----------------------
c
c   PORT  subroutine  R1MACH
c
c   REFERENCES
c   ----------
c
c   Kermack, K. A. and Haldane, J. B. S. (1950).  Organic Correlation
c      and Allometry.  Biometrika, 37, 30-41.
c
c   Pearson, K. (1901).  On lines and planes of closest fit to systems
c      of points in space.  Phil. Mag. (6), 2, 559.
c
c   Riggs, D. S., Guarnieri, J. A., and Addelman, S.  (1978).  Fitting
c      straight lines when both variables are subject to error.  Life
c      Sciences, 22, 1305-1360.
c
c   NBS CONTACT
c   -----------
c
c   Sally E. Howe
c   Scientific Computing Division
c
c   CONTRIBUTORS
c   ------------
c
c   Sally E. Howe
c   Kathryn Rensenbrink
c   Gregory S. Rhoads
c   Scientific Computing Division
c   Center for Applied Mathematics
c   National Bureau of Standards
c   Gaithersburg, MD  20899
c
c-----------------------------------------------------------------------
c
      dimension x(*), y(*), output(*)
      eps = r1mach(4)
c
      ier = 0
      do 10 i = 1 , 7
         output(i) = 0.0
 10   continue
c
      if (n .le. 2) then
         ier = 1
         return
      endif
c
c      Initialize accumulators
c
c      xsum --  sum of the x observations
c      ysum --  sum of the y observations
c      xysum--  sum of the products of the observations of x and y
c      sxx  --  sum of the squared deviations of x from the mean
c      syy  --  sum of the squared deviations of y from the mean
c      sxy  --  sum of the product of the deviations in x and in y
c
      sxx   = 0.0
      sxy   = 0.0
      syy   = 0.0
      xsum  = 0.0
      ysum  = 0.0
      xysum = 0.0
      x2sum = 0.0
      y2sum = 0.0
      rn    = float(n)
c
c      Accumulate observations, squares, and calculate the means
c
      do 20 i = 1,n
           xsum = xsum + x(i)
           ysum = ysum + y(i)
           xysum = xysum + (x(i) * y(i))
 20   continue
      output(5) = xsum / rn
      output(6) = ysum / rn
      do 30 i = 1,n
           x2sum = x2sum + x(i) * x(i)
           y2sum = y2sum + y(i) * y(i)
           xdev = x(i) - output(5)
           sxx = sxx + xdev ** 2
           ydev = y(i) - output(6)
           syy = syy + ydev ** 2
           sxy = sxy + (xdev * ydev)
 30   continue
      xsum2 = xsum * xsum
      ysum2 = ysum * ysum
      if (sxy .gt. eps) sgnsxy =  1.0
      if (abs(sxy) .le. eps) sgnsxy =  0.0
      if (sxy .lt. -eps) sgnsxy = -1.0
c
c     Set output(1) to the geometric mean slope
c
      denom = x2sum - xsum2/rn
      if (denom .le. eps) then
         ier = 2
         return
      endif
      output(1) = sgnsxy * (sqrt (abs (( y2sum - ysum2 / rn) / denom)))
c
c     Set output(2) to the y-intercept
c
      output(2) = output(6) - output(1) * output(5)
c
c     Set output(7) to the Pearson product-moment correlation coefficien
c
      cornum = xysum - ((xsum * ysum) / rn)
      corden = sqrt ((x2sum - (xsum2 / rn)) *
     1               (y2sum - (ysum2 / rn)))
      if (corden .le. eps) then
         ier = 3
         return
      endif
      output(7) = cornum / corden
c
c     Set output(3) to the standard deviation of the slope
c
      temp = sqrt ((1. - output(7) * output(7)) / rn)
      output(3) = abs(output(1)) * temp
c
c     Set output(4) to the standard deviation of the y-intercept
c
      temp = (output(5) ** 2) * (1. + output(7)) / (sxx / (rn - 1.))
      temp = (syy / (rn - 1.)) * ((1. - output(7)) / rn) * (2. + temp)
      output(4) = sqrt (temp)
      return
      end
