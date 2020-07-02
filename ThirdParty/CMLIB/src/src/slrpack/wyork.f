      subroutine wyork (n, x, y, wx, wy, bstart, mxiter, delmax, ires,
     1                  xres, yres, ier, eps, angle, u, u2, uv, v, v2,
     2                  w, w3, wwx, wwy, a, b, xbar, ybar, sdb, sda,
     3                  iter, phi)
c
c-----------------------------------------------------------------------
c
c     Package:  SLRPACK
c
c     Version:  October, 1985
c
c-----------------------------------------------------------------------
c
c   PURPOSE
c   -------
c
c   This subroutine is called by subroutine RYORK to estimate the
c   simple linear regression coefficients when both variables are
c   subject to errors which are not necessarily homogeneous in
c   variance.
c
c   DESCRIPTION
c   -----------
c
c   1. This subroutine straightforwardly implements the calculations
c      described in the reference.
c
c   2. In his paper, York discusses the problem of selecting one of the
c      at most three roots of the cubic equation as the desired slope.
c      In this implementation, at each iteration, the root selected as
c      the current slope is the closest of the at most three roots to
c      the root at the previous iteration.
c
c   INPUT PARAMETERS
c   ----------------
c
c   N, X, Y, WX, WY, BSTART, MXITER, DELMAX, IRES, XRES, YRES,
c      IER, EPS -- see subroutine RYORK.
c
c   INTERNAL PARAMETERS
c   -------------------
c
c   ANGLE   Real vector dimensioned at least 4
c           Angles of the candidate regression lines at the current
c              iteration.
c
c   U       Real vector dimensioned at least N
c           The x-deviations from the mean for each observation.
c
c   U2      Real vector dimensioned at least N
c           The square of each x-deviation.
c
c   UV      Real vector dimensioned at least N
c           The product of the x- and y-deviations.
c
c   V       Real vector dimensioned at least N
c           The y-deviations from the mean for each observation.
c
c   V2      Real vector dimensioned at least N
c           The square of each y-deviation.
c
c   W       Real vector dimensioned at least N
c           The W vector described in equation 14 in reference.
c
c   W3      Real vector dimensioned at least N
c           The product of the W vector squared and the x variance
c              for each observation.
c
c   WWX     Real vector dimensioned at least N
c           The product of the W vector and the x variance for each
c              observation.
c
c   WWY     Real vector dimensioned at least N
c           The product of the W vector and the y variance for each
c              observation.
c
c   PHI     Real vector dimensioned at least 3
c           The roots of the cubic equation.
c
c   OUTPUT PARAMETERS
c   -----------------
c
c   A       Real scalar
c           The intercept at the current iteration.
c
c   B       Real vector dimensioned at least 3
c           The slopes at the current iteration.
c
c   XBAR    Real scalar
c           The mean x-observation.
c
c   YBAR    Real scalar
c           The mean y-observation.
c
c   SDA     Real scalar
c           The standard deviation of the intercept.
c
c   SDB     Real scalar
c           The standard deviation of the slope.
c
c   ITER    Integer scalar
c           The number of iterations the routine needed.
c
c   INTERNAL PARAMETER
c   ------------------
c
c   PRKEY   Logical scalar
c
c           PRKEY = T for particular diagnostic printing
c                        (see source code)
c           PRKEY = F otherwise
c
c           PRKEY is normally set to F.
c
c-----------------------------------------------------------------------
c
      dimension angle(*), b(*), u(*), u2(*), uv(*), v(*),
     1   v2(*), w(*), w3(*), wwx(*), wwy(*), wx(*), wy(*), x(*),
     2   xres(*), y(*), yres(*), phi(*)
      logical prkey
c
c-------------
c   Initialize
c-------------
c
      ier = 0
      pi = 4.0 * atan(1.0)
      third = 1./3.
      rn = float(n)
      prkey = .false.
      iter = 1
      b(3) = bstart
      xbarp = 0.
      ybarp = 0.
      bprev = b(3)
      xbar = 0.
      ybar = 0.
      do 20 i = 1,n
          xbar = xbar + x(i)
          ybar = ybar + y(i)
 20   continue
      xbar = xbar / rn
      ybar = ybar / rn
      if (prkey) write (6,898) xbar, ybar
  898 format (' unweighted means are: '/,
     1        ' xbar = ',f8.3 / ' ybar = ',f8.3/)
c
c----------------------------
c   Iteratively compute slope
c----------------------------
c
 30   iter = iter + 1
      if (prkey) write (6,899) iter
  899 format (//' begin iteration', i4/)
      alpha = 0.
      beta = 0.
      gamma = 0.
      denom = 0.
      wbuv = 0.
      wi = 0.
      wu2 = 0.
      wx2 = 0.
c
c      Compute weighted averages of independent variable and dependent
c      variable observations.
c
      call rybar (n, x, y, wx, wy, bprev, w, prkey, xbar, ybar)
c
c      Compute alpha, beta, and gamma as in York (1966)
c
      do 40 i = 1,n
           u(i) = x(i) - xbar
           v(i) = y(i) - ybar
           u2(i) = u(i) * u(i)
           v2(i) = v(i) * v(i)
           uv(i) = u(i) * v(i)
           wwx(i) = w(i) / wx(i)
           wwy(i) = w(i) / wy(i)
           w3(i) = w(i) * wwx(i)
           denom = denom + w3(i) * u2(i)
           alpha = alpha + w3(i) * uv(i)
           beta = beta + w3(i) * v2(i) - w(i) * u2(i)
           gamma = gamma + w(i) * uv(i)
  40  continue
      if (abs(denom) .le. eps) then
         ier = 2
         b(3) = 0.
         iter = 0
         return
      endif
      alpha = alpha * 2. / (denom * 3.)
      beta = beta / (denom * 3.)
      gamma = -gamma / denom
      if (prkey) write (6,901) alpha,beta,gamma
  901 format (' alpha, beta, and gamma (York, p. 1084):',/
     1        ' alpha = ', f10.5,/,
     2        '  beta = ', f10.5,/,
     3        ' gamma = ', f10.5,//)
c
c      Compute real solution(s) to the cubic equation
c
c      (b**3) - 3 * alpha * (b**2) + 3 * beta * b - gamma = 0
c
c      where b = bprev
c
c      First convert cubic equation to normal form
c
c      (x**3) + an * (x) + bn = 0
c
c      via x = b - alpha
c
      alpha2 = alpha ** 2
      alpha3 = alpha ** 3
      an = -3. * (alpha2 - beta)
      bn = -2. * alpha3 + 3. * alpha * beta - gamma
      an3 = an ** 3
      bn2 = bn ** 2
c
c      Determine number of real roots of cubic equation
c
      testn = bn2/4. + an3/27.
      if (prkey) write (6,902) testn
  902 format (' number of real roots of cubic equation ',/,
     1        ' (York p. 1084) determined by testn = ',f8.3,//)
      if (testn.gt.0.) go to 60
      if (testn.eq.0.) go to 70
c
c---------------------------------------------
c     Case 1.  There are three real roots.
c
c              Compute roots in b(i), i = 1,3.
c---------------------------------------------
c
      a2br = sqrt (alpha2 - beta)
      term = (alpha3 - (1.5 * alpha * beta) + (.5 * gamma))
     1        / (a2br ** 3)
      phi(1) = acos (term)
      phi(2) = phi(1) + (2. * pi)
      phi(3) = phi(2) + (2. * pi)
      if (prkey) write (6,903) phi(1),phi(2),phi(3)
  903 format (' case 1 - there are three real roots -',/,
     1        ' phi-s (York, p. 1084)  are:',/
     2        ' phi(1) = ', f8.3,/,
     3        ' phi(2) = ', f8.3,/,
     4        ' phi(3) = ', f8.3,//)
c
c     Compute roots
c
      do 50 i = 1,3
         b(i) = alpha + 2. * a2br * cos(phi(i)/3.)
  50  continue
      go to 80
c
c-------------------------------------------------------
c     Case 2.  There is one real root (and two conjugate
c              imaginary roots).
c
c              Put real root in b(3).
c-------------------------------------------------------
c
  60  testnr = sqrt(testn)
      b(1) = 0.
      b(2) = 0.
      temp = -(bn/2.) + testnr
      temp3 = abs(temp) ** third
      term2 = sign (temp3,temp)
      temp = -(bn/2.) - testnr
      temp3 = abs(temp) ** third
      term3 = sign (temp3,temp)
      b(3) = alpha + term2 + term3
      a = ybar - b(3) * xbar
      if (prkey) write (6,904) iter,b(3)
  904 format (' case 2 - there is one real root, b(3,',
     1        i4,') = ', f8.3//)
      go to 100
c
c--------------------------------------------------------------
c     Case 3.  There are three real roots of which at least two
c              are equal.
c
c              Compute roots in b(i), i = 1,3.
c--------------------------------------------------------------
c
 70   temp = (- bn / 2.) ** third
      b(3) = alpha - term
      b(2) = b(3)
      b(1) = alpha + 2. * term
      if (prkey) write (6,905) (b(i), i = 1,3)
  905 format (' case 3 - there are three real roots, ',/,
     1        '          two of which are equal: ',/,
     2        ' b(1) = ',f8.3,/,
     3        ' b(2) = ',f8.3,/,
     4        ' b(3) = ',f8.3,//)
c
c--------------------------
c     Complete computations
c--------------------------
c
c      Arrange three solutions so line with slope b(3) and
c      y-intercept a lies closest to line with slope bprev
c      passing through (xbarp,ybarp)
c
  80  call ryarng (b, bprev, angle, prkey)
      a = ybar - b(3) * xbar
      if (prkey) write (6,906) iter,b(3),a
  906 format (//' iteration ', i4, ' results are:'/
     1        '    slope    intercept'/
     2        f10.5, f14.5//)
c
c----------------
c     Termination
c----------------
c
c     Check for termination
c
  100 if (iter .ge. mxiter) then
         ier = 4
         goto 140
      endif
      if (angle(3) .gt. delmax) then
         xbarp = xbar
         ybarp = ybar
         bprev = b(3)
         goto 30
      endif
c
c     Regression line slope estimates have converged.
c
c     Calculate variances and residuals if ires is not zero.
c
  140 do 150 i = 1,n
         wbuv = wbuv + w(i) * ((b(3) * u(i) - v(i)) ** 2)
         wu2 = wu2 + w(i) * u2(i)
         wx2 = wx2 + w(i) * x(i) * x(i)
         wi = wi + w(i)
         if (ires .ne. 0) then
            temp = a + b(3) * x(i) - y(i)
            xres(i) = -b(3) * wwx(i) * temp
            yres(i) = wwy(i) * temp
         endif
  150 continue
      sdb = sqrt (wbuv / ((rn - 2.) * wu2))
      sda = sdb * sqrt (wx2 / wi)
      return
      end
