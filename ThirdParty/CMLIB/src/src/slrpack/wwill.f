      subroutine wwill (n, x, y, u, v, bstart, mxiter, delmax, xres,
     1                  yres, ier, eps, iter, a, b, sda,
     2                  sdb, adjsda, adjsdb, xbar, ybar, rmsx, rmsy,
     3                  delx, dely, delz, w, z)
c
c----------------------------------------------------------------------
c
c   Package:  SLRPACK
c
c   Version:  October, 1985
c
c----------------------------------------------------------------------
c
c   PURPOSE
c   -------
c
c   This subroutine is called by user-callable subroutine RWILL and
c   straightforwardly implements the calculations described in the
c   reference to estimate the simple linear regression coefficients
c   when both variables are subject to errors which are not necessarily
c   homogeneous in variance.
c
c   INPUT AND OUTPUT PARAMETERS
c   ---------------------------
c
c   N, X, Y, U, V, BSTART, MXITER, DELMAX, XRES, YRES, IER,
c      EPS -- see subroutine RWILL.
c
c   INTERNAL PARAMETERS
c   -------------------
c
c   DELX    Real vector dimensioned at least N
c           The x deviations from the mean for each observation.
c
c   DELY    Real vector dimensioned at least N
c           The y deviations from the mean for each observation.
c
c   DELZ    Real vector dimensioned at least N
c           The deviation in z for each observation.
c
c   W       Real vector dimensioned at least N
c           The weight of each observation as in the reference.
c
c   Z       Real vector dimensioned at least N
c           The importance of each point in determining the slope
c              of the straight line as in the reference.
c
c   OUTPUT PARAMETERS
c   -----------------
c
c   ITER    Integer scalar
c           The number of iterations the routine needed.
c
c   A       Real scalar
c           The intercept estimate at the current iteration.
c
c   B       Real scalar
c           The slope estimate at the current iteration.
c
c   SDA     Real scalar
c           The standard deviation of the intercept at the current
c              iteration.
c
c   SDB     Real scalar
c           The standard deviation of the slope at the current iteration
c
c   ADJSDA  Real scalar
c           The adjusted standard deviation of the intercept at
c              the current iteration.
c
c   ADJSDB  Real scalar
c           The adjusted standard deviation of the slope at the
c              current iteration.
c
c   XBAR    Real scalar
c           The mean weighted x-observation at the current iteration
c              with the current weights.
c
c   YBAR    Real scalar
c           The mean weighted y-observation at the current iteration
c              with the current weights.
c
c   RMSX    Real scalar
c           The root mean square error in the x-direction.
c
c   RMSY    Real scalar
c           The root mean square error in the y-direction.
c
c   INTERNAL PARAMETER
c   ------------------
c
c   PRKEY   Logical scalar
c
c           PRKEY = T for particular diagnostic printing
c                     (see source code)
c           PRKEY = F otherwise
c           PRKEY is normally set to F
c
c
c----------------------------------------------------------------------
c
      dimension delx(*), dely(*), delz(*), u(*), v(*), w(*), x(*),
     1          xres(*), y(*), yres(*), z(*)
      logical prkey
      prkey = .false.
      pi = 4.0 * atan (1.0)
c
c      Initialize
c
      ier = 0
      iter = 1
      xbar = 0.
      ybar = 0.
      a = 0.
      b = bstart
      sda = 0.
      sdb = 0.
      adjsda = 0.
      adjsdb = 0.
      rmsx = 0.
      rmsy = 0.
      angle2 = atan(bstart)
      if (angle2 .le. 0.) angle2 = angle2 + pi
 10   bprev = b
      aprev = a
      xbarp = xbar
      ybarp = ybar
      angle1 = angle2
      iter = iter + 1
c
c      Compute weighted averages of independent variable and
c      dependent variable observations
c
      sumx = 0.
      sumy = 0.
      sumw = 0.
      b2 = bprev * bprev
      do 20 i = 1,n
          w(i) = 1. / (v(i) + u(i) * b2)
          sumx = sumx + w(i) * x(i)
          sumy = sumy + w(i) * y(i)
   20     sumw = sumw + w(i)
      if (prkey) write(6,900) (i,w(i), i = 1,n)
  900 format(' weight w(',i4,') = ',f8.3)
      if (abs(sumw) .le. eps) then
         ier = 4
         a = 0.
         b = 0.
         xbar = 0.
         ybar = 0.
         iter = 0
         return
      endif
      xbar = sumx / sumw
      ybar = sumy / sumw
      if (prkey) write(6,901) xbar,ybar,iter
  901 format(//' xbar = ',f8.3,', ybar = ',f8.3,' on iteration ',i4)
c
c      Compute slope and zbar as in Williamson (1968)
c
      sum1 = 0.
      sum2 = 0.
      sumz = 0.
      do 30 i = 1,n
          delx(i) = x(i) - xbar
          dely(i) = y(i) - ybar
          z(i) = w(i) * (v(i) * delx(i) + bprev * u(i) * dely(i))
          wz = w(i) * z(i)
          sumz = sumz + wz
          sum1 = sum1 + wz * dely(i)
   30     sum2 = sum2 + wz * delx(i)
      zbar = sumz / sumw
      if (abs(sum2) .le. eps) then
         ier = 2
         a = 0.
         b = 0.
         xbar = 0.
         ybar = 0.
         iter = 0
         return
      endif
      b = sum1 / sum2
      a = ybar - b * xbar
      if (prkey) write(6,902) b,a,iter
  902 format('    b = ',f8.3,',    a = ',f8.3,' on iteration ',i4//)
      if (iter .ge. mxiter) then
         ier = 7
         goto 40
      endif
c
c      Compare angle of the previous regression line (angle1) and angle
c      angle of the current regression line (angle2)
c
      angle2 = atan2(sum1, sum2)
      if (angle2 .le. 0.) angle2 = pi + angle2
      angle = abs(angle1 - angle2)
      if (angle .ge. delmax) go to 10
c
c       Compute output
c
c      Calculate uncertainties in slope and intercept
c
   40 s = 0.
      qi = 0.
      sum3 = 0.
      if (abs(b) .le. eps) then
         ier = 5
         return
      endif
      do 50 i = 1,n
          delz(i) = z(i) - zbar
          qi = qi + (w(i) * (delx(i) * dely(i) / b + 4. * delz(i)
     1       * (z(i) - delx(i))))
          s = s + w(i) * (b*delx(i) - dely(i)) ** 2
   50     sum3 = sum3 + w(i)**2 * (delx(i)**2*v(i) + dely(i)**2*u(i))
      qi2 = qi * qi
      if (abs(qi2) .le. eps) then
         ier = 6
         return
      endif
      varb = sum3 / qi2
      vara = 1./sumw + 2.*(xbar + 2.*zbar) * zbar/qi +
     1    ((xbar + 2.*zbar) ** 2) * varb
      sda = sqrt(vara)
      sdb = sqrt(varb)
      adjsda = sqrt(s/float(n - 2) * vara)
      adjsdb = sqrt(s/float(n - 2) * varb)
c
c       Compute residual information
c
      xr2 = 0.
      yr2 = 0.
      do 60 i = 1,n
          xres(i) = w(i) * (v(i)*x(i) + b*u(i)*(y(i) - a)) - x(i)
          yres(i) = (xres(i) * v(i)) / ((-1.) * b * u(i))
          xr2 = xr2 + xres(i) ** 2
   60     yr2 = yr2 + yres(i) ** 2
      rmsx = sqrt(xr2 / float(n))
      rmsy = sqrt(yr2 / float(n))
      return
      end
