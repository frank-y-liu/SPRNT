      subroutine rybar (n, x, y, wx, wy, bprev, w, prkey, xbar, ybar)
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
c     This subroutine is used by user-callable subroutine RYORK to
c     calculate the weighted averages of the observations.
c
c-----------------------------------------------------------------------
c
      dimension w(*), wx(*), wy(*), x(*), y(*)
      logical prkey
c
c      Initialize accumulators
c
      xnum = 0.
      ynum = 0.
      denom = 0.
c
c      Compute numerators and denominator of weighted averages of
c      observations
c
      b2 = bprev * bprev
      do 100 i = 1, n
           w(i) = (wx(i) * wy(i)) / ((b2 * wy(i)) + wx(i))
           xnum = xnum + (w(i) * x(i))
           ynum = ynum + (w(i) * y(i))
           denom = denom + w(i)
  100 continue
c
c      Compute weighted averages
c
      xbar = xnum / denom
      ybar = ynum / denom
      if (prkey) write (6,101) xbar,ybar
  101 format (' mean of x- and y-observations, weighted ',/,
     *        ' according to equation (19) in York, are: ',/,
     *        ' xbar = ', f8.3/ ' ybar = ', f8.3,//)
      return
      end
