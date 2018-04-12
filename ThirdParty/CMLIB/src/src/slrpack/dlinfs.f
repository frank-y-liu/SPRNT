 
      subroutine dlinfs(n, x, y, beta1, beta2, lambda, kount, p, q, r,
     1                  ifault)
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
c   This user-callable subroutine solves the simple model,
c   y = beta1 + (beta2)x, under the Chebychev norm criterion.
c
c   DESCRIPTION
c   -----------
c
c   1.  The Y's are observed values of the dependent variable and the
c       X's are the observed values of the independent variables, the
c       predictor variables.
c
c   2.  The general Chebychev problem is to minimize LAMBDA where
c
c            LAMBDA = MAX( |Y(i) - BETA1 - X(i)*BETA2| ) , over all i,
c
c       such that
c
c            Y(i)-LAMBDA <= BETA1 + X(i)*BETA2 <= Y(i)+LAMBDA
c
c   3.  The optimal value of (BETA1, BETA2)' will minimize the
c       maximum deviation and LAMBDA will be the value of this deviation
c
c   INPUT PARAMETERS
c   ----------------
c
c   N      Integer scalar (unchanged on output).
c          Number of observations.
c
c          N must be at least 3.
c
c   X      Double precision vector dimensioned at least N (unchanged
c               on output).
c          Observed values of the independent variable.
c
c   Y      Double precision vector dimensioned at least N (unchanged
c               on output).
c          Observed values of the dependent variable.
c
c   OUTPUT PARAMETERS
c   -----------------
c
c   BETA1      Double precision scalar.
c              The estimated intercept term for the fitted line.
c
c   BETA2      Double precision scalar.
c              The estimated slope term for the fitted line.
c
c   LAMBDA     Double precision scalar.
c              The least maximum absolute deviation, also
c              the objective function value.
c
c   P          Integer scalar.
c              Indicates which row of X is the p-th constraint.
c
c   Q          Integer scalar.
c              Indicates which row of X is the q-th constraint.
c
c   R          Integer scalar.
c              Indicates which row of X is the r-th constraint.
c
c   KOUNT      Integer scalar.
c              The number of iterations.
c
c   IFAULT     Integer scalar.
c              Failure indicator.
c
c              = 0    normal termination
c
c              = 1    the data set is too small
c                     cannot compute BETA1, BETA2, LAMBDA, P, Q, R
c
c              = 2    all x values are equal
c                     cannot compute BETA1, BETA2, LAMBDA, P, Q, R
c
c              = 3    the data are collinear, all constraints are feasib
c                     cannot compute P, Q, R
c
c   PRECISION
c   ---------
c
c   All calculations are done in double precision.
c
c   LANGUAGE
c   --------
c
c   All code is written in standard Fortran 77.
c
c   OTHER SUBROUTINES
c   -----------------
c
c   PORT  subroutine  R1MACH
c
c   REFERENCE
c   ---------
c
c   Sklar, M. and R. Armstrong (1983). "A Linear Programming Algorithm
c        for the Simple Model for Discrete Chebychev Curve Fitting",
c        Computers and Operations Research 10, 237-248.
c
c   Sklar, M. and R. Armstrong (1984). "An Algorithm for Discrete
c      Chebychev Curve Fitting for the Simple Model Using A Dual
c      Linear Programming Approach", Communications in Statistics:
c      Simulation and Computation, 13(4), 555-569.
c
c   AUTHORS
c   -------
c
c   Michael G. Sklar and Ronald D. Armstrong
c   University of Georgia
c
c   NBS CONTACT
c   -----------
c
c   Sally E. Howe
c   Scientific Computing Division
c   Center for Applied Mathematics
c   National Engineering Labratory
c   301-921-3395
c
c---------------------------------------------------------------------
c
      integer           ifault,irow,isave,kount,n
      integer           p,q,r,s
c
c---------------------
c   internal variables
c---------------------
c
      double precision  acu,alpha1,alpha2,beta1,beta2,big,d,delhat
      double precision  delstr,devsgn,div,dsubi,dsubr,dsubs,lambda
      double precision  ratio,r1,r2,save,signp,signq,signr
      double precision  signs,s1,s2,sumr,sums,top
      double precision  x(*), y(*), d1mach
c
c   signp, signq, signr, signs  Indicators for the p-th, q-th, r-th,
c                               and s-th constraints, respectively.
c
c   sign(k)  =  1  If the k-th constraint is at or violates its
c                  upper bound.
c            = -1  If it is at or violates its lower bound.
c
c   dsubi   Amount by which the i-th constraint is violated.
c
c   dsubr   Amount by which the r-th constraint is violated.
c
c   dsubs   Amount by which the s-th constraint is violated.
c
c   delhat  The smallest delta necessary to bring constraints into
c           feasibility.
c
c   delstr  The largest tolerable delta allowed before a feasible
c           constraint moves beyond its opposite bound.
c
c   devsgn  The signed amount violated in the s-th constraint.
c
c   alpha1, alpha2  Y-coordinates of p-th and q-th constraints.
c
c   acu     Machine-dependent constant which is used to test for zero.
c           ACU is set to 100 times the relative machine accuracy.
c
c   big     Machine-dependent constant which is used in determining the
c           minimum ratio.  BIG is set to the largest floating point
c           number available.
c
      kount=0
      beta1=0.0d0
      beta2=0.0d0
      lambda=0.0d0
      p=0
      q=0
      r=0
      acu = d1mach(4) * 100.0
      big = d1mach(2)
      ifault=0
c
c   check problem size
c
      if (n .le. 2) then
         ifault = 1
         return
      endif
c
c     get p-th and q-th constraints and
c     check that not all x-values are equal
c
      p=1
      do 110 i=2,n
          if(dabs(x(p)-x(i)).gt.acu)goto 120
  110 continue
c
      ifault=2
      p=0
      return
c
  120 q=i
c
c     step 1 - compute initial (trivial) solution
c
      s=0
      d=1./(x(q)-x(p))
      alpha1=y(p)
      alpha2=y(q)
      beta1=d*(x(q)*alpha1-x(p)*alpha2)
      beta2=d*(alpha2-alpha1)
c
c     step 2 - add a constraint (labeled r-th) to the problem
c
      irow=p
  130 if(irow.eq.n) then
         ifault=3
         p=0
         q=0
         return
      endif
      irow=irow+1
c
c     check if r-th constraint is feasible
c     (i.e. if p-th, q-th, and r-th constraints are not collinear)
c
      dsubr=beta1+x(irow)*beta2-y(irow)
      if(dabs(dsubr).lt.acu)goto 130
      signr=dsign(1.d0,dsubr)
      ifault=0
      r=irow
c
c     adjust for the r-th constraint
c
      r1=d*(x(q)-x(r))
      r2=d*(x(r)-x(p))
      signp=dsign(1.d0,-signr*r1)
      signq=dsign(1.d0,-signr*r2)
      sumr=signr-signp*r1-signq*r2
      lambda=dabs(dsubr/sumr)
      alpha1=y(p)+signp*lambda
      alpha2=y(q)+signq*lambda
      beta1=d*(x(q)*alpha1-x(p)*alpha2)
      beta2=d*(alpha2-alpha1)
c
c     step 3 - find a fourth constraint (labeled s-th) which most
c              violates the solution
c
  140 s=0
      dsubs=acu
      do 150 i=1,n
          dsubi=dabs(beta1+beta2*x(i)-y(i))-lambda
          if(dsubi.le.dsubs) go to 150
          dsubs=dsubi
          s=i
  150 continue
c
c     check for termination - no constraints are violated
c
      if(s.eq.0) go to 240
c
c     compute characteristics of most-violated (s-th) constraint
c
      devsgn=beta1+beta2*x(s)-y(s)
      signs=dsign(1.d0,devsgn)
      s1=d*(x(q)-x(s))
      s2=d*(x(s)-x(p))
      sums=-signs+signq*s2 +signp*s1
c
c     step 4 - determine which constraint will not be binding when the
c              s-th constraint is made feasible (by locating the
c              minimum ratio)
c
  160 ratio=big
c
c     if the signs are not the same, the p-th constraint is not
c     a candidate
c
      if(signs*dsign(1.d0,s1).ne.signp)goto 170
c
c     ignore ratio if denominator is effectively zero
c
      if(dabs(s1).lt.acu)goto 170
      ratio=dabs(r1/s1)
c
c     if the signs are not the same, the q-th constraint is not
c     a candidate
c
  170 if(signs*dsign(1.d0,s2).ne.signq)goto 180
c
c     ignore ratio if denominator is effectively zero
c
      if(dabs(s2).lt.acu)goto 180
c
c     if minimum is p-th, go to process p-th
c
      if(dabs(r2/s2).lt.ratio)goto 200
      go to 210
c
c     if minimum is q-th, go to switch the p-th and q-th constraints
c     so that the algorithm will always process the p-th constraint
c
  180 if(ratio.lt.big)goto 210
c
c     step 5 - no minimum was located, so process the r-th constraint
c
      delhat=dabs(dsubs/sums)
c
c     calculate the largest tolerable delta
c
      div=dabs(sumr)-2.
      if(div.lt.acu) go to 190
      delstr=2.*lambda/div
      if(delstr.ge.delhat) go to 190
c
c     step 5a - as the s-th becomes feasible, the r-th becomes
c               infeasible by moving outside a bound, so switch
c               r-th and s-th constraint indicators
c
      save=sums
      sums=-sumr+signr+signr
      sumr=-save
      save=signr
      signr=signs
      signs=-save
      dsubs=dabs(sums*delhat) -2.*lambda
      lambda=lambda+delhat
      alpha1=y(p)+signp*lambda
      alpha2=y(q)+signq*lambda
      beta1=d*(x(q)*alpha1-x(p)*alpha2)
      beta2=d*(alpha2-alpha1)
      save=r1
      r1=s1
      s1=save
      save=r2
      r2=s2
      s2=save
      isave=r
      r=s
      s=isave
      go to 160
c
c     step 5b - replace the r-th constraint with the s-th constraint,
c               since the s-th became feasible and the r-th moved
c               interior (between bounds)
c
  190 signr=signs
      r1=s1
      r2=s2
      sumr=-sums
      r=s
      go to 230
c
c     step 6 - for efficiency, the algorithm always processes the
c              p-th constraint.  if the q-th is to be processed,
c              switch the p-th and q-th constraints
c
  200 isave=p
      p=q
      q=isave
      save=signp
      signp=signq
      signq=save
      save=alpha1
      alpha1=alpha2
      alpha2=save
      save=r1
      r1=r2
      r2=save
      save=s1
      s1=s2
      s2=save
c
c     step 7 - process the movement of the p-th constraint
c
  210 delhat=dabs(r1*dsubs/(r1*sums+s1*sumr))
      top=-2.*lambda*r1
      div=r1 +r1 +signp*sumr
      if(dsign(1.d0,top).ne.dsign(1.d0,div)) go to 220
      if(dabs(div).lt.acu) go to 220
      delstr=top/div
c
c     check to see if p-th constraint swings across to its
c     opposite bound
c
      if(delstr.ge.delhat) go to 220
c
c     step 7a - the p-th constraint moves to its opposite bound,
c               and the s-th constraint is still infeasible.
c               adjust for the movement, then go to step 4 to
c               locate the next smallest ratio
c
      lambda=lambda+delstr
      dsubs=dsubs-delstr*dabs(sums+s1*sumr/r1)
      sumr=sumr+2.*signp*r1
      sums=sums-2.*signp*s1
      signp=-signp
      alpha1=y(p)+signp*lambda
      alpha2=y(q)+signq*lambda
      go to 160
c
c     step 7b - the s-th became feasible and the p-th moved interior
c               so replace the p-th constraint with the s-th constraint
c
  220 signp=signs
      r1=r1/s1
      r2=r2-s2*r1
      sumr=signr-signp*r1-signq*r2
      p=s
      kount=kount+1
c
c     update lambda, compute the new solution - go to step 3 to
c     locate a new most-violated constraint
c
  230 s=0
      lambda=lambda+delhat
      alpha1=y(p)+signp*lambda
      alpha2=y(q)+signq*lambda
      d=1./(x(q)-x(p))
      beta1=d*(x(q)*alpha1-x(p)*alpha2)
      beta2=d*(alpha2-alpha1)
      go to 140
  240 continue
      return
      end
