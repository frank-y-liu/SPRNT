      SUBROUTINE SDASSL (RES,NEQ,T,Y,YPRIME,TOUT,
     *  INFO,RTOL,ATOL,IDID,
     *  RWORK,LRW,IWORK,LIW,RPAR,IPAR,
     *  JAC)
C
C***BEGIN PROLOGUE  SDASSL
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C           Added PARAMETER statement for pointers
C           into IWORK.
C***CATEGORY NO.  D2A2
C***KEYWORDS  DIFFERENTIAL/ALGEBRAIC,BACKWARD DIFFERENTIATION FORMULAS
C             IMPLICIT DIFFERENTIAL SYSTEMS
C***AUTHOR  PETZOLD,LINDA R.
C             APPLIED MATHEMATICS DIVISION 8331
C             SANDIA NATIONAL LABORATORIES
C             LIVERMORE, CA.    94550
C***PURPOSE  DIFFERENTIAL/ALGEBRAIC SYSTEM SOLVER
C***DESCRIPTION
C  ---------------------------------------------------------------------
C
C  this code solves a system of differential/
C  algebraic equations of the form
C  g(t,y,yprime) = 0.
C
C  subroutine sdassl uses the backward
C  differentiation formulas of orders one
C  through five to solve a system of the above
C  form for y and yprime. values for y
C  and yprime at the initial time must
C  be given as input. these values must
C  be consistent, (that is. if t,y,yprime
C  are the given initial values, they must
C  satisfy g(t,y,yprime) = 0.)
C  the subroutine solves the system from t to tout. it is
C  easy to continue the solution to get results
C  at additional tout. this is the interval
C  mode of operation. intermediate results can
C  also be obtained easily by using the intermediate-
C  output capability.
C
C
C  ------------description of arguments to sdassl-----------------------
C  ------------(an overview)--------------------------------------------
C
C  the parameters are
C
C  res -- this is a subroutine which you provide
C         to define the differential/algebraic
C         system
C
C  neq -- this is the number of equations
C         to be solved
C
C  t -- this is the current value of the
C       independent variable.
C
C  tout -- this is a point at which a solution
C      is desired.
C
C  info(*) -- the basic task of the code is
C             to solve the system from t to
C             tout and return an answer at tout.
C             info(*) is an integer array which is
C             used to communicate exactly how you
C             want this task to be carried out.
C
C  y(*) -- this array contains the solution
C          components at t
C
C  yprime(*) -- this array contains the derivatives
C               of the solution components at t
C
C  rtol,atol -- these quantities represent
C               absolute and relative error
C               tolerances which you provide to indicate
C               how accurately you wish the solution
C               to be computed. you may choose them
C               to be both scalars or else both
C               vectors.
C
C  idid -- this scalar quantity is an indicator reporting
C          what the code did. you must monitor this
C          integer variable to decide what action to
C          take next.
C
C  rwork(*),lrw -- rwork(*) is a real work array of
C                  length lrw which provides the code
C                  with needed storage space.
C
C  iwork(*),liw -- iwork(*) is an integer work array
C                  of length liw which provides the code
C                  with needed storage space.
C
C  rpar,ipar -- these are real and integer parameter
C               arrays which you can use for
C               communication between your calling
C               program and the res subroutine
C               (and the jac subroutine)
C
C  jac -- this is the name of a subroutine which you
C         may choose to provide for defining
C         a matrix of partial derivatives
C         described below.
C
C  quantities which are used as input items are
C     neq,t,y(*),yprime(*),tout,info(*),
C     rtol,atol,rwork(1),rwork(2),rwork(3),lrw,iwork(1),
C     iwork(2),iwork(3),and liw.
C
C  quantities which may be altered by the code are
C     t,y(*),yprime(*),info(1),rtol,atol,
C     idid,rwork(*) and iwork(*)
C
C  ----------input-what to do on the first call to sdassl---------------
C
C
C  the first call of the code is defined to be the start of each new
C  problem. read through the descriptions of all the following items,
C  provide sufficient storage space for designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  give information about how you want the problem to be solved.
C
C
C  res -- provide a subroutine of the form
C             subroutine res(t,y,yprime,delta,ires,rpar,ipar)
C         to define the system of differential/algebraic
C         equations which is to be solved. for the given values
C         of t,y and yprime, the subroutine should
C         return the residual of the differential/algebraic
C         system
C             delta = g(t,y,yprime)
C         (delta(*) is a vector of length neq which is
C         output for res.)
C
C         subroutine res must not alter t,y or yprime.
C         you must declare the name res in an external
C         statement in your program that calls sdassl.
C         you must dimension y,yprime and delta in res.
C
C         ires is an integer flag which is always equal to
C         zero on input.  subroutine res should alter ires
C         only if it encounters an illegal value of y or
C         a stop condition.  set ires = -1 if an input value
C         is illegal, and sdassl will try to solve the problem
C         without getting ires = -1.  if ires = -2, sdassl
C         will return control to the calling program
C         with idid = -11.
C
C         rpar and ipar are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and subroutine res. they are not altered by sdassl. if you
C         do not need rpar or ipar, ignore these parameters by treat-
C         ing them as dummy arguments. if you do choose to use them,
C         dimension them in your calling program and in res as arrays
C         of appropriate length.
C
C  neq -- set it to the number of differential equations.
C         (neq .ge. 1)
C
C  t -- set it to the initial point of the integration.
C       t must be defined as a variable.
C
C  y(*) -- set this vector to the initial values of the neq solution
C          components at the initial point. you must dimension y of
C          length at least neq in your calling program.
C
C  yprime(*) -- set this vector to the initial values of
C               the neq first derivatives of the solution
C               components at the initial point. you
C               must dimension yprime at least neq
C               in your calling program.  if you do not
C               know initial values of some of the solution
C               components, see the explanation of info(11).
C
C  tout - set it to the first point at which a solution
C         is desired. you can not take tout = t.
C         integration either forward in t (tout .gt. t) or
C         backward in t (tout .lt. t) is permitted.
C
C         the code advances the solution from t to tout using
C         step sizes which are automatically selected so as to
C         achieve the desired accuracy. if you wish, the code will
C         return with the solution and its derivative at
C         intermediate steps (intermediate-output mode) so that
C         you can monitor them, but you still must provide tout in
C         accord with the basic aim of the code.
C
C         the first step taken by the code is a critical one
C         because it must reflect how fast the solution changes near
C         the initial point. the code automatically selects an
C         initial step size which is practically always suitable for
C         the problem. by using the fact that the code will not step
C         past tout in the first step, you could, if necessary,
C         restrict the length of the initial step size.
C
C         for some problems it may not be permissable to integrate
C         past a point tstop because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         tstop. when you have declared a tstop point (see info(4)
C         and rwork(1)), you have told the code not to integrate
C         past tstop. in this case any tout beyond tstop is invalid
C         input.
C
C  info(*) - use the info array to give the code more details about
C            how you want your problem solved. this array should be
C            dimensioned of length 15, though sdassl uses
C            only the first nine entries. you must respond to all of
C            the following items which are arranged as questions. the
C            simplest use of the code corresponds to answering all
C            questions as yes ,i.e. setting all entries of info to 0.
C
C       info(1) - this parameter enables the code to initialize
C              itself. you must set it to indicate the start of every
C              new problem.
C
C          **** is this the first call for this problem ...
C                yes - set info(1) = 0
C                 no - not applicable here.
C                      see below for continuation calls.  ****
C
C       info(2) - how much accuracy you want of your solution
C              is specified by the error tolerances rtol and atol.
C              the simplest use is to take them both to be scalars.
C              to obtain more flexibility, they can both be vectors.
C              the code must be told your choice.
C
C          **** are both error tolerances rtol, atol scalars ...
C                yes - set info(2) = 0
C                      and input scalars for both rtol and atol
C                 no - set info(2) = 1
C                      and input arrays for both rtol and atol ****
C
C       info(3) - the code integrates from t in the direction
C              of tout by steps. if you wish, it will return the
C              computed solution and derivative at the next
C              intermediate step (the intermediate-output mode) or
C              tout, whichever comes first. this is a good way to
C              proceed if you want to see the behavior of the solution.
C              if you must have solutions at a great many specific
C              tout points, this code will compute them efficiently.
C
C          **** do you want the solution only at
C                tout (and not at the next intermediate step) ...
C                 yes - set info(3) = 0
C                  no - set info(3) = 1 ****
C
C       info(4) - to handle solutions at a great many specific
C              values tout efficiently, this code may integrate past
C              tout and interpolate to obtain the result at tout.
C              sometimes it is not possible to integrate beyond some
C              point tstop because the equation changes there or it is
C              not defined past tstop. then you must tell the code
C              not to go past.
C
C           **** can the integration be carried out without any
C                restrictions on the independent variable t ...
C                 yes - set info(4)=0
C                  no - set info(4)=1
C                       and define the stopping point tstop by
C                       setting rwork(1)=tstop ****
C
C       info(5) - to solve differential/algebraic problems it is
C              necessary to use a matrix of partial derivatives of the
C              system of differential equations.  if you do not
C              provide a subroutine to evaluate it analytically (see
C              description of the item jac in the call list), it will
C              be approximated by numerical differencing in this code.
C              although it is less trouble for you to have the code
C              compute partial derivatives by numerical differencing,
C              the solution will be more reliable if you provide the
C              derivatives via jac. sometimes numerical differencing
C              is cheaper than evaluating derivatives in jac and
C              sometimes it is not - this depends on your problem.
C
C           **** do you want the code to evaluate the partial
C                  derivatives automatically by numerical differences ..
C                   yes - set info(5)=0
C                    no - set info(5)=1
C                  and provide subroutine jac for evaluating the
C                  matrix of partial derivatives ****
C
C       info(6) - sdassl will perform much better if the matrix of
C              partial derivatives, dg/dy + cj*dg/dyprime,
C              (here cj is a scalar determined by sdassl)
C              is banded and the code is told this. in this
C              case, the storage needed will be greatly reduced,
C              numerical differencing will be performed much cheaper,
C              and a number of important algorithms will execute much
C              faster. the differential equation is said to have
C              half-bandwidths ml (lower) and mu (upper) if equation i
C              involves only unknowns y(j) with
C                             i-ml .le. j .le. i+mu
C              for all i=1,2,...,neq. thus, ml and mu are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded. if you do not
C              indicate that the equation has a banded matrix of partial
C                 derivatives
C              the code works with a full matrix of neq**2 elements
C              (stored in the conventional way). computations with
C              banded matrices cost less time and storage than with
C              full matrices if  2*ml+mu .lt. neq.  if you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide subroutine jac to
C              compute the partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of jac.
C
C          **** do you want to solve the problem using a full
C               (dense) matrix (and not a special banded
C               structure) ...
C                yes - set info(6)=0
C                 no - set info(6)=1
C                       and provide the lower (ml) and upper (mu)
C                       bandwidths by setting
C                       iwork(1)=ml
C                       iwork(2)=mu ****
C
C
C        info(7) -- you can specify a maximum (absolute value of)
C              stepsize, so that the code
C              will avoid passing over very
C              large regions.
C
C          ****  do you want the code to decide
C                on its own maximum stepsize?
C                yes - set info(7)=0
C                 no - set info(7)=1
C                      and define hmax by setting
C                      rwork(2)=hmax ****
C
C        info(8) -- differential/algebraic problems
C              may occaisionally suffer from
C              severe scaling difficulties on the
C              first step. if you know a great deal
C              about the scaling of your problem, you can
C              help to alleviate this problem by
C              specifying an initial stepsize ho.
C
C          ****  do you want the code to define
C                its own initial stepsize?
C                yes - set info(8)=0
C                 no - set info(8)=1
C                      and define ho by setting
C                      rwork(3)=ho ****
C
C        info(9) -- if storage is a severe problem,
C              you can save some locations by
C              restricting the maximum order maxord.
C              the default value is 5. for each
C              order decrease below 5, the code
C              requires neq fewer locations, however
C              it is likely to be slower. in any
C              case, you must have 1 .le. maxord .le. 5
C          ****  do you want the maximum order to
C                default to 5?
C                yes - set info(9)=0
C                 no - set info(9)=1
C                      and define maxord by setting
C                      iwork(3)=maxord ****
C
C        info(10) --if you know that the solutions to your equations wil
C               always be nonnegative, it may help to set this
C               parameter.  however, it is probably best to
C               try the code without using this option first,
C               and only to use this option if that doesn't
C               work very well.
C           ****  do you want the code to solve the problem without
C                 invoking any special nonnegativity constraints?
C                  yes - set info(10)=0
C                   no - set info(10)=1
C
C        info(11) --sdassl normally requires the initial t,
C               y, and yprime to be consistent.  that is,
C               you must have g(t,y,yprime) = 0 at the initial
C               time.  if you do not know the initial
C               derivative precisely, you can let sdassl try
C               to compute it.
C          ****   are the initial t, y, yprime consistent?
C                 yes - set info(11) = 0
C                  no - set info(11) = 1,
C                       and set yprime to an initial approximation
C                       to yprime.  (if you have no idea what
C                       yprime should be, set it to zero. note
C                       that the initial y should be such
C                       that there must exist a yprime so that
C                       g(t,y,yprime) = 0.)
C
C   rtol, atol -- you must assign relative (rtol) and absolute (atol
C               error tolerances to tell the code how accurately you wan
C               the solution to be computed. they must be defined as
C               variables because the code may change them. you have two
C               choices --
C                     both rtol and atol are scalars. (info(2)=0)
C                     both rtol and atol are vectors. (info(2)=1)
C               in either case all components must be non-negative.
C
C               the tolerances are used by the code in a local error tes
C               at each step which requires roughly that
C                     abs(local error) .le. rtol*abs(y)+atol
C               for each vector component.
C               (more specifically, a root-mean-square norm is used to
C               measure the size of vectors, and the error test uses the
C               magnitude of the solution at the beginning of the step.)
C
C               the true (global) error is the difference between the tr
C               solution of the initial value problem and the computed
C               approximation. practically all present day codes.
C               including this one, control the local error at each step
C               and do not even attempt to control the global error
C               directly.
C               usually, but not always, the true accuracy of
C               the computed y is comparable to the error tolerances. th
C               code will usually, but not always, deliver a more accura
C               solution if you reduce the tolerances and integrate agai
C               by comparing two such solutions you can get a fairly
C               reliable idea of the true error in the solution at the
C               bigger tolerances.
C
C               setting atol=0. results in a pure relative error test on
C               that component. setting rtol=0. results in a pure absolu
C               error test on that component. a mixed test with non-zero
C               rtol and atol corresponds roughly to a relative error
C               test when the solution component is much bigger than ato
C               and to an absolute error test when the solution componen
C               is smaller than the threshold atol.
C
C               the code will not attempt to compute a solution at an
C               accuracy unreasonable for the machine being used. it wil
C               advise you if you ask for too much accuracy and inform
C               you as to the maximum accuracy it believes possible.
C
C  rwork(*) -- dimension this real work array of length lrw in your
C               calling program.
C
C  lrw -- set it to the declared length of the rwork array.
C               you must have
C                    lrw .ge. 40+(maxord+4)*neq+neq**2
C               for the full (dense) jacobian case (when info(6)=0),  or
C                    lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
C               for the banded user-defined jacobian case
C               (when info(5)=1 and info(6)=1), or
C                     lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
C                           +2*(neq/(ml+mu+1)+1)
C               for the banded finite-difference-generated jacobian case
C               (when info(5)=0 and info(6)=1)
C
C  iwork(*) -- dimension this integer work array of length liw in
C             your calling program.
C
C  liw -- set it to the declared length of the iwork array.
C               you must have liw .ge. 20+neq
C
C  rpar, ipar -- these are parameter arrays, of real and integer
C               type, respectively. you can use them for communication
C               between your program that calls sdassl and the
C               res subroutine (and the jac subroutine). they are not
C               altered by sdassl. if you do not need rpar or ipar, igno
C               these parameters by treating them as dummy arguments. if
C               you do choose to use them, dimension them in your callin
C               program and in res (and in jac) as arrays of appropriate
C               length.
C
C  jac -- if you have set info(5)=0, you can ignore this parameter
C               by treating it as a dummy argument. otherwise, you must
C               provide a subroutine of the form
C               jac(t,y,yprime,pd,cj,rpar,ipar)
C               to define the matrix of partial derivatives
C               pd=dg/dy+cj*dg/dyprime
C               cj is a scalar which is input to jac.
C               for the given values of t,y,yprime, the
C               subroutine must evaluate the non-zero partial
C               derivatives for each equation and each solution
C               compowent, and store these values in the
C               matrix pd. the elements of pd are set to zero
C               before each call to jac so only non-zero elements
C               need to be defined.
C
C               subroutine jac must not alter t,y,(*),yprime(*),or cj.
C               you must declare the name jac in an
C               external statement in your program that calls
C               sdassl. you must dimension y, yprime and pd
C               in jac.
C
C               the way you must store the elements into the pd matrix
C               depends on the structure of the matrix which you
C               indicated by info(6).
C               *** info(6)=0 -- full (dense) matrix ***
C                   when you evaluate the (non-zero) partial derivative
C                   of equation i with respect to variable j, you must
C               store it in pd according to
C                   pd(i,j) = * dg(i)/dy(j)+cj*dg(i)/dyprime(j)*
C               *** info(6)=1 -- banded jacobian with ml lower and mu
C                   upper diagonal bands (refer to info(6) description o
C                   ml and mu) ***
C                   when you evaluate the (non-zero) partial derivative
C                   of equation i with respect to variable j, you must
C                   store it in pd according to
C                   irow = i - j + ml + mu + 1
C                   pd(irow,j) = *dg(i)/dy(j)+cj*dg(i)/dyprime(j)*
C               rpar and ipar are real and integer parameter arrays whic
C               you can use for communication between your calling
C               program and your jacobian subroutine jac. they are not
C               altered by sdassl. if you do not need rpar or ipar, igno
C               these parameters by treating them as dummy arguments. if
C               you do choose to use them, dimension them in your callin
C               program and in jac as arrays of appropriate length.
C
C
C
C  optionally replaceable norm routine:
C  sdassl uses a weighted norm sdanrm to measure the size
C  of vectors such as the estimated error in each step.
C  a function subprogram
C    double precision function sdanrm(neq,v,wt,rpar,ipar)
C    dimension v(neq),wt(neq)
C  is used to define this norm.  here, v is the vector
C  whose norm is to be computed, and wt is a vector of
C  weights.  a sdanrm routine has been included with sdassl
C  which computes the weighted root-mean-square norm
C  given by
C    sdanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
C  this norm is suitable for most problems.  in some
C  special cases, it may be more convenient and/or
C  efficient to define your own norm by writing a function
C  subprogram to be called instead of sdanrm.  this should
C  however, be attempted only after careful thought and
C  consideration.
C
C
C------output-after any return from sdassl----
C
C  the principal aim of the code is to return a computed solution at
C  tout, although it is also possible to obtain intermediate results
C  along the way. to find out whether the code achieved its goal
C  or if the integration process was interrupted before the task was
C  completed, you must check the idid parameter.
C
C
C   t -- the solution was successfully advanced to the
C               output value of t.
C
C   y(*) -- contains the computed solution approximation at t.
C
C   yprime(*) -- contains the computed derivative
C               approximation at t
C
C   idid -- reports what the code did
C
C                     *** task completed ***
C                reported by positive values of idid
C
C           idid = 1 -- a step was successfully taken in the
C                   intermediate-output mode. the code has not
C                   yet reached tout.
C
C           idid = 2 -- the integration to tout was successfully
C                   completed (t=tout) by stepping exactly to tout.
C
C           idid = 3 -- the integration to tout was successfully
C                   completed (t=tout) by stepping past tout.
C                   y(*) is obtained by interpolation.
C                   yprime(*) is obtained by interpolation.
C
C                    *** task interrupted ***
C                reported by negative values of idid
C
C           idid = -1 -- a large amount of work has been expended.
C                   (about 500 steps)
C
C           idid = -2 -- the error tolerances are too stringent.
C
C           idid = -3 -- the local error test cannot be satisfied
C                   because you specified a zero component in atol
C                   and the corresponding computed solution
C                   component is zero. thus, a pure relative error
C                   test is impossible for this component.
C
C           idid = -6 -- sdassl had repeated error test
C                   failures on the last attempted step.
C
C           idid = -7 -- the corrector could not converge.
C
C           idid = -8 -- the matrix of partial derivatives
C                   is singular.
C
C           idid = -9 -- the corrector could not converge.
C                   there were repeated error test failures
C                   in this step.
C
C           idid =-10 -- the corrector could not converge
C                   because ires was equal to minus one.
C
C           idid =-11 -- ires equal to -2 was encountered
C                   and control is being returned to the
C                   calling program.
C
C           idid =-12 -- sdassl failed to compute the initial
C                   yprime.
C
C
C
C           idid = -13,..,-32 -- not applicable for this code
C
C                    *** task terminated ***
C                reported by the value of idid=-33
C
C           idid = -33 -- the code has encountered trouble from which
C                   it cannot recover. a message is printed
C                   explaining the trouble and control is returned
C                   to the calling program. for example, this occurs
C                   when invalid input is detected.
C
C   rtol, atol -- these quantities remain unchanged except when
C               idid = -2. in this case, the error tolerances have been
C               increased by the code to values which are estimated to b
C               appropriate for continuing the integration. however, the
C               reported solution at t was obtained using the input valu
C               of rtol and atol.
C
C   rwork, iwork -- contain information which is usually of no
C               interest to the user but necessary for subsequent calls.
C               however, you may find use for
C
C               rwork(3)--which contains the step size h to be
C                       attempted on the next step.
C
C               rwork(4)--which contains the current value of the
C                       independent variable, i.e. the farthest point
C                       integration has reached. this will be different
C                       from t only when interpolation has been
C                       performed (idid=3).
C
C               rwork(7)--which contains the stepsize used
C                       on the last successful step.
C
C               iwork(7)--which contains the order of the method to
C                       be attempted on the next step.
C
C               iwork(8)--which contains the order of the method used
C                       on the last step.
C
C               iwork(11)--which contains the number of steps taken so f
C
C               iwork(12)--which contains the number of calls to res
C                        so far.
C
C               iwork(13)--which contains the number of evaluations of
C                        the matrix of partial derivatives needed so far
C
C               iwork(14)--which contains the total number
C                        of error test failures so far.
C
C               iwork(15)--which contains the total number
C                        of convergence test failures so far.
C                        (includes singular iteration matrix
C                        failures.)
C
C
C
C   input -- what to do to continue the integration
C            (calls after the first)                **
C
C     this code is organized so that subsequent calls to continue the
C     integration involve little (if any) additional effort on your
C     part. you must monitor the idid parameter in order to determine
C     what to do next.
C
C     recalling that the principal task of the code is to integrate
C     from t to tout (the interval mode), usually all you will need
C     to do is specify a new tout upon reaching the current tout.
C
C     do not alter any quantity not specifically permitted below,
C     in particular do not alter neq,t,y(*),yprime(*),rwork(*),iwork(*)
C     or the differential equation in subroutine res. any such
C     alteration constitutes a new problem and must be treated as such,
C     i.e. you must start afresh.
C
C     you cannot change from vector to scalar error control or vice
C     versa (info(2)) but you can change the size of the entries of
C     rtol, atol. increasing a tolerance makes the equation easier
C     to integrate. decreasing a tolerance will make the equation
C     harder to integrate and should generally be avoided.
C
C     you can switch from the intermediate-output mode to the
C     interval mode (info(3)) or vice versa at any time.
C
C     if it has been necessary to prevent the integration from going
C     past a point tstop (info(4), rwork(1)), keep in mind that the
C     code will not integrate to any tout beyound the currently
C     specified tstop. once tstop has been reached you must change
C     the value of tstop or set info(4)=0. you may change info(4)
C     or tstop at any time but you must supply the value of tstop in
C     rwork(1) whenever you set info(4)=1.
C
C     do not change info(5), info(6), iwork(1), or iwork(2)
C     unless you are going to restart the code.
C
C                    *** following a completed task ***
C     if
C     idid = 1, call the code again to continue the integration
C                  another step in the direction of tout.
C
C     idid = 2 or 3, define a new tout and call the code again.
C                  tout must be different from t. you cannot change
C                  the direction of integration without restarting.
C
C                    *** following an interrupted task ***
C                  to show the code that you realize the task was
C                  interrupted and that you want to continue, you
C                  must take appropriate action and set info(1) = 1
C     if
C     idid = -1, the code has taken about 500 steps.
C                  if you want to continue, set info(1) = 1 and
C                  call the code again. an additional 500 steps
C                  will be allowed.
C
C
C     idid = -2, the error tolerances rtol, atol have been
C                  increased to values the code estimates appropriate
C                  for continuing. you may want to change them
C                  yourself. if you are sure you want to continue
C                  with relaxed error tolerances, set info(1)=1 and
C                  call the code again.
C
C     idid = -3, a solution component is zero and you set the
C                  corresponding component of atol to zero. if you
C                  are sure you want to continue, you must first
C                  alter the error criterion to use positive values
C                  for those components of atol corresponding to zero
C                  solution components, then set info(1)=1 and call
C                  the code again.
C
C     idid = -4,-5  --- cannot occur with this code
C
C     idid = -6, repeated error test failures occurred on the
C                  last attempted step in sdassl. a singularity in the
C                  solution may be present. if you are absolutely
C                  certain you want to continue, you should restart
C                  the integration.(provide initial values of y and
C                  yprime which are consistent)
C
C     idid = -7, repeated convergence test failures occurred
C                  on the last attempted step in sdassl. an inaccurate o
C                  illconditioned jacobian may be the problem. if you
C                  are absolutely certain you want to continue, you
C                  should restart the integration.
C
C     idid = -8, the matrix of partial derivatives is singular.
C                  some of your equations may be redundant.
C                  sdassl cannot solve the problem as stated.
C                  it is possible that the redundant equations
C                  could be removed, and then sdassl could
C                  solve the problem. it is also possible
C                  that a solution to your problem either
C                  does not exist or is not unique.
C
C     idid = -9, sdassl had multiple convergence test
C                  failures, preceeded by multiple error
C                  test failures, on the last attempted step.
C                  it is possible that your problem
C                  is ill-posed, and cannot be solved
C                  using this code.  or, there may be a
C                  discontinuity or a singularity in the
C                  solution.  if you are absolutely certain
C                  you want to continue, you should restart
C                  the integration.
C
C    idid =-10, sdassl had multiple convergence test failures
C                  because ires was equal to minus one.
C                  if you are absolutely certain you want
C                  to continue, you should restart the
C                  integration.
C
C    idid =-11, ires=-2 was encountered, and control is being
C                  returned to the calling program.
C
C    idid =-12, sdassl failed to compute the initial yprime.
C               this could happen because the initial
C               approximation to yprime was not very good, or
C               if a yprime consistent with the initial y
C               does not exist.  the problem could also be caused
C               by an inaccurate or singular iteration matrix.
C
C
C
C     idid = -13,..,-32 --- cannot occur with this code
C
C                       *** following a terminated task ***
C     if idid= -33, you cannot continue the solution of this
C                  problem. an attempt to do so will result in your
C                  run being terminated.
C
C  ---------------------------------------------------------------------
C
C***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
C                  SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
C                  SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
C***ROUTINES CALLED  SDASTP,SDAINI,SDANRM,SDAWTS,SDATRP,XERRWV,R1MACH
C***COMMON BLOCKS    SDA001
C
C   000330  Modified code to use PARAMETER statement to define pointers
C           in IWORK that are not in COMMON.
C
C***END PROLOGUE SDASSL
C
C
      LOGICAL DONE
      EXTERNAL RES,JAC
      DIMENSION Y(*),YPRIME(*)
      DIMENSION INFO(15)
      DIMENSION RWORK(*),IWORK(*)
      DIMENSION RTOL(*),ATOL(*)
      DIMENSION RPAR(*),IPAR(*)
C
C     SET POINTERS INTO IWORK   000330
C
      PARAMETER (LNPD=16, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *  LNS=9, LNSTL=10, LIWM=1)

      COMMON/SDA001/NPD,NTEMP,
     *   LML,LMU,LMXORD,LMTYPE,
     *   LNST,LNRE,LNJE,LETF,LCTF,LIPVT
      DATA LTSTOP,LHMAX,LH,LTN,
     *   LCJ,LCJOLD,LHOLD,LS,LROUND,
     *   LALPHA,LBETA,LGAMMA,
     *   LPSI,LSIGMA,LDELTA
     *   /1,2,3,4,
     *   5,6,7,8,9,
     *   11,17,23,
     *   29,35,41/
      IF(INFO(1).NE.0)GO TO 100
C
C-----------------------------------------------------------------------
C     this block is executed for the initial call only.
C     it contains checking of inputs and initializations.
C-----------------------------------------------------------------------
C
C     first check info array to make sure all elements of info
C     are either zero or one.
      DO 10 I=2,11
         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
10       CONTINUE
C
      IF(NEQ.LE.0)GO TO 702
C
C     set pointers into iwork
      LML=1
      LMU=2
      LMXORD=3
      LMTYPE=4
      LNST=11
      LNRE=12
      LNJE=13
      LETF=14
      LCTF=15
      LIPVT=21
C
C     check and compute maximum order
      MXORD=5
      IF(INFO(9).EQ.0)GO TO 20
         MXORD=IWORK(LMXORD)
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
20       IWORK(LMXORD)=MXORD
C
C     compute mtype,lenpd,lenrw.check ml and mu.
      IF(INFO(6).NE.0)GO TO 40
         LENPD=NEQ**2
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
         IF(INFO(5).NE.0)GO TO 30
            IWORK(LMTYPE)=2
            GO TO 60
30          IWORK(LMTYPE)=1
            GO TO 60
40    IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
      IF(INFO(5).NE.0)GO TO 50
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NEQ/MBAND)+1
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
         GO TO 60
50       IWORK(LMTYPE)=4
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
C
C     check lengths of rwork and iwork
60    LENIW=20+NEQ
      IWORK(LNPD)=LENPD
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     check to see that tout is different from t
      IF(TOUT .EQ. T)GO TO 719
C
C     check hmax
      IF(INFO(7).EQ.0)GO TO 70
         HMAX=RWORK(LHMAX)
         IF(HMAX.LE.0.0E0)GO TO 710
70    CONTINUE
C
C     initialize counters
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
C
      IWORK(LNSTL)=0
      IDID=1
      GO TO 200
C
C-----------------------------------------------------------------------
C     this block is for continuation calls
C     only. here we check info(1),and if the
C     last step was interrupted we check whether
C     appropriate action was taken.
C-----------------------------------------------------------------------
C
100   CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      IF(INFO(1).NE.-1)GO TO 701
C     if we are here, the last step was interrupted
C     by an error condition from sdastp,and
C     appropriate action was not taken. this
C     is a fatal error.
      CALL XERRWV(
     *'DASSL--  THE LAST STEP TERMINATED WITH A NEGATIVE',
     *49,201,0,0,0,0,0,0.0E0,0.0E0)
      CALL XERRWV(
     *'DASSL--  VALUE (=I1) OF IDID AND NO APPROPRIATE',
     *47,202,0,1,IDID,0,0,0.0E0,0.0E0)
      CALL XERRWV(
     *'DASSL--  ACTION WAS TAKEN. RUN TERMINATED',
     *41,203,1,0,0,0,0,0.0E0,0.0E0)
      RETURN
110   CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
C
C-----------------------------------------------------------------------
C     this block is executed on all calls.
C     the error tolerance parameters are
C     checked, and the work array pointers
C     are set.
C-----------------------------------------------------------------------
C
200   CONTINUE
C     check rtol,atol
      NZFLG=0
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 210 I=1,NEQ
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
         IF(RTOLI.GT.0.0E0.OR.ATOLI.GT.0.0E0)NZFLG=1
         IF(RTOLI.LT.0.0E0)GO TO 706
         IF(ATOLI.LT.0.0E0)GO TO 707
210      CONTINUE
      IF(NZFLG.EQ.0)GO TO 708
C
C     set up rwork storage.iwork storage is fixed
C     in data statement.
      LE=LDELTA+NEQ
      LWT=LE+NEQ
      LPHI=LWT+NEQ
      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
      LWM=LPD
      NPD=1
      NTEMP=NPD+IWORK(LNPD)
      IF(INFO(1).EQ.1)GO TO 400
C
C-----------------------------------------------------------------------
C     this block is executed on the initial call
C     only. set the initial step size, and
C     the error weight vector, and phi.
C     compute initial yprime, if necessary.
C-----------------------------------------------------------------------
C
300   CONTINUE
      TN=T
      IDID=1
C
C     set error weight vector wt
      CALL SDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      DO 305 I = 1,NEQ
         IF(RWORK(LWT+I-1).LE.0.0E0) GO TO 713
305      CONTINUE
C
C     compute unit roundoff and hmin
      UROUND = R1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0E0*UROUND*AMAX1(ABS(T),ABS(TOUT))
C
C     check initial interval to see that it is long enough
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     check ho, if this was input
      IF (INFO(8) .EQ. 0) GO TO 310
         HO = RWORK(LH)
         IF ((TOUT - T)*HO .LT. 0.0E0) GO TO 711
         IF (HO .EQ. 0.0E0) GO TO 712
         GO TO 320
310    CONTINUE
C
C     compute initial stepsize, to be used by either
C     sdastp or sdaini, depending on info(11)
      HO = 0.001E0*TDIST
      YPNORM = SDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5E0/HO) HO = 0.5E0/YPNORM
      HO = SIGN(HO,TOUT-T)
C     adjust ho if necessary to meet hmax bound
320   IF (INFO(7) .EQ. 0) GO TO 330
         RH = ABS(HO)/HMAX
         IF (RH .GT. 1.0E0) HO = HO/RH
C     compute tstop, if applicable
330   IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0E0) GO TO 715
         IF ((T + HO - TSTOP)*HO .GT. 0.0E0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0E0) GO TO 709
C
C     compute initial derivative, if applicable
340   IF (INFO(11) .EQ. 0) GO TO 350
      CALL SDAINI(T,Y,YPRIME,NEQ,
     *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,
     *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),INFO(10))
      IF (IDID .LT. 0) GO TO 390
C
C     load h with ho.  store h in rwork(lh)
350   H = HO
      RWORK(LH) = H
C
C     load y and h*yprime into phi(*,1) and phi(*,2)
360   ITEMP = LPHI + NEQ
      DO 370 I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
C
390   GO TO 500
C
C-------------------------------------------------------
C     this block is for continuation calls only. its
C     purpose is to check stop conditions before
C     taking a step.
C     adjust h if necessary to meet hmax bound
C-------------------------------------------------------
C
400   CONTINUE
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = ABS(H)/HMAX
         IF(RH .GT. 1.0E0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0E0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0E0)GO TO 490
      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0E0) GO TO 490
      IF((TN - TOUT)*H .GT. 0.0E0) GO TO 425
      CALL SDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
425   CONTINUE
      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
430   IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
      IF((TN-TSTOP)*H.GT.0.0E0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0E0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0E0)GO TO 450
      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GT. 0.0E0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0E0) GO TO 709
      IF((TN-T)*H .LE. 0.0E0) GO TO 450
      IF((TN - TOUT)*H .GT. 0.0E0) GO TO 445
      CALL SDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
445   CONTINUE
      CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
450   CONTINUE
C     check whether we are with in roundoff of tstop
      IF(ABS(TN-TSTOP).GT.100.0E0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 460
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H*(1.0E0+4.0E0*UROUND)
      IF((TNEXT-TSTOP)*H.LE.0.0E0)GO TO 490
      H=(TSTOP-TN)*(1.0E0-4.0E0*UROUND)
      RWORK(LH)=H
C
490   IF (DONE) GO TO 590
C
C-------------------------------------------------------
C     the next block contains the call to the
C     one-step integrator sdastp.
C     this is a looping point for the integration
C     steps.
C     check for too many steps.
C     update wt.
C     check for too much accuracy requested.
C     compute minimum stepsize.
C-------------------------------------------------------
C
500   CONTINUE
C     check for failure to compute initial yprime
      IF (IDID .EQ. -12) GO TO 527
C
C     check for too many steps
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)
     *   GO TO 510
           IDID=-1
           GO TO 527
C
C     update wt
510   CALL SDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),
     *  RWORK(LWT),RPAR,IPAR)
      DO 520 I=1,NEQ
         IF(RWORK(I+LWT-1).GT.0.0E0)GO TO 520
           IDID=-3
           GO TO 527
520   CONTINUE
C
C     test for too much accuracy requested.
      R=SDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*
     *   100.0E0*UROUND
      IF(R.LE.1.0E0)GO TO 525
C     multiply rtol and atol by r and return
      IF(INFO(2).EQ.1)GO TO 523
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
           IDID=-2
           GO TO 527
523   DO 524 I=1,NEQ
           RTOL(I)=R*RTOL(I)
524        ATOL(I)=R*ATOL(I)
      IDID=-2
      GO TO 527
525   CONTINUE
C
C     compute minimum stepsize
      HMIN=4.0E0*UROUND*AMAX1(ABS(TN),ABS(TOUT))
C
      CALL SDASTP(TN,Y,YPRIME,NEQ,
     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,
     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *   RWORK(LWM),IWORK(LIWM),
     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *   RWORK(LPSI),RWORK(LSIGMA),
     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),
     *   RWORK(LS),HMIN,RWORK(LROUND),
     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),
     *   IWORK(LKOLD),IWORK(LNS),INFO(10))
527   IF(IDID.LT.0)GO TO 600
C
C------------------------------------------------------
C     this block handles the case of a successful
C     return from sdastp (idid=1) test for
C     stop conditions.
C--------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0E0)GO TO 500
             CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
530          IF((TN-TOUT)*H.GE.0.0E0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
535          CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
540   IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0E0)GO TO 542
         CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
542   IF(ABS(TN-TSTOP).LE.100.0E0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H*(1.0E0+4.0E0*UROUND)
      IF((TNEXT-TSTOP)*H.LE.0.0E0)GO TO 500
      H=(TSTOP-TN)*(1.0E0-4.0E0*UROUND)
      GO TO 500
545   IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GE.0.0E0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0E0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   IDID=2
      T=TSTOP
      GO TO 580
555   CALL SDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
580   CONTINUE
C
C--------------------------------------------------------
C     all successful returns from sdassl are made from
C     this block.
C--------------------------------------------------------
C
590   CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     this block handles all unsuccessful
C     returns other than for illegal input.
C-----------------------------------------------------------------------
C
600   CONTINUE
      ITEMP=-IDID
      GO TO (610,620,630,690,690,640,650,660,670,675,
     *  680,685), ITEMP
C
C     the maximum number of steps was taken before
C     reaching tout
610   CALL XERRWV(
     *'DASSL--  AT CURRENT T (=R1)  500 STEPS',
     *38,610,0,0,0,0,1,TN,0.0E0)
      CALL XERRWV('DASSL--  TAKEN ON THIS CALL BEFORE REACHING TOUT',
     *48,611,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     too much accuracy for machine precision
620   CALL XERRWV(
     *'DASSL--  AT T (=R1) TOO MUCH ACCURACY REQUESTED',
     *47,620,0,0,0,0,1,TN,0.0E0)
      CALL XERRWV(
     *'DASSL--  FOR PRECISION OF MACHINE. RTOL AND ATOL',
     *48,621,0,0,0,0,0,0.0E0,0.0E0)
      CALL XERRWV(
     *'DASSL--  WERE INCREASED TO APPROPRIATE VALUES',
     *45,622,0,0,0,0,0,0.0E0,0.0E0)
C
      GO TO 690
C     wt(i) .le. 0.0E0 for some i (not at start of problem)
630   CALL XERRWV(
     *'DASSL--  AT T (=R1) SOME ELEMENT OF WT',
     *38,630,0,0,0,0,1,TN,0.0E0)
      CALL XERRWV('DASSL--  HAS BECOME .LE. 0.0',
     *28,631,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     error test failed repeatedly or with h=hmin
640   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *44,640,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',
     *57,641,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     corrector convergence failed repeatedly or with h=hmin
650   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *44,650,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  CORRECTOR FAILED TO CONVERGE REPEATEDLY',
     *48,651,0,0,0,0,0,0.0E0,0.0E0)
      CALL XERRWV(
     *'DASSL--  OR WITH ABS(H)=HMIN',
     *28,652,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     the iteration matrix is singular
660   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *44,660,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  ITERATION MATRIX IS SINGULAR',
     *37,661,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     corrector failure preceeded by error test failures.
670   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *44,670,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  CORRECTOR COULD NOT CONVERGE.  ALSO, THE',
     *49,671,0,0,0,0,0,0.0E0,0.0E0)
      CALL XERRWV(
     *'DASSL--  ERROR TEST FAILED REPEATEDLY.',
     *38,672,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     corrector failure because ires = -1
675   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *44,675,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  CORRECTOR COULD NOT CONVERGE BECAUSE',
     *45,676,0,0,0,0,0,0.0E0,0.0E0)
      CALL XERRWV(
     *'DASSL--  IRES WAS EQUAL TO MINUS ONE',
     *36,677,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     failure because ires = -2
680   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2)',
     *40,680,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  IRES WAS EQUAL TO MINUS TWO',
     *36,681,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
C
C     failed to compute initial yprime
685   CALL XERRWV(
     *'DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *44,685,0,0,0,0,2,TN,HO)
      CALL XERRWV(
     *'DASSL--  INITIAL YPRIME COULD NOT BE COMPUTED',
     *45,686,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 690
690   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C-----------------------------------------------------------------------
C     this block handles all error returns due
C     to illegal input, as detected before calling
C     sdastp. first the error message routine is
C     called. if this happens twice in
C     succession, execution is terminated
C
C-----------------------------------------------------------------------
701   CALL XERRWV(
     *'DASSL--  SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE',
     *55,1,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 750
702   CALL XERRWV('DASSL--  NEQ (=I1) .LE. 0',
     *25,2,0,1,NEQ,0,0,0.0E0,0.0E0)
      GO TO 750
703   CALL XERRWV('DASSL--  MAXORD (=I1) NOT IN RANGE',
     *34,3,0,1,MXORD,0,0,0.0E0,0.0E0)
      GO TO 750
704   CALL XERRWV(
     *'DASSL--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)',
     *60,4,0,2,LENRW,LRW,0,0.0E0,0.0E0)
      GO TO 750
705   CALL XERRWV(
     *'DASSL--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)',
     *60,5,0,2,LENIW,LIW,0,0.0E0,0.0E0)
      GO TO 750
706   CALL XERRWV(
     *'DASSL--  SOME ELEMENT OF RTOL IS .LT. 0',
     *39,6,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 750
707   CALL XERRWV(
     *'DASSL--  SOME ELEMENT OF ATOL IS .LT. 0',
     *39,7,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 750
708   CALL XERRWV(
     *'DASSL--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO',
     *47,8,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 750
709   CALL XERRWV(
     *'DASSL--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)',
     *54,9,0,0,0,0,2,TSTOP,TOUT)
      GO TO 750
710   CALL XERRWV('DASSL--  HMAX (=R1) .LT. 0.0',
     *28,10,0,0,0,0,1,HMAX,0.0E0)
      GO TO 750
711   CALL XERRWV('DASSL--  TOUT (=R1) BEHIND T (=R2)',
     *34,11,0,0,0,0,2,TOUT,T)
      GO TO 750
712   CALL XERRWV('DASSL--  INFO(8)=1 AND H0=0.0',
     *29,12,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 750
713   CALL XERRWV('DASSL--  SOME ELEMENT OF WT IS .LE. 0.0',
     *39,13,0,0,0,0,0,0.0E0,0.0E0)
      GO TO 750
714   CALL XERRWV(
     *'DASSL--  TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION',
     *61,14,0,0,0,0,2,TOUT,T)
      GO TO 750
715   CALL XERRWV(
     *'DASSL--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)',
     *49,15,0,0,0,0,2,TSTOP,T)
      GO TO 750
717   CALL XERRWV(
     *'DASSL--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ',
     *52,17,0,1,IWORK(LML),0,0,0.0E0,0.0E0)
      GO TO 750
718   CALL XERRWV(
     *'DASSL--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ',
     *52,18,0,1,IWORK(LMU),0,0,0.0E0,0.0E0)
      GO TO 750
719   CALL XERRWV(
     *'DASSL--  TOUT (=R1) IS EQUAL TO T (=R2)',
     *39,19,0,0,0,0,2,TOUT,T)
      GO TO 750
750   IF(INFO(1).EQ.-1) GO TO 760
      INFO(1)=-1
      IDID=-33
      RETURN
760   CALL XERRWV(
     *'DASSL--  REPEATED OCCURRENCES OF ILLEGAL INPUT',
     *46,801,0,0,0,0,0,0.0E0,0.0E0)
770   CALL XERRWV(
     *'DASSL--  RUN TERMINATED. APPARENT INFINITE LOOP',
     *47,802,1,0,0,0,0,0.0E0,0.0E0)
      RETURN
C-----------end of subroutine sdassl------------------------------------
      END
