      SUBROUTINE QAWFE(F,A,OMEGA,INTEGR,EPSABS,LIMLST,LIMIT,MAXP1,
     1   RESULT,ABSERR,NEVAL,IER,RSLST,ERLST,IERLST,LST,ALIST,BLIST,
     2   RLIST,ELIST,IORD,NNLOG,CHEBMO)
C***BEGIN PROLOGUE  QAWFE
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A3A1
C***KEYWORDS  AUTOMATIC INTEGRATOR,CONVERGENCE ACCELERATION WITH DQELG,
C             FOURIER INTEGRALS,INTEGRATION BETWEEN ZEROS WITH DQAWOE,
C             SPECIAL-PURPOSE
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a
C            given Fourier integal
C            I = Integral of F(X)*W(X) over (A,INFINITY)
C             where W(X) = COS(OMEGA*X) or W(X) = SIN(OMEGA*X),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.EPSABS.
C***DESCRIPTION
C
C        Computation of Fourier integrals
C        Standard fortran subroutine
C        Real version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Real
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to
C                     be declared E X T E R N A L in the driver program.
C
C            A      - Real
C                     Lower limit of integration
C
C            OMEGA  - Real
C                     Parameter in the WEIGHT function
C
C            INTEGR - Integer
C                     Indicates which WEIGHT function is used
C                     INTEGR = 1      W(X) = COS(OMEGA*X)
C                     INTEGR = 2      W(X) = SIN(OMEGA*X)
C                     If INTEGR.NE.1.AND.INTEGR.NE.2, the routine will
C                     end with IER = 6.
C
C            EPSABS - Real
C                     absolute accuracy requested, EPSABS.GT.0
C                     If EPSABS.LE.0, the routine will end with IER = 6.
C
C            LIMLST - Integer
C                     LIMLST gives an upper bound on the number of
C                     cycles, LIMLST.GE.1.
C                     If LIMLST.LT.3, the routine will end with IER = 6.
C
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     allowed in the partition of each cycle, LIMIT.GE.1
C                     each cycle, LIMIT.GE.1.
C
C            MAXP1  - Integer
C                     Gives an upper bound on the number of
C                     Chebyshev moments which can be stored, I.E.
C                     for the intervals of lengths ABS(B-A)*2**(-L),
C                     L=0,1, ..., MAXP1-2, MAXP1.GE.1
C
C         ON RETURN
C            RESULT - Real
C                     Approximation to the integral X
C
C            ABSERR - Real
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - IER = 0 Normal and reliable termination of
C                             the routine. It is assumed that the
C                             requested accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine. The
C                             estimates for integral and error are less
C                             reliable. It is assumed that the requested
C                             accuracy has not been achieved.
C            ERROR MESSAGES
C                    If OMEGA.NE.0
C                     IER = 1 Maximum number of  cycles  allowed
C                             Has been achieved., i.e. of subintervals
C                             (A+(K-1)C,A+KC) where
C                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
C                             for K = 1, 2, ..., LST.
C                             One can allow more cycles by increasing
C                             the value of LIMLST (and taking the
C                             according dimension adjustments into
C                             account).
C                             Examine the array IWORK which contains
C                             the error flags on the cycles, in order to
C                             look for eventual local integration
C                             difficulties. If the position of a local
C                             difficulty can be determined (e.g.
C                             SINGULARITY, DISCONTINUITY within the
C                             interval) one will probably gain from
C                             splitting up the interval at this point
C                             and calling appropriate integrators on
C                             the subranges.
C                         = 4 The extrapolation table constructed for
C                             convergence acceleration of the series
C                             formed by the integral contributions over
C                             the cycles, does not converge to within
C                             the requested accuracy. As in the case of
C                             IER = 1, it is advised to examine the
C                             array IWORK which contains the error
C                             flags on the cycles.
C                         = 6 The input is invalid because
C                             (INTEGR.NE.1 AND INTEGR.NE.2) or
C                              EPSABS.LE.0 or LIMLST.LT.3.
C                              RESULT, ABSERR, NEVAL, LST are set
C                              to zero.
C                         = 7 Bad integrand behaviour occurs within one
C                             or more of the cycles. Location and type
C                             of the difficulty involved can be
C                             determined from the vector IERLST. Here
C                             LST is the number of cycles actually
C                             needed (see below).
C                             IERLST(K) = 1 The maximum number of
C                                           subdivisions (= LIMIT) has
C                                           been achieved on the K th
C                                           cycle.
C                                       = 2 Occurrence of roundoff error
C                                           is detected and prevents the
C                                           tolerance imposed on the
C                                           K th cycle, from being
C                                           achieved.
C                                       = 3 Extremely bad integrand
C                                           behaviour occurs at some
C                                           points of the K th cycle.
C                                       = 4 The integration procedure
C                                           over the K th cycle does
C                                           not converge (to within the
C                                           required accuracy) due to
C                                           roundoff in the
C                                           extrapolation procedure
C                                           invoked on this cycle. It
C                                           is assumed that the result
C                                           on this interval is the
C                                           best which can be obtained.
C                                       = 5 The integral over the K th
C                                           cycle is probably divergent
C                                           or slowly convergent. It
C                                           must be noted that
C                                           divergence can occur with
C                                           any other value of
C                                           IERLST(K).
C                    If OMEGA = 0 and INTEGR = 1,
C                    The integral is calculated by means of DQAGIE
C                    and IER = IERLST(1) (with meaning as described
C                    for IERLST(K), K = 1).
C
C            RSLST  - Real
C                     Vector of dimension at least LIMLST
C                     RSLST(K) contains the integral contribution
C                     over the interval (A+(K-1)C,A+KC) where
C                     C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
C                     K = 1, 2, ..., LST.
C                     Note that, if OMEGA = 0, RSLST(1) contains
C                     the value of the integral over (A,INFINITY).
C
C            ERLST  - Real
C                     Vector of dimension at least LIMLST
C                     ERLST(K) contains the error estimate corresponding
C                     with RSLST(K).
C
C            IERLST - Integer
C                     Vector of dimension at least LIMLST
C                     IERLST(K) contains the error flag corresponding
C                     with RSLST(K). For the meaning of the local error
C                     flags see description of output parameter IER.
C
C            LST    - Integer
C                     Number of subintervals needed for the integration
C                     If OMEGA = 0 then LST is set to 1.
C
C            ALIST, BLIST, RLIST, ELIST - Real
C                     vector of dimension at least LIMIT,
C
C            IORD, NNLOG - Integer
C                     Vector of dimension at least LIMIT, providing
C                     space for the quantities needed in the subdivision
C                     process of each cycle
C
C            CHEBMO - Real
C                     Array of dimension at least (MAXP1,25), providing
C                     space for the Chebyshev moments needed within the
C                     cycles
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QAGIE,QAWOE,QELG,R1MACH
C***END PROLOGUE  QAWFE
C
      REAL A,ABSEPS,ABSERR,ALIST,BLIST,CHEBMO,CORREC,CYCLE,
     1  C1,C2,DL,DLA,DRL,ELIST,EP,EPS,EPSA,EPSABS,ERLST,
     2  ERRSUM,FACT,OMEGA,P,PI,P1,PSUM,RESEPS,RESULT,RES3LA,RLIST,RSLST
     3  ,R1MACH,UFLOW
      INTEGER IER,IERLST,INTEGR,IORD,KTMIN,L,LST,LIMIT,LL,MAXP1,
     1    NEV,NEVAL,NNLOG,NRES,NUMRL2
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),CHEBMO(MAXP1,25),ELIST(LIMIT),
     1  ERLST(LIMLST),IERLST(LIMLST),IORD(LIMIT),NNLOG(LIMIT),PSUM(52),
     2  RES3LA(3),RLIST(LIMIT),RSLST(LIMLST)
C
      EXTERNAL F
C
C
C            THE DIMENSION OF  PSUM  IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE QELG (PSUM MUST BE
C            OF DIMENSION (LIMEXP+2) AT LEAST).
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           C1, C2    - END POINTS OF SUBINTERVAL (OF LENGTH
C                       CYCLE)
C           CYCLE     - (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA)
C           PSUM      - VECTOR OF DIMENSION AT LEAST (LIMEXP+2)
C                       (SEE ROUTINE QELG)
C                       PSUM CONTAINS THE PART OF THE EPSILON
C                       TABLE WHICH IS STILL NEEDED FOR FURTHER
C                       COMPUTATIONS.
C                       EACH ELEMENT OF PSUM IS A PARTIAL SUM OF
C                       THE SERIES WHICH SHOULD SUM TO THE VALUE OF
C                       THE INTEGRAL.
C           ERRSUM    - SUM OF ERROR ESTIMATES OVER THE
C                       SUBINTERVALS, CALCULATED CUMULATIVELY
C           EPSA      - ABSOLUTE TOLERANCE REQUESTED OVER CURRENT
C                       SUBINTERVAL
C           CHEBMO    - ARRAY CONTAINING THE MODIFIED CHEBYSHEV
C                       MOMENTS (SEE ALSO ROUTINE QC25F)
C
      DATA P/0.9E+00/,PI/0.31415926535897932E+01/
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
C***FIRST EXECUTABLE STATEMENT  QAWFE
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      NEVAL = 0
      LST = 0
      IER = 0
      IF((INTEGR.NE.1.AND.INTEGR.NE.2).OR.EPSABS.LE.0.0E+00.OR.
     1  LIMLST.LT.3) IER = 6
      IF(IER.EQ.6) GO TO 999
      IF(OMEGA.NE.0.0E+00) GO TO 10
C
C           INTEGRATION BY QAGIE IF OMEGA IS ZERO
C           --------------------------------------
C
      IF(INTEGR.EQ.1) CALL QAGIE(F,0.0E+00,1,EPSABS,0.0E+00,LIMIT,
     1  RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
      RSLST(1) = RESULT
      ERLST(1) = ABSERR
      IERLST(1) = IER
      LST = 1
      GO TO 999
C
C           INITIALIZATIONS
C           ---------------
C
   10 L = ABS(OMEGA)
      DL = 2*L+1
      CYCLE = DL*PI/ABS(OMEGA)
      IER = 0
      KTMIN = 0
      NEVAL = 0
      NUMRL2 = 0
      NRES = 0
      C1 = A
      C2 = CYCLE+A
      P1 = 0.1E+01-P
      EPS = EPSABS
      UFLOW = R1MACH(1)
      IF(EPSABS.GT.UFLOW/P1) EPS = EPSABS*P1
      EP = EPS
      FACT = 0.1E+01
      CORREC = 0.0E+00
      ABSERR = 0.0E+00
      ERRSUM = 0.0E+00
C
C           MAIN DO-LOOP
C           ------------
C
      DO 50 LST = 1,LIMLST
C
C           INTEGRATE OVER CURRENT SUBINTERVAL.
C
        DLA = LST
        EPSA = EPS*FACT
        CALL QAWOE(F,C1,C2,OMEGA,INTEGR,EPSA,0.0E+00,LIMIT,LST,MAXP1,
     1  RSLST(LST),ERLST(LST),NEV,IERLST(LST),LAST,ALIST,BLIST,RLIST,
     2  ELIST,IORD,NNLOG,MOMCOM,CHEBMO)
        NEVAL = NEVAL+NEV
        FACT = FACT*P
        ERRSUM = ERRSUM+ERLST(LST)
        DRL = 0.5E+02*ABS(RSLST(LST))
C
C           TEST ON ACCURACY WITH PARTIAL SUM
C
        IF(ERRSUM+DRL.LE.EPSABS.AND.LST.GE.6) GO TO 80
        CORREC = AMAX1(CORREC,ERLST(LST))
        IF(IERLST(LST).NE.0) EPS = AMAX1(EP,CORREC*P1)
        IF(IERLST(LST).NE.0) IER = 7
        IF(IER.EQ.7.AND.(ERRSUM+DRL).LE.CORREC*0.1E+02.AND.
     1  LST.GT.5) GO TO 80
        NUMRL2 = NUMRL2+1
        IF(LST.GT.1) GO TO 20
        PSUM(1) = RSLST(1)
        GO TO 40
   20   PSUM(NUMRL2) = PSUM(LL)+RSLST(LST)
        IF(LST.EQ.2) GO TO 40
C
C           TEST ON MAXIMUM NUMBER OF SUBINTERVALS
C
        IF(LST.EQ.LIMLST) IER = 1
C
C           PERFORM NEW EXTRAPOLATION
C
        CALL QELG(NUMRL2,PSUM,RESEPS,ABSEPS,RES3LA,NRES)
C
C           TEST WHETHER EXTRAPOLATED RESULT IS INFLUENCED BY
C           ROUNDOFF
C
        KTMIN = KTMIN+1
        IF(KTMIN.GE.15.AND.ABSERR.LE.0.1E-02*(ERRSUM+DRL)) IER = 4
        IF(ABSEPS.GT.ABSERR.AND.LST.NE.3) GO TO 30
        ABSERR = ABSEPS
        RESULT = RESEPS
        KTMIN = 0
C
C           IF IER IS NOT 0, CHECK WHETHER DIRECT RESULT (PARTIAL
C           SUM) OR EXTRAPOLATED RESULT YIELDS THE BEST INTEGRAL
C           APPROXIMATION
C
        IF((ABSERR+0.1E+02*CORREC).LE.EPSABS.OR.
     1  (ABSERR.LE.EPSABS.AND.0.1E+02*CORREC.GE.EPSABS)) GO TO 60
   30   IF(IER.NE.0.AND.IER.NE.7) GO TO 60
   40   LL = NUMRL2
        C1 = C2
        C2 = C2+CYCLE
   50 CONTINUE
C
C         SET FINAL RESULT AND ERROR ESTIMATE
C         -----------------------------------
C
   60 ABSERR = ABSERR+0.1E+02*CORREC
      IF(IER.EQ.0) GO TO 999
      IF(RESULT.NE.0.0E+00.AND.PSUM(NUMRL2).NE.0.0E+00) GO TO 70
      IF(ABSERR.GT.ERRSUM) GO TO 80
      IF(PSUM(NUMRL2).EQ.0.0E+00) GO TO 999
   70 IF(ABSERR/ABS(RESULT).GT.(ERRSUM+DRL)/ABS(PSUM(NUMRL2)))
     1  GO TO 80
      IF(IER.GE.1.AND.IER.NE.7) ABSERR = ABSERR+DRL
      GO TO 999
   80 RESULT = PSUM(NUMRL2)
      ABSERR = ERRSUM+DRL
  999 RETURN
      END
