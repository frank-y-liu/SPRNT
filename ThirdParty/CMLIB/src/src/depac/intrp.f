      SUBROUTINE INTRP(X,Y,XOUT,YOUT,YPOUT,NEQN,KOLD,PHI,PSI)
C***BEGIN PROLOGUE  INTRP
C***DATE WRITTEN   740101   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ADAMS METHOD,DEPAC,INITIAL VALUE PROBLEMS,ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS,PREDICTOR-CORRECTOR
C***AUTHOR  SHAMPINE, L. F., (SNLA)
C           GORDON, M. K., (SNLA)
C***PURPOSE  Approximates the solution at XOUT by evaluating the
C            polynomial computed in STEP2 at XOUT.  Must be used in
C            conjunction with STEP2.
C***DESCRIPTION
C
C   Written by L. F. Shampine and M. K. Gordon
C
C   Abstract
C
C
C   The methods in subroutine  STEP2  approximate the solution near  X
C   by a polynomial.  Subroutine  INTRP  approximates the solution at
C   XOUT  by evaluating the polynomial there.  Information defining this
C   polynomial is passed from  STEP2  so  INTRP  cannot be used alone.
C
C   This code is completely explained and documented in the text,
C   "Computer Solution of Ordinary Differential Equations, the Initial
C   Value Problem"  by L. F. Shampine and M. K. Gordon.
C   Further details on use of this code are available in "Solving
C   Ordinary Differential Equations with ODE, STEP, and INTRP",
C   by L. F. Shampine and M. K. Gordon, SLA-73-1060.
C
C   Input to INTRP --
C
C   The user provides storage in the calling program for the arrays in
C   the call list
C      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),PSI(12)
C   and defines
C      XOUT -- point at which solution is desired.
C   The remaining parameters are defined in  STEP2  and passed to
C   INTRP  from that subroutine
C
C   Output from  INTRP --
C
C      YOUT(*) -- solution at  XOUT
C      YPOUT(*) -- derivative of solution at  XOUT
C   The remaining parameters are returned unaltered from their input
C   values.  Integration with  STEP2  may be continued.
C***REFERENCES  SHAMPINE, L. F., GORDON, M. K., *SOLVING
C                 ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP,
C                 AND INTRP*, SLA-73-1060, SANDIA LABORATORIES,
C                 1973.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  INTRP
C
C
      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),PSI(12)
      DIMENSION G(13),W(13),RHO(13)
      DATA G(1)/1.0/,RHO(1)/1.0/
C
C***FIRST EXECUTABLE STATEMENT  INTRP
      HI = XOUT - X
      KI = KOLD + 1
      KIP1 = KI + 1
C
C   INITIALIZE W(*) FOR COMPUTING G(*)
C
      DO 5 I = 1,KI
        TEMP1 = I
 5      W(I) = 1.0/TEMP1
      TERM = 0.0
C
C   COMPUTE G(*)
C
      DO 15 J = 2,KI
        JM1 = J - 1
        PSIJM1 = PSI(JM1)
        GAMMA = (HI + TERM)/PSIJM1
        ETA = HI/PSIJM1
        LIMIT1 = KIP1 - J
        DO 10 I = 1,LIMIT1
 10       W(I) = GAMMA*W(I) - ETA*W(I+1)
        G(J) = W(1)
        RHO(J) = GAMMA*RHO(JM1)
 15     TERM = PSIJM1
C
C   INTERPOLATE
C
      DO 20 L = 1,NEQN
        YPOUT(L) = 0.0
 20     YOUT(L) = 0.0
      DO 30 J = 1,KI
        I = KIP1 - J
        TEMP2 = G(I)
        TEMP3 = RHO(I)
        DO 25 L = 1,NEQN
          YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
 25       YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
 30     CONTINUE
      DO 35 L = 1,NEQN
 35     YOUT(L) = Y(L) + HI*YOUT(L)
      RETURN
      END
