      SUBROUTINE DDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)
C***BEGIN PROLOGUE  DDZRO
C***REFER TO  DDRIV3
C     This is a special purpose version of ZEROIN, modified for use with
C     the DDRIV1 package.
C
C     Sandia Mathematical Program Library
C     Mathematical Computing Services Division 5422
C     Sandia Laboratories
C     P. O. Box 5800
C     Albuquerque, New Mexico  87115
C     Control Data 6600 Version 4.5, 1 November 1971
C
C     ABSTRACT
C        ZEROIN searches for a zero of a function F(N, T, Y, IROOT)
C        between the given values B and C until the width of the
C        interval (B, C) has collapsed to within a tolerance specified
C        by the stopping criterion, ABS(B - C) .LE. 2.*(RW*ABS(B) + AE).
C
C     Description of parameters
C        F     - Name of the external function, which returns a
C                double precision result.  This name must be in an
C                EXTERNAL statement in the calling program.
C        B     - One end of the interval (B, C).  The value returned for
C                B usually is the better approximation to a zero of F.
C        C     - The other end of the interval (B, C).
C        RE    - Relative error used for RW in the stopping criterion.
C                If the requested RE is less than machine precision,
C                then RW is set to approximately machine precision.
C        AE    - Absolute error used in the stopping criterion.  If the
C                given interval (B, C) contains the origin, then a
C                nonzero value should be chosen for AE.
C
C     REFERENCES
C       1.  L F Shampine and H A Watts, ZEROIN, A Root-Solving Routine,
C           SC-TM-70-631, Sept 1970.
C       2.  T J Dekker, Finding a Zero by Means of Successive Linear
C           Interpolation, "Constructive Aspects of the Fundamental
C           Theorem of Algebra", edited by B Dejon and P Henrici, 1969.
C***ROUTINES CALLED  DDNTP
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDZRO
      DOUBLE PRECISION A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
     8     H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  DDZRO
      ER = 4.D0*UROUND
      RW = MAX(RE, ER)
      IC = 0
      ACBS = ABS(B - C)
      A = C
      FA = FC
      KOUNT = 0
C                                                    Perform interchange
 10   IF (ABS(FC) .LT. ABS(FB)) THEN
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
      END IF
      CMB = 0.5D0*(C - B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AE
C                                                Test stopping criterion
      IF (ACMB .LE. TOL) RETURN
      IF (KOUNT .GT. 50) RETURN
C                                    Calculate new iterate implicitly as
C                                    B + P/Q, where we arrange P .GE. 0.
C                         The implicit form is used to prevent overflow.
      P = (B - A)*FB
      Q = FA - FB
      IF (P .LT. 0.D0) THEN
        P = -P
        Q = -Q
      END IF
C                          Update A and check for satisfactory reduction
C                          in the size of our bounding interval.
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.D0*ACMB .GE. ACBS) THEN
C                                                                 Bisect
          B = 0.5D0*(C + B)
          GO TO 20
        END IF
        IC = 0
      END IF
      ACBS = ACMB
C                                            Test for too small a change
      IF (P .LE. ABS(Q)*TOL) THEN
C                                                 Increment by tolerance
        B = B + SIGN(TOL, CMB)
C                                               Root ought to be between
C                                               B and (C + B)/2.
      ELSE IF (P .LT. CMB*Q) THEN
C                                                            Interpolate
        B = B + P/Q
      ELSE
C                                                                 Bisect
        B = 0.5D0*(C + B)
      END IF
C                                             Have completed computation
C                                             for new iterate B.
 20   CALL DDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (FB .EQ. 0.D0) RETURN
      KOUNT = KOUNT + 1
C
C             Decide whether next step is interpolation or extrapolation
C
      IF (SIGN(1.0D0, FB) .EQ. SIGN(1.0D0, FC)) THEN
        C = A
        FC = FA
      END IF
      GO TO 10
      END
