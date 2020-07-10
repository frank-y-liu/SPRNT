      SUBROUTINE QMOMO(ALFA,BETA,RI,RJ,RG,RH,INTEGR)
C***BEGIN PROLOGUE  QMOMO
C***DATE WRITTEN   810101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A2A1,C3A2
C***KEYWORDS  MODIFIED CHEBYSHEV MOMENTS
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  This routine computes modified CHEBSYSHEV moments. The K-th
C            modified CHEBYSHEV moment is defined as the integral over
C            (-1,1) of W(X)*T(K,X), where T(K,X) is the CHEBYSHEV
C            POLYNOMIAL of degree K.
C***DESCRIPTION
C
C        MODIFIED CHEBYSHEV MOMENTS
C        STANDARD FORTRAN SUBROUTINE
C        REAL VERSION
C
C        PARAMETERS
C           ALFA   - Real
C                    Parameter in the weight function W(X), ALFA.GT.(-1)
C
C           BETA   - Real
C                    Parameter in the weight functioN W(X), BETA.GT.(-1)
C
C           RI     - Real
C                    Vector of dimension 25
C                    RI(K) is the integral over (-1,1) of
C                    (1+X)**ALFA*T(K-1,X), K = 1, ..., 25.
C
C           RJ     - Real
C                    Vector of dimension 25
C                    RJ(K) is the integral over (-1,1) of
C                    (1-X)**BETA*T(K-1,X), K = 1, ..., 25.
C
C           RG     - Real
C                    Vector of dimension 25
C                    RG(K) is the integral over (-1,1) of
C                    (1+X)**ALFA*LOG((1+X)/2)*T(K-1,X), K = 1, ..., 25.
C
C           RH     - Real
C                    Vector of dimension 25
C                    RH(K) is the integral over (-1,1) of
C                    (1-X)**BETA*LOG((1-X)/2)*T(K-1,X), K = 1, ..., 25.
C
C           INTEGR - Integer
C                    Input parameter indicating the modified
C                    Moments to be computed
C                    Integr = 1 compute RI, RJ
C                           = 2 compute RI, RJ, RG
C                           = 3 compute RI, RJ, RH
C                           = 4 compute RI, RJ, RG, RH
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  QMOMO
C
      REAL ALFA,ALFP1,ALFP2,AN,ANM1,BETA,BETP1,
     1  BETP2,RALF,RBET,RG,RH,RI,RJ
      INTEGER I,IM1,INTEGR
C
      DIMENSION RG(25),RH(25),RI(25),RJ(25)
C
C
C***FIRST EXECUTABLE STATEMENT  QMOMO
      ALFP1 = ALFA+0.1E+01
      BETP1 = BETA+0.1E+01
      ALFP2 = ALFA+0.2E+01
      BETP2 = BETA+0.2E+01
      RALF = 0.2E+01**ALFP1
      RBET = 0.2E+01**BETP1
C
C           COMPUTE RI, RJ USING A FORWARD RECURRENCE RELATION.
C
      RI(1) = RALF/ALFP1
      RJ(1) = RBET/BETP1
      RI(2) = RI(1)*ALFA/ALFP2
      RJ(2) = RJ(1)*BETA/BETP2
      AN = 0.2E+01
      ANM1 = 0.1E+01
      DO 20 I=3,25
        RI(I) = -(RALF+AN*(AN-ALFP2)*RI(I-1))/
     1  (ANM1*(AN+ALFP1))
        RJ(I) = -(RBET+AN*(AN-BETP2)*RJ(I-1))/
     1  (ANM1*(AN+BETP1))
        ANM1 = AN
        AN = AN+0.1E+01
   20 CONTINUE
      IF(INTEGR.EQ.1) GO TO 70
      IF(INTEGR.EQ.3) GO TO 40
C
C           COMPUTE RG USING A FORWARD RECURRENCE RELATION.
C
      RG(1) = -RI(1)/ALFP1
      RG(2) = -(RALF+RALF)/(ALFP2*ALFP2)-RG(1)
      AN = 0.2E+01
      ANM1 = 0.1E+01
      IM1 = 2
      DO 30 I=3,25
        RG(I) = -(AN*(AN-ALFP2)*RG(IM1)-AN*RI(IM1)+ANM1*RI(I))/
     1  (ANM1*(AN+ALFP1))
        ANM1 = AN
        AN = AN+0.1E+01
        IM1 = I
   30 CONTINUE
      IF(INTEGR.EQ.2) GO TO 70
C
C           COMPUTE RH USING A FORWARD RECURRENCE RELATION.
C
   40 RH(1) = -RJ(1)/BETP1
      RH(2) = -(RBET+RBET)/(BETP2*BETP2)-RH(1)
      AN = 0.2E+01
      ANM1 = 0.1E+01
      IM1 = 2
      DO 50 I=3,25
        RH(I) = -(AN*(AN-BETP2)*RH(IM1)-AN*RJ(IM1)+
     1  ANM1*RJ(I))/(ANM1*(AN+BETP1))
        ANM1 = AN
        AN = AN+0.1E+01
        IM1 = I
   50 CONTINUE
      DO 60 I=2,25,2
        RH(I) = -RH(I)
   60 CONTINUE
   70 DO 80 I=2,25,2
        RJ(I) = -RJ(I)
   80 CONTINUE
   90 RETURN
      END
