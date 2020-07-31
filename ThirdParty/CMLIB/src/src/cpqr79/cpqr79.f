      SUBROUTINE CPQR79(NDEG,COEFF,ROOT,IERR,WORK)
C***BEGIN PROLOGUE  CPQR79
C***DATE WRITTEN   791201   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  F1A1B
C***KEYWORDS  COMPLEX POLYNOMIAL,POLYNOMIAL ZEROS,ZEROS
C***AUTHOR  VANDEVENDER, W. H., (SNLA)
C***PURPOSE  Find the zeros of a polynomial with complex coefficients
C***DESCRIPTION
C
C     Abstract
C         This routine computes all zeros of a polynomial
C         of degree NDEG with complex coefficients
C         by computing the eigenvalues of the companion matrix.
C
C     Description of Parameters
C         The user must dimension all arrays appearing in the call list
C              COEFF(NDEG+1), ROOT(NDEG), WORK(*)
C
C      --Input--
C        NDEG    degree of polynomial
C
C        COEFF   complex coefficients in descending order.  i.e.,
C               P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)
C
C        WORK    real work array of dimension at least 2*NDEG*(NDEG+1)
C
C     --Output--
C        ROOT    computed complex roots
C
C        IERR    Output Error Code
C             - Normal Code
C            0  means the roots were computed.
C             - Abnormal Codes
C            1  more than 30 QR interations on some
C               eigenvalue of the companion matrix
C            2  COEFF(1)=0.0
C            3  NDEG is invalid( .LE. 0)
C
C       Externals Called
C
C         Fortran
C          COMQR--EISPACK routine
C          XERROR--SLATEC error handler
C          CABS--Fortran library function
C         Intrinsic
C          REAL
C          AIMAG
C          CMPLX
C     This routine is an interface to an eigenvalue routine in EISPACK.
C     This interface was written by Walter H. Vandevender.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COMQR,XERROR
C***END PROLOGUE  CPQR79
C
      COMPLEX COEFF(*),ROOT(*),SCALE,C
      REAL WORK(*)
      INTEGER NDEG,IERR,K,KHR,KHI,KWR,KWI,KAD,KJ
C***FIRST EXECUTABLE STATEMENT  CPQR79
      IERR = 0
      IF(CABS(COEFF(1)).EQ.0.0) GO TO 5
      IF(NDEG.LE.0) GO TO 6
      IF(NDEG.EQ.1)  GOTO 7
      SCALE=1.0E0/COEFF(1)
      KHR=1
      KHI=KHR+NDEG*NDEG
      KWR=KHI+KHI-KHR
      KWI=KWR+NDEG
      DO 1 K=1,KWR
1     WORK(K)=0.0E0
C
      DO 2 K =1,NDEG
       KAD=(K-1)*NDEG+1
       C=SCALE*COEFF(K+1)
       WORK(KAD)= -REAL(C)
       KJ=KHI+KAD-1
       WORK(KJ)= -AIMAG(C)
       IF(K.EQ.NDEG) GO TO 2
         WORK(KAD+K)=1.0E0
2     CONTINUE
C
      CALL COMQR(NDEG,NDEG,1,NDEG,WORK(KHR),WORK(KHI),WORK(KWR),WORK(KWI
     1),IERR)
      IF(IERR.NE.0)  GO TO 4
      DO 3 K=1,NDEG
      KM1=K-1
3     ROOT(K)=CMPLX(WORK(KWR+KM1),WORK(KWI+KM1))
      RETURN
C
4     IERR = 1
      CALL XERROR( 'CPQR - NO CONVERGENCE IN 30 QR ITERATIONS.',42,8,1)
      RETURN
    5 IERR = 2
      CALL XERROR( 'CPQR - LEADING COEFFICIENT IS ZERO.',35,2,1)
      RETURN
6     IERR=3
      CALL XERROR( 'CPQR - DEGREE INVALID.' ,22,2,1)
      RETURN
7     ROOT(1)= -COEFF(2)/COEFF(1)
      RETURN
      END
