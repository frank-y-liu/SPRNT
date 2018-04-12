      SUBROUTINE RPQR79(NDEG,COEFF,ROOT,IERR,WORK)
C***BEGIN PROLOGUE  RPQR79
C***DATE WRITTEN   800601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  F1A1A
C***KEYWORDS  POLYNOMIAL ROOTS,REAL,ROOTS,ZEROES,ZEROS
C***AUTHOR  VANDEVENDER, W. H., (SNLA)
C***PURPOSE  To find the zeros of a polynomial with real coefficients.
C***DESCRIPTION
C
C     This routine is an interface to an eigenvalue routine in EISPACK.
C     This interface was written by Walter H. Vandevender.
C
C     Abstract
C
C         This routine computes all roots of a polynomial with real
C         coefficients by computing the eigenvalues of the
C         companion matrix.
C
C     Description of Parameters
C         The user must dimension all arrays appearing in the call list
C              COEFF(NDEG+1),ROOT(NDEG),WORK(NDEG*NDEG+2*NDEG)
C
C      Input -
C       NDEG   degree of polynomial
C
C       COEFF  array of coeffficients in order of descending powers of
C              Z,   i.e. COEFF(1)*(Z**NDEG) + ... + COEFF(NDEG+1)
C
C       WORK   real work vector of length at least NDEG*(NDEG+2)
C
C     Output-
C       ROOT   complex vector of computed roots.
C
C       IERR   output error code
C             - normal code
C              0  means the roots were computed.
C             - abnormal codes
C              1  more than 30 QR iterations on some eigenvalue
C                 of the companion matrix
C              2  COEF(1) = 0.0
C              3  NDEG invalid (less than 0).
C***REFERENCES  (NONE)
C***ROUTINES CALLED  HQR,XERROR
C***END PROLOGUE  RPQR79
      REAL COEFF(*),WORK(*),SCALE
      COMPLEX ROOT(*)
      INTEGER NDEG,IERR,K,J,KH,KWR,KWI,KCOL
C***FIRST EXECUTABLE STATEMENT  RPQR79
      IERR = 0
      IF(COEFF(1).EQ.0.0)   GOTO 5
      IF(NDEG.LT.0) GOTO 6
      IF(NDEG.EQ.1) GO TO 7
      SCALE=1.E0/COEFF(1)
      KH=1
      KWR=KH+NDEG*NDEG
      KWI=KWR+NDEG
      KWEND=KWI+NDEG-1
      DO 1 K=1,KWEND
1     WORK(K)=0.E0
      DO 2 K=1,NDEG
       KCOL=(K-1)*NDEG+1
       WORK(KCOL)=-COEFF(K+1)*SCALE
       IF(K.EQ.NDEG) GO TO 2
        WORK(KCOL+K)=1.E0
2     CONTINUE
      CALL HQR(NDEG,NDEG,1,NDEG,WORK(KH),WORK(KWR),WORK(KWI),IERR)
      IF(IERR.NE.0) GO TO 4
      DO 3 K=1,NDEG
       KM1=K-1
3     ROOT(K)=CMPLX(WORK(KWR+KM1),WORK(KWI+KM1))
      RETURN
C
4     IERR = 1
      CALL XERROR ( 'RPQR - NO CONVERGENCE IN 30 QR ITERATIONS.', 42,
     1 8,1)
      RETURN
5     IERR = 2
      CALL XERROR ( 'RPQR - LEADING COEFFICIENT IS ZERO.',35,2,1)
      RETURN
6     IERR = 3
      CALL XERROR ( 'RPQR - DEGREE INVALID.' ,22,2,1)
      RETURN
7      ROOT(1)=CMPLX(-COEFF(2)/COEFF(1),0.)
      RETURN
      END
