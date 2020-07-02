      SUBROUTINE RPZERO(N,A,R,T,IFLG,S)
C***BEGIN PROLOGUE  RPZERO
C***DATE WRITTEN   810223   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  F1A1A
C***KEYWORDS  POLYNOMIAL ZEROS,REAL ROOTS,ROOTS,ZEROES,ZEROS
C***AUTHOR  KAHANER, D. K., (NBS)
C***PURPOSE  Find the zeros of a polynomial with real coefficients.
C***DESCRIPTION
C
C      Find the zeros of the real polynomial
C         P(X)= A(1)*X**N + A(2)*X**(N-1) +...+ A(N+1)
C
C    Input...
C       N = degree of P(X)
C       A = real vector containing coefficients of P(X),
C            A(I) = coefficient of X**(N+1-I)
C       R = N word complex vector containing initial estimates for zeros
C            if these are known.
C       T = 6(N+1) word array used for temporary storage
C       IFLG = flag to indicate if initial estimates of
C              zeros are input.
C            If IFLG .EQ. 0, no estimates are input.
C            If IFLG .NE. 0, the vector R contains estimates of
C               the zeros
C       ** Warning ****** If estimates are input, they must
C                         be separated; that is, distinct or
C                         not repeated.
C       S = an N word array
C
C    Output...
C       R(I) = ith zero,
C       S(I) = bound for R(I) .
C       IFLG = error diagnostic
C    Error Diagnostics...
C       If IFLG .EQ. 0 on return, all is well.
C       If IFLG .EQ. 1 on return, A(1)=0.0 or N=0 on input.
C       If IFLG .EQ. 2 on return, the program failed to coverge
C                after 25*N iterations.  Best current estimates of the
C                zeros are in R(I).  Error bounds are not calculated.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CPZERO
C***END PROLOGUE  RPZERO
C
      COMPLEX T(*),R(*)
      REAL A(*),S(*)
C***FIRST EXECUTABLE STATEMENT  RPZERO
      N1=N+1
      DO 1 I=1,N1
      T(I)= CMPLX(A(I),0.0)
    1 CONTINUE
      CALL CPZERO(N,T,R,T(N+2),IFLG,S)
      RETURN
      END
