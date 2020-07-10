      SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)
C***BEGIN PROLOGUE  R1MPYQ
C***REFER TO  SNSQ,SNSQE
C
C     **********
C     Subroutine R1MPYQ
C
C     Given an M by N matrix A, this subroutine computes A*Q where
C     Q is the product of 2*(N - 1) transformations
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
C     eliminate elements in the I-th and N-th planes, respectively.
C     Q itself is not given, rather the information to recover the
C     GV, GW rotations is supplied.
C
C     The subroutine statement is
C
C       SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A.
C
C       N is a positive integer input variable set to the number
C         of columns of A.
C
C       A is an M by N ARRAY. On input A must contain the matrix
C         to be postmultiplied by the orthogonal matrix Q
C         described above. On output A*Q has replaced A.
C
C       LDA is a positive integer input variable not less than M
C         which specifies the leading dimension of the array A.
C
C       V is an input array of length N. V(I) must contain the
C         information necessary to recover the Givens rotation GV(I)
C         described above.
C
C       W is an input array of length N. W(I) must contain the
C         information necessary to recover the Givens rotation GW(I)
C         described above.
C
C     Subroutines called
C
C       FORTRAN-Supplied ... abs,sqrt
C
C     MINPACK. Version of December 1978.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C     **********
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  R1MPYQ
      INTEGER M,N,LDA
      REAL A(LDA,N),V(N),W(N)
      INTEGER I,J,NMJ,NM1
      REAL COS,ONE,SIN,TEMP
      DATA ONE /1.0E0/
C***FIRST EXECUTABLE STATEMENT  R1MPYQ
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 50
      DO 20 NMJ = 1, NM1
         J = N - NMJ
         IF (ABS(V(J)) .GT. ONE) COS = ONE/V(J)
         IF (ABS(V(J)) .GT. ONE) SIN = SQRT(ONE-COS**2)
         IF (ABS(V(J)) .LE. ONE) SIN = V(J)
         IF (ABS(V(J)) .LE. ONE) COS = SQRT(ONE-SIN**2)
         DO 10 I = 1, M
            TEMP = COS*A(I,J) - SIN*A(I,N)
            A(I,N) = SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
C
      DO 40 J = 1, NM1
         IF (ABS(W(J)) .GT. ONE) COS = ONE/W(J)
         IF (ABS(W(J)) .GT. ONE) SIN = SQRT(ONE-COS**2)
         IF (ABS(W(J)) .LE. ONE) SIN = W(J)
         IF (ABS(W(J)) .LE. ONE) COS = SQRT(ONE-SIN**2)
         DO 30 I = 1, M
            TEMP = COS*A(I,J) + SIN*A(I,N)
            A(I,N) = -SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE R1MPYQ.
C
      END
