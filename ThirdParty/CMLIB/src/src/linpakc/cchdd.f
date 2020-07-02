      SUBROUTINE CCHDD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S,INFO)
C***BEGIN PROLOGUE  CCHDD
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D7B
C***KEYWORDS  CHOLESKY DECOMPOSITION,COMPLEX,DOWNDATE,LINEAR ALGEBRA,
C             LINPACK,MATRIX
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Downdates an augmented Cholesky decomposition or the
C            triangular factor of an augmented QR decomposition.
C***DESCRIPTION
C
C     CCHDD downdates an augmented Cholesky decomposition or the
C     triangular factor of an augmented QR decomposition.
C     Specifically, given an upper triangular matrix R of order P,  a
C     row vector X, a column vector Z, and a scalar Y, CCHDD
C     determines a unitary matrix U and a scalar ZETA such that
C
C                        (R   Z )     (RR  ZZ)
C                    U * (      )  =  (      ) ,
C                        (0 ZETA)     ( X   Y)
C
C     where RR is upper triangular.  If R and Z have been obtained
C     from the factorization of a least squares problem, then
C     RR and ZZ are the factors corresponding to the problem
C     with the observation (X,Y) removed.  In this case, if RHO
C     is the norm of the residual vector, then the norm of
C     the residual vector of the downdated problem is
C     SQRT(RHO**2 - ZETA**2).  CCHDD will simultaneously downdate
C     several triplets (Z,Y,RHO) along with R.
C     For a less terse description of what CCHDD does and how
C     it may be applied, see the LINPACK Guide.
C
C     The matrix U is determined as the product U(1)*...*U(P)
C     where U(I) is a rotation in the (P+1,I)-plane of the
C     form
C
C                       ( C(I)  -CONJG(S(I)) )
C                       (                    ) .
C                       ( S(I)       C(I)    )
C
C     the rotations are chosen so that C(I) is real.
C
C     The user is warned that a given downdating problem may
C     be impossible to accomplish or may produce
C     inaccurate results.  For example, this can happen
C     if X is near a vector whose removal will reduce the
C     rank of R.  Beware.
C
C     On Entry
C
C         R      COMPLEX(LDR,P), where LDR .GE. P.
C                R contains the upper triangular matrix
C                that is to be downdated.  The part of R
C                below the diagonal is not referenced.
C
C         LDR    INTEGER.
C                LDR is the leading dimension of the array R.
C
C         p      INTEGER.
C                P is the order of the matrix R.
C
C         X      COMPLEX(P).
C                X contains the row vector that is to
C                be removed from R.  X is not altered by CCHDD.
C
C         Z      COMPLEX(LDZ,NZ), where LDZ .GE. P.
C                Z is an array of NZ P-vectors which
C                are to be downdated along with R.
C
C         LDZ    INTEGER.
C                LDZ is the leading dimension of the array Z.
C
C         NZ     INTEGER.
C                NZ is the number of vectors to be downdated
C                NZ may be zero, in which case Z, Y, and RHO
C                are not referenced.
C
C         Y      COMPLEX(NZ).
C                Y contains the scalars for the downdating
C                of the vectors Z.  Y is not altered bY CCHDD.
C
C         RHO    REAL(NZ).
C                RHO contains the norms of the residual
C                vectors that are to be downdated.
C
C     On Return
C
C         R
C         Z      contain the downdated quantities.
C         RHO
C
C         C      REAL(P).
C                C contains the cosines of the transforming
C                rotations.
C
C         S      COMPLEX(P).
C                S contains the sines of the transforming
C                rotations.
C
C         INFO   INTEGER.
C                INFO is set as follows.
C
C                   INFO = 0  if the entire downdating
C                             was successful.
C
C                   INFO =-1  if R could not be downdated.
C                             in this case, all quantities
C                             are left unaltered.
C
C                   INFO = 1  if some RHO could not be
C                             downdated.  The offending RHO's are
C                             set to -1.
C
C     LINPACK.  This version dated 08/14/78 .
C     Stewart, G. W., University of Maryland, Argonne National Lab.
C
C     CCHDD uses the following functions and subprograms.
C
C     Fortran AIMAG,CABS,CONJG
C     BLAS CDOTC, SCNRM2
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  CDOTC,SCNRM2
C***END PROLOGUE  CCHDD
      INTEGER LDR,P,LDZ,NZ,INFO
      COMPLEX R(LDR,*),X(*),Z(LDZ,*),Y(*),S(*)
      REAL RHO(*),C(*)
C
      INTEGER I,II,J
      REAL A,ALPHA,AZETA,NORM,SCNRM2
      COMPLEX CDOTC,T,ZETA,B,XX
C
C     SOLVE THE SYSTEM CTRANS(R)*A = X, PLACING THE RESULT
C     IN THE ARRAY S.
C
C***FIRST EXECUTABLE STATEMENT  CCHDD
      INFO = 0
      S(1) = CONJG(X(1))/CONJG(R(1,1))
      IF (P .LT. 2) GO TO 20
      DO 10 J = 2, P
         S(J) = CONJG(X(J)) - CDOTC(J-1,R(1,J),1,S,1)
         S(J) = S(J)/CONJG(R(J,J))
   10 CONTINUE
   20 CONTINUE
      NORM = SCNRM2(P,S,1)
      IF (NORM .LT. 1.0E0) GO TO 30
         INFO = -1
      GO TO 120
   30 CONTINUE
         ALPHA = SQRT(1.0E0-NORM**2)
C
C        DETERMINE THE TRANSFORMATIONS.
C
         DO 40 II = 1, P
            I = P - II + 1
            SCALE = ALPHA + CABS(S(I))
            A = ALPHA/SCALE
            B = S(I)/SCALE
            NORM = SQRT(A**2+REAL(B)**2+AIMAG(B)**2)
            C(I) = A/NORM
            S(I) = CONJG(B)/NORM
            ALPHA = SCALE*NORM
   40    CONTINUE
C
C        APPLY THE TRANSFORMATIONS TO R.
C
         DO 60 J = 1, P
            XX = (0.0E0,0.0E0)
            DO 50 II = 1, J
               I = J - II + 1
               T = C(I)*XX + S(I)*R(I,J)
               R(I,J) = C(I)*R(I,J) - CONJG(S(I))*XX
               XX = T
   50       CONTINUE
   60    CONTINUE
C
C        IF REQUIRED, DOWNDATE Z AND RHO.
C
         IF (NZ .LT. 1) GO TO 110
         DO 100 J = 1, NZ
            ZETA = Y(J)
            DO 70 I = 1, P
               Z(I,J) = (Z(I,J) - CONJG(S(I))*ZETA)/C(I)
               ZETA = C(I)*ZETA - S(I)*Z(I,J)
   70       CONTINUE
            AZETA = CABS(ZETA)
            IF (AZETA .LE. RHO(J)) GO TO 80
               INFO = 1
               RHO(J) = -1.0E0
            GO TO 90
   80       CONTINUE
               RHO(J) = RHO(J)*SQRT(1.0E0-(AZETA/RHO(J))**2)
   90       CONTINUE
  100    CONTINUE
  110    CONTINUE
  120 CONTINUE
      RETURN
      END
