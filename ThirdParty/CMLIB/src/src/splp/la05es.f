      SUBROUTINE LA05ES(A, IRN, IP, N, IW, IA, REALS)
C***BEGIN PROLOGUE  LA05ES
C***REFER TO  SPLP
C     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
C     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
C     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
C     THE FINAL LETTER =S= IN THE NAMES USED HERE.
C     REVISED SEP. 13, 1979.
C
C     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
C     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
C     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
C     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
C     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    LA05DS
C***END PROLOGUE  LA05ES
      LOGICAL REALS
      REAL A(IA)
      INTEGER IRN(IA), IW(N)
      INTEGER IP(N)
      COMMON /LA05DS/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
C***FIRST EXECUTABLE STATEMENT  LA05ES
      NCP = NCP + 1
C     COMPRESS FILE OF POSITIVE INTEGERS. ENTRY J STARTS AT IRN(IP(J))
C  AND CONTAINS IW(J) INTEGERS,J=1,N. OTHER COMPONENTS OF IRN ARE ZERO.
C  LENGTH OF COMPRESSED FILE PLACED IN LROW IF REALS IS .TRUE. OR LCOL
C  OTHERWISE.
C  IF REALS IS .TRUE. ARRAY A CONTAINS A REAL FILE ASSOCIATED WITH IRN
C  AND THIS IS COMPRESSED TOO.
C  A,IRN,IP,IW,IA ARE INPUT/OUTPUT VARIABLES.
C  N,REALS ARE INPUT/UNCHANGED VARIABLES.
C
      DO 10 J=1,N
C STORE THE LAST ELEMENT OF ENTRY J IN IW(J) THEN OVERWRITE IT BY -J.
         NZ = IW(J)
         IF (NZ.LE.0) GO TO 10
         K = IP(J) + NZ - 1
         IW(J) = IRN(K)
         IRN(K) = -J
   10 CONTINUE
C KN IS THE POSITION OF NEXT ENTRY IN COMPRESSED FILE.
      KN = 0
      IPI = 0
      KL = LCOL
      IF (REALS) KL = LROW
C LOOP THROUGH THE OLD FILE SKIPPING ZERO (DUMMY) ELEMENTS AND
C     MOVING GENUINE ELEMENTS FORWARD. THE ENTRY NUMBER BECOMES
C     KNOWN ONLY WHEN ITS END IS DETECTED BY THE PRESENCE OF A NEGATIVE
C     INTEGER.
      DO 30 K=1,KL
         IF (IRN(K).EQ.0) GO TO 30
         KN = KN + 1
         IF (REALS) A(KN) = A(K)
         IF (IRN(K).GE.0) GO TO 20
C END OF ENTRY. RESTORE IRN(K), SET POINTER TO START OF ENTRY AND
C     STORE CURRENT KN IN IPI READY FOR USE WHEN NEXT LAST ENTRY
C     IS DETECTED.
         J = -IRN(K)
         IRN(K) = IW(J)
         IP(J) = IPI + 1
         IW(J) = KN - IPI
         IPI = KN
   20    IRN(KN) = IRN(K)
   30 CONTINUE
      IF (REALS) LROW = KN
      IF (.NOT.REALS) LCOL = KN
      RETURN
      END