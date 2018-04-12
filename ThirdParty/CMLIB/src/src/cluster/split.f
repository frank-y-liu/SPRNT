      SUBROUTINE SPLIT(MM, M, N, A, CLAB, RLAB, DMW, W, IL, IU, DMU, U,
     *                 WCLAB, IM, DM)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      SPLITS A ROW CLUSTER ON SELECTED VARIABLES
C
C   DESCRIPTION
C   -----------
C
C   1.  INITIALLY, THE FIRST CLUSTER CONSISTS OF ALL CASES BETWEEN IL
C       AND IU AND THE SECOND CLUSTER IS EMPTY.  THE WEIGHTED MEANS ARE
C       DETERMINED AND USED TO FIND THE REDUCTION IN THE WITHIN-CLUSTER
C       SUM OF SQUARES FOR MOVING EACH CASE FROM THE FIRST CLUSTER TO
C       THE SECOND.  THE OBJECT THAT REDUCES THE SUM OF SQUARES THE
C       MOST IS MOVED AND THIS CONTINUES UNTIL ALL OBJECTS ARE MOVED
C       WITH EACH REDUCTION STORED.  THEN THE SPLIT THAT HAD THE
C       SMALLEST REDUCTION OF ALL IS RETURNED AS THE OPTIMUM SPLIT.
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM, M, N, A, CLAB, RLAB, DMW, W -- SEE SUBROUTINE SPLIT1
C
C   IL, IU INTEGER SCALARS (UNCHANGED ON OUTPUT).
C         THE FIRST AND LAST OBJECTS IN THE BLOCK TO BE SPLIT.
C
C   DMU   INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE LEADING DIMENSION OF MATRIX U.  MUST BE AT LEAST 4.
C
C   U     REAL MATRIX WHOSE FIRST DIMENSION MUST BE DMU AND SECOND
C            DIMENSION MUST BE AT LEAST N (CHANGED ON OUTPUT).
C         MATRIX OF CLUSTER MEANS.
C
C   OUTPUT PARAMETERS
C   -----------------
C
C   WCLAB INTEGER VECTOR DIMENSIONED AT LEAST N.
C         WCLAB(I) WILL STORE THE CLUSTER (EITHER 1 OR 2) OBJECT I WAS
C            ASSIGNED TO.
C
C   IM    INTEGER SCALAR.
C         THE BORDER OF THE SPLIT.  OBJECTS IL,...,IM WERE ASSIGNED TO
C            CLUSTER 1 AND OBJECTS IM+1,...,IU WERE ASSIGNED TO CLUSTER
C            2.
C
C   DM    INTEGER SCALAR.
C         THE REDUCTION IN THE WITHIN-CLUSTER SUM OF SQUARES.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 272.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      INTEGER DMW, DMU, WCLAB(*)
      DIMENSION W(DMW,*), A(MM,*), U(DMU,*)
      CHARACTER*4 CLAB(*), RLAB(*), CTEMP
C
C     FIND WEIGHTED MEAN OF ALL CASES
C
      TH=R1MACH(4)
      DO 10 J=1,N
         U(1,J)=0.
         U(3,J)=0.
         U(2,J)=TH
         U(4,J)=TH
   10 CONTINUE
      DO 30 J=1,N
         IF(WCLAB(J).NE.0) THEN
            DO 20 I=IL,IU
               U(1,J)=U(1,J)+A(I,J)*W(I,J)
   20          U(2,J)=U(2,J)+WCLAB(J)
            U(1,J)=U(1,J)/U(2,J)
         ENDIF
   30 CONTINUE
      DM=0.
      DD=0.
      DO 80 IC=IL,IU
         II=IU-IC+IL
         DMAX=-R1MACH(2)
         IMAX=II
C
C     DETERMINE THE EFFECT OF MOVING ITH CASE
C
         DO 50 I=IL,II
            D=0.
            DO 40 J=1,N
               IF(WCLAB(J).NE.0) THEN
                 IF(U(2,J).EQ.W(I,J)) U(2,J)=W(I,J)+TH
                 D=D+W(I,J)*U(2,J)*(A(I,J)-U(1,J))**2/(U(2,J)-W(I,J))
                 D=D-W(I,J)*U(4,J)*(A(I,J)-U(3,J))**2/(U(4,J)+W(I,J))
               ENDIF
   40       CONTINUE
C
C     STORE THE LARGEST
C
            IF(D.GT.DMAX) THEN
               IMAX=I
               DMAX=D
            ENDIF
   50    CONTINUE
         DD=DD+DMAX
         IF(DD.GT.DM) IM=II-1
         IF(DD.GT.DM) DM=DD
C
C     UPDATE MEANS OF THE TWO CLUSTERS
C
         I=IMAX
         DO 60 J=1,N
            IF(WCLAB(J).NE.0) THEN
               U(2,J)=U(2,J)-W(I,J)
               IF(U(2,J).LT.TH) U(2,J)=TH
               U(1,J)=U(1,J)+(U(1,J)-A(I,J))*W(I,J)/U(2,J)
               U(4,J)=U(4,J)+W(I,J)
               U(3,J)=U(3,J)-(U(3,J)-A(I,J))*W(I,J)/U(4,J)
            ENDIF
   60    CONTINUE
C
C     INTERCHANGE SELECTED ROW WITH LAST FEASIBLE ROW
C
         DO 70 J=1,N
            C=A(I,J)
            A(I,J)=A(II,J)
            A(II,J)=C
            C=W(I,J)
            W(I,J)=W(II,J)
   70       W(II,J)=C
         CTEMP = RLAB(I)
         RLAB(I) = RLAB(II)
         RLAB(II) = CTEMP
   80 CONTINUE
      RETURN
      END
