      SUBROUTINE OHTROL(Q,N,NRDA,DIAG,IRANK,DIV,TD)
C***BEGIN PROLOGUE  OHTROL
C***REFER TO  SUDS
C
C     For a rank deficient problem,additional orthogonal
C     HOUSEHOLDER transformations are applied to the left side
C     of Q to futher reduce the triangular form.
C     Thus, after application of the routines ORTHOR and OHTROL
C     to the original matrix,the result is a nonsingular
C     triangular matrix while the remainder of the matrix
C     has been zeroed out.
C***ROUTINES CALLED  SDOT
C***END PROLOGUE  OHTROL
      DIMENSION Q(NRDA,N),DIAG(N),DIV(N),TD(N)
C***FIRST EXECUTABLE STATEMENT  OHTROL
      NMIR=N-IRANK
      IRP=IRANK+1
      DO 30 K=1,IRANK
         KIR=IRP-K
         DIAGK=DIAG(KIR)
         SIG=(DIAGK*DIAGK)+SDOT(NMIR,Q(IRP,KIR),1,Q(IRP,KIR),1)
         DD=SIGN(SQRT(SIG),-DIAGK)
         DIV(KIR)=DD
         TDV=DIAGK-DD
         TD(KIR)=TDV
         IF (K .EQ. IRANK) GO TO 30
         KIRM=KIR-1
         SQD=DD*DIAGK-SIG
         DO 20 J=1,KIRM
            QS=((TDV*Q(KIR,J))+SDOT(NMIR,Q(IRP,J),1,Q(IRP,KIR),1))
     1               /SQD
            Q(KIR,J)=Q(KIR,J)+QS*TDV
            DO 10 L=IRP,N
   10          Q(L,J)=Q(L,J)+QS*Q(L,KIR)
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
