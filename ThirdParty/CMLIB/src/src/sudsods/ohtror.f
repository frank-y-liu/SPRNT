      SUBROUTINE OHTROR(Q,N,NRDA,DIAG,IRANK,DIV,TD)
C***BEGIN PROLOGUE  OHTROR
C***REFER TO  SODS
C
C     For a rank deficient problem,additional orthogonal
C     HOUSEHOLDER transformations are applied to the right side
C     of Q to futher reduce the triangular form.
C     Thus, after application of the routines ORTHOL and OHTROR
C     to the original matrix,the result is a nonsingular
C     triangular matrix while the remainder of the matrix
C     has been zeroed out.
C***ROUTINES CALLED  SDOT
C***END PROLOGUE  OHTROR
      DIMENSION Q(NRDA,N),DIAG(N),DIV(N),TD(N)
C***FIRST EXECUTABLE STATEMENT  OHTROR
      NMIR=N-IRANK
      IRP=IRANK+1
      DO 30 K=1,IRANK
         KIR=IRP-K
         DIAGK=DIAG(KIR)
         SIG=(DIAGK*DIAGK)+SDOT(NMIR,Q(KIR,IRP),NRDA,Q(KIR,IRP),NRDA)
         DD=SIGN(SQRT(SIG),-DIAGK)
         DIV(KIR)=DD
         TDV=DIAGK-DD
         TD(KIR)=TDV
         IF (K .EQ. IRANK) GO TO 30
         KIRM=KIR-1
         SQD=DD*DIAGK-SIG
         DO 20 J=1,KIRM
            QS=((TDV*Q(J,KIR))+SDOT(NMIR,Q(J,IRP),NRDA,Q(KIR,IRP),NRDA))
     1               /SQD
            Q(J,KIR)=Q(J,KIR)+QS*TDV
            DO 10 L=IRP,N
   10          Q(J,L)=Q(J,L)+QS*Q(KIR,L)
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
