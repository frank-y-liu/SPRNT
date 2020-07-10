      SUBROUTINE PSTG3D (LDIMF,MDIMF,F,W)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
      DIMENSION       F(LDIMF,MDIMF,*),W(*)
      L=W(1)
      LP=W(2)
      M=W(3)
      MP=W(4)
      N=W(5)
      NP=W(6)
C
C     ALLOCATION OF WORK ARRAY W
C
      IA=7
      IC=IA+2*L
      ICFY=IC+L
      ICFZ=ICFY+4*M
      IFCTRD=ICFZ+4*N
      IWSY=IFCTRD+L*M*N
      IWSZ=IWSY+M+15
      IFT=IWSZ+N+15
C     IEND=IFT+L*M*N
      GO TO (105,114),LP
C
C     REORDER UNKNOWNS WHEN LPEROD = 0.
C
  105 LH = (L+1)/2
      LODD = 1
      IF (2*LH .EQ. L) LODD = 2
      DO 111 J=1,M
         DO 110 K=1,N
           DO 106 I=1,LH-1
               W(I+IFT) = F(LH-I,J,K)-F(LH+I,J,K)
               W(LH+I+IFT) = F(LH-I,J,K)+F(LH+I,J,K)
  106       CONTINUE
            W(LH+IFT) = 2.*F(LH,J,K)
            GO TO (108,107),LODD
  107       W(L+IFT) = 2.*F(L,J,K)
  108       DO 109 I=1,L
               F(I,J,K) = W(I+IFT)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
  114 CONTINUE
      IF (LDIMF.EQ.L .AND. MDIMF.EQ.M) GO TO 300
        CALL P3PACK(F,LDIMF,MDIMF,L,M,N,W(IFT))
      CALL PST3D1 (L,M,MP,N,NP,W(IFT),W(ICFY),W(ICFZ),F,
     1             W(IA),W(IC),W(IFCTRD),W(IWSY),W(IWSZ))
        CALL P3UNPK(F,LDIMF,MDIMF,L,M,N,W(IFT))
       GO TO 400
  300 CALL PST3D1 (L,M,MP,N,NP,F,W(ICFY),W(ICFZ),W(IFT),
     *             W(IA),W(IC),W(IFCTRD),W(IWSY),W(IWSZ))
  400 CONTINUE
      GO TO (115,122),LP
  115 DO 121 J=1,M
         DO 120 K=1,N
            DO 116 I=1,LH-1
               W(LH-I+IFT) = .5*(F(LH+I,J,K)+F(I,J,K))
               W(LH+I+IFT) = .5*(F(LH+I,J,K)-F(I,J,K))
  116       CONTINUE
            W(LH+IFT) = .5*F(LH,J,K)
            GO TO (118,117),LODD
  117       W(L+IFT) = .5*F(L,J,K)
  118       DO 119 I=1,L
               F(I,J,K) = W(I+IFT)
  119       CONTINUE
  120    CONTINUE
  121 CONTINUE
  122 CONTINUE
      RETURN
      END
