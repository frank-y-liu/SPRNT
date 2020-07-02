      SUBROUTINE DPARCK(ALG, D, IV, LIV, LV, N, V)
C
C  ***  CHECK ***SOL (VERSION 2.3) PARAMETERS, PRINT CHANGED VALUES  ***
C
C  ***  ALG = 1 FOR REGRESSION, ALG = 2 FOR GENERAL UNCONSTRAINED OPT.
C
      INTEGER ALG, LIV, LV, N
      INTEGER IV(LIV)
      DOUBLE PRECISION D(N), V(LV)
C
      EXTERNAL  DVDFLT
      DOUBLE PRECISION D1MACH
C DVDFLT  -- SUPPLIES DEFAULT PARAMETER VALUES TO V ALONE.
C/+
      INTEGER MAX0
C/
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, II, IV1, J, K, L, M, MIV1, MIV2, NDFALT, PARSV1, PU
      INTEGER IJMP, JLIM(2), MINIV(2), NDFLT(2)
C/6
      INTEGER VARNM(2), SH(2)
      REAL CNGD(3), DFLT(3), VN(2,34), WHICH(3)
C/7
C     CHARACTER*1 VARNM(2), SH(2)
C     CHARACTER*4 CNGD(3), DFLT(3), VN(2,34), WHICH(3)
C/
      DOUBLE PRECISION BIG, MACHEP, TINY, VK, VM(34), VX(34), ZERO
C
C  ***  IV AND V SUBSCRIPTS  ***
C
      INTEGER ALGSAV, DINIT, DTYPE, DTYPE0, EPSLON, INITS, IVNEED,
     1        LASTIV, LASTV, LMAT, NEXTIV, NEXTV, NVDFLT, OLDN,
     2        PARPRT, PARSAV, PERM, PRUNIT, VNEED
C
C
C/6
      DATA ALGSAV/51/, DINIT/38/, DTYPE/16/, DTYPE0/54/, EPSLON/19/,
     1     INITS/25/, IVNEED/3/, LASTIV/44/, LASTV/45/, LMAT/42/,
     2     NEXTIV/46/, NEXTV/47/, NVDFLT/50/, OLDN/38/, PARPRT/20/,
     3     PARSAV/49/, PERM/58/, PRUNIT/21/, VNEED/4/
C/7
C     PARAMETER (ALGSAV=51, DINIT=38, DTYPE=16, DTYPE0=54, EPSLON=19,
C    1           INITS=25, IVNEED=3, LASTIV=44, LASTV=45, LMAT=42,
C    2           NEXTIV=46, NEXTV=47, NVDFLT=50, OLDN=38, PARPRT=20,
C    3           PARSAV=49, PERM=58, PRUNIT=21, VNEED=4)
C     SAVE BIG, MACHEP, TINY
C/
C
      DATA BIG/0.D+0/, MACHEP/-1.D+0/, TINY/1.D+0/, ZERO/0.D+0/
C/6
      DATA VN(1,1),VN(2,1)/4HEPSL,4HON../
      DATA VN(1,2),VN(2,2)/4HPHMN,4HFC../
      DATA VN(1,3),VN(2,3)/4HPHMX,4HFC../
      DATA VN(1,4),VN(2,4)/4HDECF,4HAC../
      DATA VN(1,5),VN(2,5)/4HINCF,4HAC../
      DATA VN(1,6),VN(2,6)/4HRDFC,4HMN../
      DATA VN(1,7),VN(2,7)/4HRDFC,4HMX../
      DATA VN(1,8),VN(2,8)/4HTUNE,4HR1../
      DATA VN(1,9),VN(2,9)/4HTUNE,4HR2../
      DATA VN(1,10),VN(2,10)/4HTUNE,4HR3../
      DATA VN(1,11),VN(2,11)/4HTUNE,4HR4../
      DATA VN(1,12),VN(2,12)/4HTUNE,4HR5../
      DATA VN(1,13),VN(2,13)/4HAFCT,4HOL../
      DATA VN(1,14),VN(2,14)/4HRFCT,4HOL../
      DATA VN(1,15),VN(2,15)/4HXCTO,4HL.../
      DATA VN(1,16),VN(2,16)/4HXFTO,4HL.../
      DATA VN(1,17),VN(2,17)/4HLMAX,4H0.../
      DATA VN(1,18),VN(2,18)/4HLMAX,4HS.../
      DATA VN(1,19),VN(2,19)/4HSCTO,4HL.../
      DATA VN(1,20),VN(2,20)/4HDINI,4HT.../
      DATA VN(1,21),VN(2,21)/4HDTIN,4HIT../
      DATA VN(1,22),VN(2,22)/4HD0IN,4HIT../
      DATA VN(1,23),VN(2,23)/4HDFAC,4H..../
      DATA VN(1,24),VN(2,24)/4HDLTF,4HDC../
      DATA VN(1,25),VN(2,25)/4HDLTF,4HDJ../
      DATA VN(1,26),VN(2,26)/4HDELT,4HA0../
      DATA VN(1,27),VN(2,27)/4HFUZZ,4H..../
      DATA VN(1,28),VN(2,28)/4HRLIM,4HIT../
      DATA VN(1,29),VN(2,29)/4HCOSM,4HIN../
      DATA VN(1,30),VN(2,30)/4HHUBE,4HRC../
      DATA VN(1,31),VN(2,31)/4HRSPT,4HOL../
      DATA VN(1,32),VN(2,32)/4HSIGM,4HIN../
      DATA VN(1,33),VN(2,33)/4HETA0,4H..../
      DATA VN(1,34),VN(2,34)/4HBIAS,4H..../
C/7
C     DATA VN(1,1),VN(2,1)/'EPSL','ON..'/
C     DATA VN(1,2),VN(2,2)/'PHMN','FC..'/
C     DATA VN(1,3),VN(2,3)/'PHMX','FC..'/
C     DATA VN(1,4),VN(2,4)/'DECF','AC..'/
C     DATA VN(1,5),VN(2,5)/'INCF','AC..'/
C     DATA VN(1,6),VN(2,6)/'RDFC','MN..'/
C     DATA VN(1,7),VN(2,7)/'RDFC','MX..'/
C     DATA VN(1,8),VN(2,8)/'TUNE','R1..'/
C     DATA VN(1,9),VN(2,9)/'TUNE','R2..'/
C     DATA VN(1,10),VN(2,10)/'TUNE','R3..'/
C     DATA VN(1,11),VN(2,11)/'TUNE','R4..'/
C     DATA VN(1,12),VN(2,12)/'TUNE','R5..'/
C     DATA VN(1,13),VN(2,13)/'AFCT','OL..'/
C     DATA VN(1,14),VN(2,14)/'RFCT','OL..'/
C     DATA VN(1,15),VN(2,15)/'XCTO','L...'/
C     DATA VN(1,16),VN(2,16)/'XFTO','L...'/
C     DATA VN(1,17),VN(2,17)/'LMAX','0...'/
C     DATA VN(1,18),VN(2,18)/'LMAX','S...'/
C     DATA VN(1,19),VN(2,19)/'SCTO','L...'/
C     DATA VN(1,20),VN(2,20)/'DINI','T...'/
C     DATA VN(1,21),VN(2,21)/'DTIN','IT..'/
C     DATA VN(1,22),VN(2,22)/'D0IN','IT..'/
C     DATA VN(1,23),VN(2,23)/'DFAC','....'/
C     DATA VN(1,24),VN(2,24)/'DLTF','DC..'/
C     DATA VN(1,25),VN(2,25)/'DLTF','DJ..'/
C     DATA VN(1,26),VN(2,26)/'DELT','A0..'/
C     DATA VN(1,27),VN(2,27)/'FUZZ','....'/
C     DATA VN(1,28),VN(2,28)/'RLIM','IT..'/
C     DATA VN(1,29),VN(2,29)/'COSM','IN..'/
C     DATA VN(1,30),VN(2,30)/'HUBE','RC..'/
C     DATA VN(1,31),VN(2,31)/'RSPT','OL..'/
C     DATA VN(1,32),VN(2,32)/'SIGM','IN..'/
C     DATA VN(1,33),VN(2,33)/'ETA0','....'/
C     DATA VN(1,34),VN(2,34)/'BIAS','....'/
C/
C
      DATA VM(1)/1.0D-3/, VM(2)/-0.99D+0/, VM(3)/1.0D-3/, VM(4)/1.0D-2/,
     1     VM(5)/1.2D+0/, VM(6)/1.D-2/, VM(7)/1.2D+0/, VM(8)/0.D+0/,
     2     VM(9)/0.D+0/, VM(10)/1.D-3/, VM(11)/-1.D+0/, VM(15)/0.D+0/,
     3     VM(16)/0.D+0/, VM(19)/0.D+0/, VM(20)/-10.D+0/, VM(21)/0.D+0/,
     4     VM(22)/0.D+0/, VM(23)/0.D+0/, VM(27)/1.01D+0/,
     5     VM(28)/1.D+10/, VM(30)/0.D+0/, VM(31)/0.D+0/, VM(32)/0.D+0/,
     6     VM(34)/0.D+0/
      DATA VX(1)/0.9D+0/, VX(2)/-1.D-3/, VX(3)/1.D+1/, VX(4)/0.8D+0/,
     1     VX(5)/1.D+2/, VX(6)/0.8D+0/, VX(7)/1.D+2/, VX(8)/0.5D+0/,
     2     VX(9)/0.5D+0/, VX(10)/1.D+0/, VX(11)/1.D+0/, VX(14)/0.1D+0/,
     3     VX(15)/1.D+0/, VX(16)/1.D+0/, VX(19)/1.D+0/, VX(23)/1.D+0/,
     4     VX(24)/1.D+0/, VX(25)/1.D+0/, VX(26)/1.D+0/, VX(27)/1.D+10/,
     5     VX(29)/1.D+0/, VX(31)/1.D+0/, VX(32)/1.D+0/, VX(33)/1.D+0/,
     6     VX(34)/1.D+0/
C
C/6
      DATA VARNM(1)/1HP/, VARNM(2)/1HN/, SH(1)/1HS/, SH(2)/1HH/
      DATA CNGD(1),CNGD(2),CNGD(3)/4H---C,4HHANG,4HED V/,
     1     DFLT(1),DFLT(2),DFLT(3)/4HNOND,4HEFAU,4HLT V/
C/7
C     DATA VARNM(1)/'P'/, VARNM(2)/'N'/, SH(1)/'S'/, SH(2)/'H'/
C     DATA CNGD(1),CNGD(2),CNGD(3)/'---C','HANG','ED V'/,
C    1     DFLT(1),DFLT(2),DFLT(3)/'NOND','EFAU','LT V'/
C/
      DATA IJMP/33/, JLIM(1)/0/, JLIM(2)/24/, NDFLT(1)/32/, NDFLT(2)/25/
      DATA MINIV(1)/80/, MINIV(2)/59/
C
C...............................  BODY  ................................
C
      IF (ALG .LT. 1 .OR. ALG .GT. 2) GO TO 330
      IF (IV(1) .EQ. 0) CALL DDEFLT(ALG, IV, LIV, LV, V)
      PU = IV(PRUNIT)
      MIV1 = MINIV(ALG)
      IF (PERM .LE. LIV) MIV1 = MAX0(MIV1, IV(PERM) - 1)
      IF (IVNEED .LE. LIV) MIV2 = MIV1 + MAX0(IV(IVNEED), 0)
      IF (LASTIV .LE. LIV) IV(LASTIV) = MIV2
      IF (LIV .LT. MIV1) GO TO 290
      IV(IVNEED) = 0
      IV(LASTV) = MAX0(IV(VNEED), 0) + IV(LMAT) - 1
      IF (LIV .LT. MIV2) GO TO 290
      IF (LV .LT. IV(LASTV)) GO TO 310
      IV(VNEED) = 0
      IF (ALG .EQ. IV(ALGSAV)) GO TO 20
         IF (PU .NE. 0) WRITE(PU,10) ALG, IV(ALGSAV)
 10      FORMAT(/' THE FIRST PARAMETER TO DDEFLT SHOULD BE',I3,
     1          12H RATHER THAN,I3)
         IV(1) = 82
         GO TO 999
 20   IV1 = IV(1)
      IF (IV1 .LT. 12 .OR. IV1 .GT. 14) GO TO 50
         IF (N .GE. 1) GO TO 40
              IV(1) = 81
              IF (PU .EQ. 0) GO TO 999
              WRITE(PU,30) VARNM(ALG), N
 30           FORMAT(/8H /// BAD,A1,2H =,I5)
              GO TO 999
 40      IF (IV1 .NE. 14) IV(NEXTIV) = IV(PERM)
         IF (IV1 .NE. 14) IV(NEXTV) = IV(LMAT)
         IF (IV1 .EQ. 13) GO TO 999
         K = IV(PARSAV) - EPSLON
         CALL DVDFLT(ALG, LV-K, V(K+1))
         IV(DTYPE0) = 2 - ALG
         IV(OLDN) = N
         WHICH(1) = DFLT(1)
         WHICH(2) = DFLT(2)
         WHICH(3) = DFLT(3)
         GO TO 100
 50   IF (N .EQ. IV(OLDN)) GO TO 70
         IV(1) = 17
         IF (PU .EQ. 0) GO TO 999
         WRITE(PU,60) VARNM(ALG), IV(OLDN), N
 60      FORMAT(/5H /// ,1A1,14H CHANGED FROM ,I5,4H TO ,I5)
         GO TO 999
C
 70   IF (IV1 .LE. 11 .AND. IV1 .GE. 1) GO TO 90
         IV(1) = 80
         IF (PU .NE. 0) WRITE(PU,80) IV1
 80      FORMAT(/13H ///  IV(1) =,I5,28H SHOULD BE BETWEEN 0 AND 14.)
         GO TO 999
C
 90   WHICH(1) = CNGD(1)
      WHICH(2) = CNGD(2)
      WHICH(3) = CNGD(3)
C
 100  IF (IV1 .EQ. 14) IV1 = 12
      IF (BIG .GT. TINY) GO TO 110
         TINY = D1MACH(1)
         MACHEP = D1MACH(4)
         BIG = D1MACH(2)
         VM(12) = MACHEP
         VX(12) = BIG
         VM(13) = TINY
         VX(13) = BIG
         VM(14) = MACHEP
         VM(17) = TINY
         VX(17) = BIG
         VM(18) = TINY
         VX(18) = BIG
         VX(20) = BIG
         VX(21) = BIG
         VX(22) = BIG
         VM(24) = MACHEP
         VM(25) = MACHEP
         VM(26) = MACHEP
         VX(28) = DSQRT(D1MACH(2))*16.
         VM(29) = MACHEP
         VX(30) = BIG
         VM(33) = MACHEP
 110  M = 0
      I = 1
      J = JLIM(ALG)
      K = EPSLON
      NDFALT = NDFLT(ALG)
      DO 140 L = 1, NDFALT
         VK = V(K)
         IF (VK .GE. VM(I) .AND. VK .LE. VX(I)) GO TO 130
              M = K
              IF (PU .NE. 0) WRITE(PU,120) VN(1,I), VN(2,I), K, VK,
     1                                    VM(I), VX(I)
 120          FORMAT(/6H ///  ,2A4,5H.. V(,I2,3H) =,D11.3,7H SHOULD,
     1               11H BE BETWEEN,D11.3,4H AND,D11.3)
 130     K = K + 1
         I = I + 1
         IF (I .EQ. J) I = IJMP
 140     CONTINUE
C
      IF (IV(NVDFLT) .EQ. NDFALT) GO TO 160
         IV(1) = 51
         IF (PU .EQ. 0) GO TO 999
         WRITE(PU,150) IV(NVDFLT), NDFALT
 150     FORMAT(/13H IV(NVDFLT) =,I5,13H RATHER THAN ,I5)
         GO TO 999
 160  IF ((IV(DTYPE) .GT. 0 .OR. V(DINIT) .GT. ZERO) .AND. IV1 .EQ. 12)
     1                  GO TO 190
      DO 180 I = 1, N
         IF (D(I) .GT. ZERO) GO TO 180
              M = 18
              IF (PU .NE. 0) WRITE(PU,170) I, D(I)
 170     FORMAT(/8H ///  D(,I3,3H) =,D11.3,19H SHOULD BE POSITIVE)
 180     CONTINUE
 190  IF (M .EQ. 0) GO TO 200
         IV(1) = M
         GO TO 999
C
 200  IF (PU .EQ. 0 .OR. IV(PARPRT) .EQ. 0) GO TO 999
      IF (IV1 .NE. 12 .OR. IV(INITS) .EQ. ALG-1) GO TO 220
         M = 1
         WRITE(PU,210) SH(ALG), IV(INITS)
 210     FORMAT(/22H NONDEFAULT VALUES..../5H INIT,A1,14H..... IV(25) =,
     1          I3)
 220  IF (IV(DTYPE) .EQ. IV(DTYPE0)) GO TO 240
         IF (M .EQ. 0) WRITE(PU,250) WHICH
         M = 1
         WRITE(PU,230) IV(DTYPE)
 230     FORMAT(20H DTYPE..... IV(16) =,I3)
 240  I = 1
      J = JLIM(ALG)
      K = EPSLON
      L = IV(PARSAV)
      NDFALT = NDFLT(ALG)
      DO 280 II = 1, NDFALT
         IF (V(K) .EQ. V(L)) GO TO 270
              IF (M .EQ. 0) WRITE(PU,250) WHICH
 250          FORMAT(/1H ,3A4,9HALUES..../)
              M = 1
              WRITE(PU,260) VN(1,I), VN(2,I), K, V(K)
 260          FORMAT(1X,2A4,5H.. V(,I2,3H) =,D15.7)
 270     K = K + 1
         L = L + 1
         I = I + 1
         IF (I .EQ. J) I = IJMP
 280     CONTINUE
C
      IV(DTYPE0) = IV(DTYPE)
      PARSV1 = IV(PARSAV)
      CALL DCOPY(IV(NVDFLT), V(EPSLON),1,V(PARSV1),1)
      GO TO 999
C
 290  IV(1) = 15
      IF (PU .EQ. 0) GO TO 999
      WRITE(PU,300) LIV, MIV2
 300  FORMAT(/10H /// LIV =,I5,17H MUST BE AT LEAST,I5)
      IF (LIV .LT. MIV1) GO TO 999
      IF (LV .LT. IV(LASTV)) GO TO 310
      GO TO 999
C
 310  IV(1) = 16
      IF (PU .EQ. 0) GO TO 999
      WRITE(PU,320) LV, IV(LASTV)
 320  FORMAT(/9H /// LV =,I5,17H MUST BE AT LEAST,I5)
      GO TO 999
C
 330  IV(1) = 67
C
 999  RETURN
C  ***  LAST CARD OF DPARCK FOLLOWS  ***
      END
