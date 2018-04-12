      SUBROUTINE DDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,Y,YWT,H,MNTOLD,MTROLD,NFE,RC,YH,
     8   A,CONVRG,EL,IER,IPVT,NQ,NWAIT,RH,RMAX,SAVE2,TQ,TREND,ISWFLG)
C***BEGIN PROLOGUE  DDNTL
C***REFER TO  DDRIV3
C  Subroutine DDNTL is called to set parameters on the first call
C  to DDSTP, on an internal restart, or when the user has altered
C  MINT, MITER, and/or H.
C  On the first call, the order is set to 1 and the initial derivatives
C  are calculated.  RMAX is the maximum ratio by which H can be
C  increased in one step.  It is initially RMINIT to compensate
C  for the small initial H, but then is normally equal to RMNORM.
C  If a failure occurs (in corrector convergence or error test), RMAX
C  is set at RMFAIL for the next increase.
C  If the caller has changed MINT, or if JTASK = 0, DDCST is called
C  to set the coefficients of the method.  If the caller has changed H,
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is
C  reset to NQ + 2 to prevent further increases in H for that many
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is
C  set to 0 to force the partials to be updated, if partials are used.
C***ROUTINES CALLED  DDCST,SDSCL,DGEFA,DGESL,DGBFA,DGBSL,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDNTL
      DOUBLE PRECISION A(MATDIM,*), EL(13,12), EPS, H, HMAX, HNEW, HOLD,
     8     OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX, SMIN,
     8     DNRM2, SUM, SUM0, T, TQ(3,12), TREND, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.D0)
C***FIRST EXECUTABLE STATEMENT  DDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
        RC = 0.D0
        CONVRG = .FALSE.
        TREND = 1.D0
        RMAX = RMINIT
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              CALL DGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              CALL DGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            DO 150 I = 1,NDE
              IF(A(I,1) .EQ. 0.D0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = 0.D0
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/YWT(I)
        SUM = DNRM2(NDE, SAVE1, 1)
        SUM0 = 1.D0/MAX(1.D0, ABS(T))
        SMAX = MAX(SUM0, SUM)
        SMIN = MIN(SUM0, SUM)
        SUM = SMAX*SQRT(1.D0 + (SMIN/SMAX)**2)/SQRT(DBLE(NDE))
        H = SIGN(MIN(2.D0*EPS/SUM, ABS(H)), H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.D0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          HNEW = H
          RH = H/HOLD
          H = HOLD
          CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          H = SIGN(MIN(ABS(H), ABS(HNEW)), H)
        END IF
      END IF
      END
