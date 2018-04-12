      DOUBLE PRECISION FUNCTION DB3VAL(XVAL,YVAL,ZVAL,IDX,IDY,IDZ,
     *  TX,TY,TZ,NX,NY,NZ,KX,KY,KZ,BCOEF,WORK)
C***BEGIN PROLOGUE  DB3VAL
C***DATE WRITTEN   25 MAY 1982
C***REVISION DATE  25 MAY 1982
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  INTERPOLATION, THREE-DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  BOISVERT, RONALD, NBS
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             WASHINGTON, DC 20234
C***PURPOSE  DB3VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
C            FUNCTION CONSTRUCTED BY THE ROUTINE B3INK OR ONE OF ITS
C            PARTIAL DERIVATIVES.
C            DOUBLE PRECISION VERSION OF B3VAL.
C***DESCRIPTION
C
C   DB3VAL  evaluates   the   tensor   product   piecewise   polynomial
C   interpolant constructed  by  the  routine  DB3INK  or  one  of  its
C   derivatives  at  the  point  (XVAL,YVAL,ZVAL).  To   evaluate   the
C   interpolant  itself,  set  IDX=IDY=IDZ=0,  to  evaluate  the  first
C   partial with respect to x, set IDX=1,IDY=IDZ=0, and so on.
C
C   DB3VAL returns 0.0D0 if (XVAL,YVAL,ZVAL) is out of range. That is,
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
C            ZVAL.LT.TZ(1) .OR. ZVAL.GT.TZ(NZ+KZ)
C   If the knots TX, TY, and TZ were chosen by  DB3INK,  then  this  is
C   equivalent to
C            XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
C            YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY .OR.
C            ZVAL.LT.Z(1) .OR. ZVAL.GT.Z(NZ)+EPSZ
C   where EPSX = 0.1*(X(NX)-X(NX-1)), EPSY =  0.1*(Y(NY)-Y(NY-1)),  and
C   EPSZ = 0.1*(Z(NZ)-Z(NZ-1)).
C
C   The input quantities TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and  BCOEF
C   should remain unchanged since the last call of DB3INK.
C
C
C   I N P U T
C   ---------
C
C   XVAL    Double precision scalar
C           X coordinate of evaluation point.
C
C   YVAL    Double precision scalar
C           Y coordinate of evaluation point.
C
C   ZVAL    Double precision scalar
C           Z coordinate of evaluation point.
C
C   IDX     Integer scalar
C           X derivative of piecewise polynomial to evaluate.
C
C   IDY     Integer scalar
C           Y derivative of piecewise polynomial to evaluate.
C
C   IDZ     Integer scalar
C           Z derivative of piecewise polynomial to evaluate.
C
C   TX      Double precision 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Same as in last call to DB3INK.)
C
C   TY      Double precision 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Same as in last call to DB3INK.)
C
C   TZ      Double precision 1D array (size NZ+KZ)
C           Sequence of knots defining the piecewise polynomial in
C           the z direction.  (Same as in last call to DB3INK.)
C
C   NX      Integer scalar
C           The number of interpolation points in x.
C           (Same as in last call to DB3INK.)
C
C   NY      Integer scalar
C           The number of interpolation points in y.
C           (Same as in last call to DB3INK.)
C
C   NZ      Integer scalar
C           The number of interpolation points in z.
C           (Same as in last call to DB3INK.)
C
C   KX      Integer scalar
C           Order of polynomial pieces in x.
C           (Same as in last call to DB3INK.)
C
C   KY      Integer scalar
C           Order of polynomial pieces in y.
C           (Same as in last call to DB3INK.)
C
C   KZ      Integer scalar
C           Order of polynomial pieces in z.
C           (Same as in last call to DB3INK.)
C
C   BCOEF   Double precision 2D array (size NX by NY by NZ)
C           The B-spline coefficients computed by DB3INK.
C
C   WORK    Double precision 1D array (size KY*KZ+3*max(KX,KY,KZ)+KZ)
C           A working storage array.
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C***ROUTINES CALLED  DINTRV,DBVALU
C***END PROLOGUE  DB3VAL
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   MODIFICATION
C   ------------
C
C   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
C
C   R.F. BOISVERT, NIST
C   22 FEB 00
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        IDX, IDY, IDZ, NX, NY, NZ, KX, KY, KZ
      DOUBLE PRECISION
     *     XVAL, YVAL, ZVAL, TX(*), TY(*), TZ(*), BCOEF(NX,NY,NZ),
     *     WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        ILOY, ILOZ, INBVX, INBV1, INBV2, LEFTY, LEFTZ, MFLAG,
     *        KCOLY, KCOLZ, IZ, IZM1, IW, I, J, K
      DOUBLE PRECISION
     *     DBVALU
C
      DATA ILOY /1/,  ILOZ /1/,  INBVX /1/
C     SAVE ILOY    ,  ILOZ    ,  INBVX
C
C
C***FIRST EXECUTABLE STATEMENT
      DB3VAL = 0.0D0
C  NEXT STATEMENT - RFB MOD
      IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
     +    YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
     +    ZVAL.LT.TZ(1) .OR. ZVAL.GT.TZ(NZ+KZ)) GO TO 100
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0)  GO TO 100
      CALL DINTRV(TZ,NZ+KZ,ZVAL,ILOZ,LEFTZ,MFLAG)
      IF (MFLAG .NE. 0)  GO TO 100
         IZ = 1 + KY*KZ
         IW = IZ + KZ
         KCOLZ = LEFTZ - KZ
         I = 0
         DO 50 K=1,KZ
            KCOLZ = KCOLZ + 1
            KCOLY = LEFTY - KY
            DO 50 J=1,KY
               I = I + 1
               KCOLY = KCOLY + 1
               WORK(I) = DBVALU(TX,BCOEF(1,KCOLY,KCOLZ),NX,KX,IDX,XVAL,
     *                           INBVX,WORK(IW))
   50    CONTINUE
         INBV1 = 1
         IZM1 = IZ - 1
         KCOLY = LEFTY - KY + 1
         DO 60 K=1,KZ
            I = (K-1)*KY + 1
            J = IZM1 + K
            WORK(J) = DBVALU(TY(KCOLY),WORK(I),KY,KY,IDY,YVAL,
     *                           INBV1,WORK(IW))
  60     CONTINUE
         INBV2 = 1
         KCOLZ = LEFTZ - KZ + 1
         DB3VAL = DBVALU(TZ(KCOLZ),WORK(IZ),KZ,KZ,IDZ,ZVAL,INBV2,
     *                  WORK(IW))
  100 CONTINUE
      RETURN
      END
