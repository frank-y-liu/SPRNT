      INTEGER FUNCTION ISAMIN(N,SX,INCX)
C                                                                       
C     FINDS THE 1ST INDEX OF ELEMENT HAVING MIN. ABSOLUTE VALUE.        
C     D. KAHANER,LASL,6/14/79                                           
C                                                                       
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
      REAL SX(*),SMIN                                                   
      INTEGER I,INCX,IX,N                                               
C                                                                       
      ISAMIN = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      ISAMIN = 1                                                        
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      IX = 1                                                            
      SMIN = ABS(SX(1))                                                 
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF(ABS(SX(IX)).GE.SMIN) GO TO 5                                
         ISAMIN = I                                                     
         SMIN = ABS(SX(IX))                                             
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 SMIN = ABS(SX(1))                                                 
      DO 30 I = 2,N                                                     
         IF(ABS(SX(I)).GE.SMIN) GO TO 30                                
         ISAMIN = I                                                     
         SMIN = ABS(SX(I))                                              
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
