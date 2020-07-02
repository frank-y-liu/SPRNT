      INTEGER FUNCTION ISMAX(N,SX,INCX)
C                                                                       
C     FINDS THE 1ST INDEX OF ELEMENT HAVING MAX. VALUE.                 
C     D. KAHANER LASL 6/14/79                                           
C                                                                       
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
      REAL SX(*),SMAX                                                   
      INTEGER I,INCX,IX,N                                               
C                                                                       
      ISMAX = 0                                                         
      IF( N .LT. 1 ) RETURN                                             
      ISMAX = 1                                                         
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      IX = 1                                                            
      SMAX = (SX(1))                                                    
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF((SX(IX)).LE.SMAX) GO TO 5                                   
         ISMAX = I                                                      
         SMAX = (SX(IX))                                                
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 SMAX = (SX(1))                                                    
      DO 30 I = 2,N                                                     
         IF((SX(I)).LE.SMAX) GO TO 30                                   
         ISMAX = I                                                      
         SMAX = (SX(I))                                                 
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
