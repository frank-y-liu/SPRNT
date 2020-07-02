      INTEGER FUNCTION ISMIN(N,SX,INCX)
C                                                                       
C     FINDS THE 1ST INDEX OF ELEMENT HAVING MIN. VALUE.                 
C     D. KAHANER,LASL,6/14/79                                           
C                                                                       
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
      REAL SX(*),SMIN                                                   
      INTEGER I,INCX,IX,N                                               
C                                                                       
      ISMIN = 0                                                         
      IF( N .LT. 1 ) RETURN                                             
      ISMIN = 1                                                         
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      IX = 1                                                            
      SMIN = (SX(1))                                                    
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF((SX(IX)).GE.SMIN) GO TO 5                                   
         ISMIN = I                                                      
         SMIN = (SX(IX))                                                
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 SMIN = (SX(1))                                                    
      DO 30 I = 2,N                                                     
         IF((SX(I)).GE.SMIN) GO TO 30                                   
         ISMIN = I                                                      
         SMIN = (SX(I))                                                 
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
