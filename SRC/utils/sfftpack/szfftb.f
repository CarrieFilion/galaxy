      SUBROUTINE SZFFTB (N,C,WSAVE)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION       C(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL SCFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
