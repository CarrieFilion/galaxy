      SUBROUTINE RZFFTI (N,WSAVE)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION       WSAVE(1)
      IF (N .EQ. 1) RETURN
      CALL SEZFFT1 (N,WSAVE(2*N+1),WSAVE(3*N+1))
      RETURN
      END
