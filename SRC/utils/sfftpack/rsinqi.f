      SUBROUTINE RSINQI (N,WSAVE)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION       WSAVE(1)
      CALL RCOSQI (N,WSAVE)
      RETURN
      END
