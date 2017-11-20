      SUBROUTINE RCOSQI (N,WSAVE)
      IMPLICIT real*4 (A-H,O-Z)
      DIMENSION       WSAVE(1)
      DATA PIH /1.57079632679489661923/
      DT = PIH/FLOAT(N)
      FK = 0.0D0
      DO 101 K=1,N
         FK = FK+1.0
         WSAVE(K) = COS(FK*DT)
  101 CONTINUE
      CALL SFFTI (N,WSAVE(N+1))
      RETURN
      END