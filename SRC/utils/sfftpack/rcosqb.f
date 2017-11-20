      SUBROUTINE RCOSQB (N,X,WSAVE)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION       X(*)       ,WSAVE(*)
      DATA TSQRT2 /2.82842712474619009760/
      IF (N-2) 101,102,103
  101 X(1) = 4.0*X(1)
      RETURN
  102 X1 = 4.0*(X(1)+X(2))
      X(2) = TSQRT2*(X(1)-X(2))
      X(1) = X1
      RETURN
  103 CALL SCOSQB1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
