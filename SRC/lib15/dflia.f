            REAL*8 FUNCTION DFLIA(E,H)
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c unreconstructed routine written by Lia for 2-parameter DFs for the KT disk
C
C CALLING ARGUMENTS
      REAL*8 E,H
C
C COMMON BLOCK
c
      include 'inc/params.f'
C
      include 'inc/model.f'
C
C LOCAL VARIABLES
      INTEGER I,J,K
      REAL*8 BI,BIN,E2,FEK,FI,FI2,FK2,FMH,FXI,GI,SUMI,SUMK
      REAL*8 T,TERM,TK,TTK,TOI,X,X2
      real*8 dfb, dfm
      include 'inc/pi.f'
C
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
C COMPUTE X
      X=-H*SQRT(-2.*E) / ( phi0 * rstar )
C USEFUL FACTORS
      FMH=.5*dfm
      E2=E*E
      X2=X*X
C START FROM K=0
      K=0
      SUMK=0.
      FEK=1.
      TK=1.
      TTK=1.
C NEW I SUM FOR EACH K
    2 I=0
      SUMI=TTK
      BI=TTK
      TOI=1.
      FXI=1.
C NEXT I
    3 I=I+1
      FI=I
C TOI = BINOM((3/2-M/2),I)
      TOI=TOI*(2.5-FMH-FI)/FI
C GI = SUM( dfb**J/J! * BINOM((3/2-M/2),I-J) )
      GI=0.
      TERM=1.
      BIN=TOI
      DO 12 J=1,I
      GI=GI+BIN*TERM
      TERM=TERM*dfb/REAL(J)
      BIN=BIN*REAL(I-(J-1))/(2.5-FMH-REAL(I-(J-1)))
   12 CONTINUE
      GI=GI+TERM
C BI = (-.25)**I * GAMMA(2K+2I+M+1)/(GAMMA(I+.5)*GAMMA(2K+M+I))
      FI2=2*I
      FK2=2*K
      BI=-0.25*BI*((FK2+FI2+dfm)/(FI-.5))
      BI=BI*((FK2+FI2+dfm-1.)/(FK2+FI+dfm-1.))
C NEXT I TERM NOW READY
      FXI=FXI*X2
      T=GI*BI*FXI
      SUMI=SUMI+T
C I SUM CONVERGED?
      IF(ABS(T).GT.1.E-8*ABS(SUMI))GO TO 3
      T=FEK*SUMI*TK
      SUMK=SUMK+T
C K SUM CONVERGED?
      IF(ABS(T).LE.1.E-8*ABS(SUMK))GO TO 4
C NEXT K
      K=K+1
      FEK=FEK*E2
      TK=-dfb*TK/REAL(K)
      FK2=2*K
      TTK=(FK2+dfm)*TTK/(FK2+dfm-2.)
      GO TO 2
C OTHER FACTORS
    4 DFLIA=dfm*EXP(dfb)*(E/PHI0)**(dfm-1.)*SUMK/(2.*PI**2)
      RETURN
      END
