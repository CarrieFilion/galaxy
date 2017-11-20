      REAL*8 FUNCTION SGULIA(R)
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
C RETURNS RADIAL VELOCITY DISPERSION FOR LIA'S DIST. FN. FOR KT DISC
C
C CALLING ARGUMENT
      REAL*8 R
C
C COMMON BLOCK
c
      include 'inc/params.f'
c
      include 'inc/model.f'
C
C LOCAL VARIABLES
      INTEGER K,L
      REAL*8 FK,FK2,FL,FL2,FMH,PHI,PHI2,RR,R1,SUMK,SUML,T,TK,TL
      REAL*8 TL1,TL2,X2,X3
      real*8 dfb, dfm
C
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
      RR=R*R
      PHI2=1./(1.+RR)
      PHI=SQRT(PHI2)
      X2=PHI2*RR
      X3=-PHI2*(1.+RR)*dfb
      FMH=dfm/2.
      K=0
      SUML=0.
      SUMK=0.
      FK=K
      FK2=2*K
      TK=1.
    1 L=0
      TL=1.
      TL1=1./(FK2+dfm+1.)
      SUML=TL1
    2 L=L+1
      FL=L
      FL2=2*L
      TL=TL*X3/FL
      TL2=TL/(FK2+FL2+1.+dfm)
      SUML=SUML+TL2
      IF(ABS(TL2).GT.1.E-3*ABS(SUML)) GO TO 2
      T=TK*SUML
      SUMK=SUMK+T
      IF(ABS(T).LE.1.E-3*ABS(SUMK)) GO TO 3
      K=K+1
      FK=K
      FK2=2*K
      TK=TK*(FK+FMH-2.5)/FK
      TK=TK*X2
      GO TO 1
    3 SGULIA=PHI**(dfm+1.)*EXP(dfb)*SUMK
      R1=1.+RR
      SGULIA=SGULIA*R1**1.5
      SGULIA=SQRT(SGULIA)
      RETURN
      END
