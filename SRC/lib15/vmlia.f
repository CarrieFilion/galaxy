      REAL*8 FUNCTION VMLIA(R)
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
C RETURNS VALUE OF MEAN ORBITAL VELOCITY FOR LIA'S DIST. FN. IN K/T DISC
C
C CALLING ARGUMENT
      REAL*8 R
C
C COMMON BLOCK
c
      include 'inc/params.f'
C
      include 'inc/model.f'
C
C EXTERNAL
      REAL*8 GAMMAF
C
C LOCAL VARIABLES
      INTEGER IFAIL,J,JFAIL,K,KFAIL,L
      REAL*8 EBETA,FJ,FJ2,FK,FK2,FL,FL2,FMHH,FMJK2,FMJ2,FMK2,FMK21
      REAL*8 R2,SUMJ,SUMK,SUML,T,TJ,TK,TL,TL1,TT,W,W2,Y2
      real*8 dfb, dfm
C
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
      R2=R*R
      Y2=R2/(1.+R2)
      W=1./SQRT(1.+R2)
      W2=W*W
      J=0
      FJ=0.
      FJ2=0.
      FMJ2=dfm
      FMHH=.5*(dfm-5.)
      EBETA=EXP(dfb)
      SUMJ=0.
      IFAIL=0
      JFAIL=0
      KFAIL=0
      TL1=GAMMAF(dfm+1.,IFAIL)/
     +               (GAMMAF(dfm+1.5,JFAIL)*GAMMAF(.5d0,KFAIL))
      TJ=1.
    1 K=0
      FK=0.
      FK2=0.
      FMK2=dfm
      FMK21=FMK2-1.
      FMJK2=FMJ2
      TK=1.
      SUMK=0.
    2 L=0
      TL=TL1
      SUML=TL1
    3 L=L+1
      FL=L
      FL2=2*L
      TL=TL*(FL+FMHH)/FL
      T=FMJK2+FL2
      TT=FL+FJ
      TL=TL*TT/(T+.5)
      TL=TL*(T-1.)/(TT-0.5)
      TL=TL*T/(T-.5)
      TL=TL*Y2
      SUML=SUML+TL
      IF(ABS(TL).GT.2.E-5*ABS(SUML))GO TO 3
      T=TK*SUML
      SUMK=SUMK+T
      IF(ABS(T).LE.4.E-5*ABS(SUMK))GO TO 4
      K=K+1
      FK=K
      FK2=2*K
      TK=-TK*dfb*W2/FK
      FMK2=dfm+FK2
      FMK21=FMK2-1.
      FMJK2=FMJ2+FK2
      TK=TK*(FMJK2-1.)/(FMJK2-.5)
      TK=TK*FMJK2/(FMJK2+.5)
      GO TO 2
    4 T=TJ*SUMK
      SUMJ=SUMJ+T
      IF(ABS(T).LE.1.E-4*ABS(SUMJ))GO TO 5
      J=J+1
      FJ=J
      FJ2=2*J
      FMJ2=dfm+FJ2
      TJ=-TJ*dfb*Y2/(FMJ2+.5)
      TJ=TJ*(FMJ2-1.)/(FJ-.5)
      TJ=TJ*(FMJ2)/(FMJ2-.5)
      GO TO 1
    5 VMLIA=EBETA*W**dfm*SUMJ
      VMLIA=VMLIA*SQRT(2.)*(1.+R2)**1.5*SQRT(W)
      RETURN
      END
