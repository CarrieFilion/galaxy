      subroutine euler( alpha, beta, gamma, a )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to set up a rotation matrix for the three Euler angles defined
c   in Arfken pp 199-200
c
c calling arguments
      real*8 alpha, beta, gamma, a( 3, 3 )
c
c local variables
      real*8 ca, cb, cg, sa, sb, sg
c
      ca = cos( alpha )
      sa = sin( alpha )
      cb = cos( beta )
      sb = sin( beta )
      cg = cos( gamma )
      sg = sin( gamma )
c construct rotation matrix - my convention is 1st index for the column
      a( 1, 1 ) = cg * cb * ca - sg * sa
      a( 1, 2 ) = -sg * cb * ca - cg * sa
      a( 1, 3 ) = sb * ca
      a( 2, 1 ) = cg * cb * sa + sg * ca
      a( 2, 2 ) = -sg * cb * sa + cg * ca
      a( 2, 3 ) = sb * sa
      a( 3, 1 ) = -cg * sb
      a( 3, 2 ) = sg * sb
      a( 3, 3 ) = cb
      return
      end
