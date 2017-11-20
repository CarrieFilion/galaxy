      subroutine rotate( alpha, beta, gamma, x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to rotate the 6 phase space coordinates of a particle through the
c   Euler angles defined in Arfken pp 199-200
c
c calling arguments
      real x( 6 )
      real*8 alpha, beta, gamma
c
c local array
      real*8 e( 3, 3 )
c
c set matrix
      call euler( alpha, beta, gamma, e )
c rotate vectors
      call rotmat( x, e )
      return
      end
