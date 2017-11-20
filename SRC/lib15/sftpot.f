      real function sftpot( d2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns potential from a softened point mass at squared distance d2
c
c calling argument
      real d2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variable
      real x
c
      if( tsoft .eq. 1 )then
c Plummer softening
        sftpot = -1. / sqrt( d2 + softl2 )
      else if( tsoft .eq. 2 )then
c cubic spline softening
        x = sqrt( d2 ) / softl
        if( x .le. 0. )then
          sftpot = -1.4 / softl
        else if( x .lt. 1. )then
          sftpot = -1.4 + x**2 * ( 0.6666667 - .3 * x**2 + .1 * x**3 )
          sftpot = sftpot / softl
        else if( x .lt. 2. )then
          sftpot = -1.6 + 1. / ( 15. * x ) + x**2 *
     +               ( 1.33333333 - x + 0.3 * x**2 - .033333333 * x**3 )
          sftpot = sftpot / softl
        else
          sftpot = -1. / sqrt( d2 )
        end if
      else
        call crash( 'SFTPOT', 'Unrecognized softening rule' )
      end if
      return
      end
