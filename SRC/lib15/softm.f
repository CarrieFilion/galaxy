      real function softm( ds )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns fraction of a softened point mass interior to scaled distance ds
c
c calling argument
      real ds
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      if( tsoft .eq. 2 )then
c cubic spline softening
        if( ds .le. 0. )then
          softm = 0
        else if( ds .lt. 1. )then
          softm = 4. * ds**3 *
     +                        ( .333333333 - .3 * ds**2 + .125 * ds**3 )
        else if( ds .lt. 2. )then
          softm = -1. / 15. + ds**3 *
     +        ( 2.6666667 - 3. * ds + 1.2 * ds**2 - .166666667 * ds**3 )
        else
          softm = 1
        end if
      else
        call crash( 'SOFTM', 'Unrecognized softening rule' )
      end if
      return
      end
