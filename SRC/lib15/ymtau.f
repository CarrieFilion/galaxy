      real*8 function ymtau( y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Evaluates y**m * tau(y) needed for Kalnajs type DFs.  The expression
c   can be recast as r**m * sigma( r ) - see Kalnajs Ap.J. 205, p 755
c
c calling argument
      real*8 y
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gsigma, rofy
c
c local variables
      real*8 dfm, r, y1
      include 'inc/pi.f'
c
      ymtau = 0
      if( y .gt. 0. )then
        if( ( ncmp .gt. 1 ) .or. ( .not. disc( 1 ) )
     +        )call crash( 'YMTAU', 'called for a not-disc only model' )
        dfm = dfcns( 1, icmp )
c Kuz'min/Toomre disc
        if( ctype( 1 ) .eq. 'KT  ' )then
          y1 = 1. - 1.e-12
          y1 = min( y1, y )
          ymtau = y1**dfm / ( 1. - y1 * y1 )**( .5 * ( dfm - 3. ) )
          ymtau = ymtau / ( 2. * pi )
        else
c general expression
          r = rofy( y )
          ymtau = ( r / rstar )**dfm * gsigma( r )
        end if
      end if
      return
      end
