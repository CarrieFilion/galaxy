      real*8 function Lzmin( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c except for Kalnajs functions for the Maclaurin disk, this function returns
c   the specific minimum angular momentum a particle could have to prevent
c   any part of its orbit from lying within an inner hole in the mass
c   distribution.  rhole is already defined in the common block and if
c   this radius is zero, the value returned is, of course, zero.
c
c For Kalnajs function, the value returned is the most negative angular
c   momentum present for the selected value of Omega.
c
c calling arguments
      real*8 E
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 Phitot
c
c local variables
      real*8 bdy, v2
      include 'inc/pi.f'
c
c allow leeway for precision
      if( E - Emine .lt. -1.d-5 )then
        print *, E, Emine, E - Emine
        call crash( 'LZMIN', 'Impossible energy' )
      end if
c Omega models
      if( cdft( icmp ) .eq. 'OMEG' )then
        Lzmin = -sqrt( 4. / ( 3. * pi ) ) * ( E - Phi0 )
        if( Omegak .gt. 0. )then
          bdy = ( E - Phi0 - .5 * Omeg21 ) / Omegak
          Lzmin = max( Lzmin, bdy )
        end if
      else if( cdft( icmp ) .eq. 'USPS' )then
c Polyachenko-Shukhman DF for uniform sphere
        Lzmin = 0.d0
        if( E - Phimax .gt. 0.d0 )Lzmin =
     +                         rscale( 1 ) * sqrt( 2. * ( E - Phimax ) )
      else
c finds minimum angular momentum that will just avoid central hole
        Lzmin = 0.
        if( rhole .gt. 0. )then
          Lzmin = 2. * ( E - Phitot( rhole ) )
          Lzmin = rhole * sqrt( Lzmin )
        end if
c if E is near to Emax, Lzmin is further restricted
        v2 = 2. * ( E - Phimax )
        if( v2 .gt. 0. )Lzmin = max( Lzmin, rmax * sqrt( v2 ) )
      end if
      return
      end
