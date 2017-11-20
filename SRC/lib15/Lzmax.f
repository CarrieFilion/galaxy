      real*8 function Lzmax( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the maximum allowed angular momentum (of a circuar orbit) of
c   a particle of energy E in the current theoretical potential well
c
c calling argument
      real*8 E
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 rcofe, vcirc
c
c local variables
      real*8 cc, r, r2
      include 'inc/pi.f'
c
      Lzmax = -1
c use analytic expression when possible
      if( ncmp .eq. 1 )then
        if( disc( 1 )  )then
c Maclaurin/Freeman/Kalnajs disc - r2 is square of rad of circ orb
          if( ctype( 1 ) .eq. 'MFK ' )then
            r2 = 1. + 4. * E / ( 3. * pi )
            Lzmax = sqrt( 3. * pi / 4. ) * r2
c Mestel/Toomre/Zang disc - circular orbit given by E = .5 + ln( Lzmax )
          else if( ctype( 1 ) .eq. 'MTZ ' )then
            Lzmax = exp( E - .5d0 )
c power law disc
          else if( ctype( 1 ) .eq. 'POWE' )then
            cc = dfcns( 3, 1 )
            r2 = ( 1.d0 + cc ) / ( 2. * cc )
            Lzmax = ( 2.d0 * E / ( 1.d0 + 1.d0 / cc ) )**r2
          end if
        else
          if( ctype( 1 ) .eq. 'KEPL' )then
c Kepler potential
            Lzmax = cmpmas( 1 ) / sqrt( -2. * E )
          else if( ctype( 1 ) .eq. 'UNIS' )then
c uniform sphere - PS DF extends to circular orbits at the outer edge
            if( E .le. .5 * Phimax )then
              Lzmax = ( E - Phi0 ) * sqrt( rscale( 1 ) / cmpmas( 1 ) )
            else
              Lzmax = cmpmas( 1 ) / sqrt( -2. * E )
            end if
          end if
        end if
      end if
c none of the above - find radius of circular orbit with this energy
      if( Lzmax .lt. 0.d0 )then
        r = rcofe( E )
        r = max( r, rhole )
        r = min( r, rmax )
        Lzmax = r * vcirc( r )
      end if
      return
      end
