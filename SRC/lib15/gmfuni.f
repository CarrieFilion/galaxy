      real*8 function gmfuni( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for mass as a function of r in the disc component.  Uses the
c   integrated mass from the DF (through tabulated values for greater
c   efficiency) where this is available, and otherwise substitutes the
c   tapered surface density function for cold discs
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gsigma, gsigmt, taper, vcirc
c
c local variables
      real*8 h
      include 'inc/pi.f'
c
      gmfuni = 0
      if( ( r .gt. rhole ) .and. ( disc( icmp ) ) )then
        if( ( .not. dist( icmp ) ) .or. kcold )then
          gmfuni = 2. * pi * r * gsigma( r )
          if( Lztapr( icmp ) )then
            h = r * vcirc( r )
            gmfuni = gmfuni * taper( h )
          end if
        else
          gmfuni = 2. * pi * r * gsigmt( r )
        end if
      end if
      return
      end
