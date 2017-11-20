      real*8 function vcirc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns circular velocity in mid-plane for a composite disc + halo model
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
      real*8 frtot, vbss
c
c local variables
      integer ip
      real*8 haloc, hrad, r2, v2
c
      if( r .eq. 0.d0 )then
        vcirc = 0.
      else
        if( fixed )then
          do ip = 1, ncmp
            if( .not. disc( ip ) )then
              haloc = dfcns( 3, ip )
              if( ctype( ip ) .eq. 'FE  ' )then
c Fall/Efstathiou model
                hrad = rscale( ip )
                r2 = r / hrad
                vcirc = haloc * r2 / sqrt( 1. + r2 * r2 )
              else if( ctype( ip ) .eq. 'BSS ' )then
c Bahcall/Schmidt/Soniera model
                vcirc = vbss( abs( r ) )
              else if( ctype( ip ) .eq. '3198' )then
c NGC 3198 model
                call crash( 'VCIRC', '3198 option not available' )
c                r2 = r * r
c                v2 = .63**2 * r2 / ( r2 + 4.7**2 ) + vdisc( r )**2 / dmass
c                vcirc = sqrt( v2 )
              else if( ctype( ip ) .eq. 'POWE' )then
c power law rotation curve
                vcirc = r**haloc
              else if( ctype( ip ) .eq. 'POWC' )then
c power law with a core
                hrad = rscale( ip )
                if( abs( haloc ) .lt. .01 )then
                  vcirc = r / sqrt( hrad * ( r**2 + hrad**2 ) )
                else
                  haloc = dfcns( 3, ip )
                  vcirc = r /
     +                   ( 1. +  ( r / hrad )**2 )**( .25 * haloc + .5 )
                  vcirc = vcirc * sqrt( abs( haloc ) / hrad**3 )
                end if
              end if
            end if
          end do
        else
c arbitrary model
          r2 = abs( r )
          v2 = max( -r2 * frtot( r2 ), 0.d0 )
          vcirc = sqrt( v2 )
        end if
      end if
      vcirc = sign( vcirc, r )
      return
      end
