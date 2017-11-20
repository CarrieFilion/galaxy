      real function sbulge( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the projected surface density of the currently selected "halo"
c
c calling argument
      real r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real hrad
      include 'inc/pi.f'
c
      if( .not. disc( icmp ) )then
c Plummer sphere
        if( ctype( icmp ) .eq. 'PLUM' )then
          hrad = rscale( icmp )
          sbulge = cmpmas( icmp ) * hrad * hrad /
     +                               ( pi * ( r * r + hrad * hrad )**2 )
        end if
      else
        call crash( 'SBULGE', 'Unrecognised halo' )
      end if
      return
      end
