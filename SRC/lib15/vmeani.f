      real*8 function vmeani( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes vmean from a double integral over distribution fn
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c external
      real*8 velint
c
c mean speed must be zero at the center, unless the disc has no retro ptcls
      vmeani = 0
      if( ( r .gt. 0. ) .or. ( Lzrmn( icmp ) .eq. 0. ) )then
c set power factors for velocity weighting
        iu = 0
        iv = 1
        vmeani = velint( r )
c normalize
        if( vmeani .gt. 0. )then
          iu = 0
          iv = 0
          vmeani = vmeani / velint( r )
        end if
      end if
      return
      end
