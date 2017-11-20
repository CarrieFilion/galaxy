      real*8 function dfmiym( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c DF given by Miyamoto (1971 PASJ v23 p21) for the Kuz'min-Toomre disk
c   given by the sum of two terminating power series for the even and odd
c   parts separately
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      integer is
      real*8 E2, feven, fodd, Lz2, t
c
      E2 = -2. * E
      if( Lz .eq. 0.d0 )then
c zero angular momentum value
        fodd = 0.
        feven = amiy( 1, mmiy + 1 ) * E2**( 2 * mmiy + 2 )
      else
c even part
        feven = 0.
        t = E2**( 2 + mmiy ) * Lz**( 2 * mmiy )
        Lz2 = Lz * Lz
        do is = 0, mmiy
          feven = feven + amiy( 1, is + 1 ) * t
          t = t * E2 / Lz2
        end do
c odd part
        fodd = 0.
        t = E2**( real( nmiy ) + 2.5 ) * Lz**( 2 * nmiy + 1 )
        do is = 0, nmiy
          fodd = fodd + amiy( 2, is + 1 ) * t
          t = t * E2 / Lz2
        end do
      end if
      dfmiym = feven + fodd
      return
      end
