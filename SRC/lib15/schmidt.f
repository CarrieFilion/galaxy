      real*8 function schmidt( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the shperical density of Schmidt's model for the Milky Way (I think)
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
c external
      real*8 algrng2
c
c local arrays
      real*8 an( 15 ), rhof( 15 )
c
c local variables
      integer iuse
      real*8 facrho, haloc
      save iuse, facrho
c
      data iuse / 0 /
      data an /  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,
     +           2.2,  2.3,  2.4,  2.5,  2.6,  2.7 /
      data rhof / 126.89528,  64.99177,  44.77646,  34.89432,  29.10037,
     +             25.32542,  22.68915,  20.75503,  19.28271,  18.12926,
     +             17.20464,  16.44943,  15.82259,  15.29553,  14.84720/
c
      haloc = dfcns( 3, icmp )
      if( iuse .eq. 0 )then
        facrho = algrng2( an, rhof, 15, haloc )
        iuse = 1
      end if
      schmidt = 1. / ( facrho * r**1.8 * ( 1. + r**haloc ) )
      return
      end
