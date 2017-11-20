      real*8 function vbss( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the circular speed at radius R in the Bahcall-Schmidt-Soniera
c   model for the Milky Way
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local arrays
      integer m
      parameter ( m = 22 )
      real rad( m ), vdisc( m ), vcen( m ), vsph( m ), vhalo( m )
      real*8 c( m + 4 ), lamda( m + 4 ), x( m ), y( m )
c
c local variables
      integer i, iuse
      real*8 rs
      save c, iuse, lamda
c
      data iuse / 0 /
c
c rotation curve data from Bahcall, Shcmidt & Soniera paper
      data rad /   0.0,   0.25,  0.5,   0.75,  1.0,   2.0,   3.0,   4.0,
     +             5.0,   6.0,   7.0,   7.7,   8.0,   8.3,   9.0,  10.0,
     +            12.0,  14.0,  16.0,  18.0,  20.0,  25.0 /
      data vdisc / 0.0,  22.0,  38.5,  52.5,  64.7, 101.2, 124.9, 140.2,
     +           149.6, 154.9, 157.2, 157.4, 157.2, 156.9, 155.8, 153.2,
     +           146.1, 138.0, 129.7, 121.9, 114.8, 100.4 /
      data vcen /  0.0, 231.3, 243.2, 241.0, 227.8, 159.7, 128.3, 110.6,
     +            98.7,  90.0,  83.3,  79.3,  77.8,  76.4,  73.4,  69.6,
     +            63.5,  58.8,  55.0,  51.8,  49.2,  44.0 /
      data vsph /  0.0,  41.9,  45.4,  46.2,  46.2,  44.0,  41.4,  39.1,
     +            37.1,  35.3,  33.8,  32.8,  32.4,  32.0,  31.2,  30.1,
     +            28.2,  26.6,  25.2,  24.0,  23.0,  20.9 /
      data vhalo / 0.0,   7.7,  15.0,  21.7,  27.9,  49.1,  66.1,  80.2,
     +            92.4, 103.1, 112.8, 119.0, 121.6, 124.0, 129.6, 137.1,
     +           150.7, 162.8, 173.8, 184.0, 193.4, 214.4 /
c
      if( iuse .eq. 0 )then
c initialisation
        do i = 1, m
          x( i ) = rad( i ) / 3.5
          y( i ) = sqrt( vdisc( i )**2 + vcen( i )**2 + vsph( i )**2
     +                                                + vhalo( i )**2 )
          y( i ) = y( i ) / 262.1
        end do
c fit a cubic spline
        call splint2( x, y, m, 1.d0, lamda, c, .true. )
        iuse = 1
      end if
c interpolate value at requested point - flat rotation curve beyond rmax
      rs = min( r, x( m ) )
      call splint2( x, y, m, rs, lamda, c, .false. )
      return
      end
