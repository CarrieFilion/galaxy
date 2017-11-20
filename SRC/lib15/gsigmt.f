      real*8 function GSigmt( r )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c tabulated version of GSigmi - the disc surface density
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      real*8 algrng2, GSigma, GSigmi, taper
c
c local arrays
      integer ntab
      parameter ( ntab = 101 )
      real*8 x( ntab + 1 ), y( ntab + 1 )
      save x, y
c
c local variables
      integer iuse, j
      real*8 Lz, rad, radm
      save iuse, radm
      data iuse / 0 /
c
      GSigmt = 0
      if( r .le. rmax )then
c use analytic disc surface density if no DF type is defined
        if( .not. dist( icmp ) )then
          GSigmt = GSigma( r )
          if( Lztapr( icmp ) )GSigmt = GSigmt * taper( r )
        else
c create table if needed
          if( iuse .ne. icmp )then
            if( master )print *,
     +                        'Building a table of disc surface density'
            radm = rtrunc( icmp )
            Lz = 0
            if( singlr )Lz = 1.d-3
            do j = 1, ntab
             rad = Lz + ( radm - Lz ) * dble( j - 1 ) / dble( ntab - 1 )
              x( j ) = rad
              y( j ) = GSigmi( rad )
            end do
            if( master )print *, 'Surface density table ready'
            iuse = icmp
          end if
c look up value in table
          GSigmt = 0
          if( r .gt. x( 1 ) )then
            GSigmt = algrng2( x, y, ntab, r )
          else if( r .lt. radm )then
            GSigmt = y( 2 ) + ( y( 1 ) - y( 2 ) ) * ( x( 2 ) - r ) /
     +                                            ( x( 2 ) - x( 1 ) )
          end if
        end if
      end if
      return
      end
