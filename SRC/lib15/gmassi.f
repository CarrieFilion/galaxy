      real*8 function gmassi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns disc mass interior to r by integrating over DF
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
      include 'inc/model.f'
c
c externals
      real*8 gmfuni, gmassd
c
c local arrays
      integer npts
      parameter ( npts = 32 )
      real*8 abs( npts ), wt( npts )
c
c local variables
      integer i
      real*8 r1, r2
c
      gmassi = 0
      if( ( r .gt. rhole ) .and. disc( icmp ) )then
c M(R) is unknown when density is tapered
        if( Lztapr( icmp ) .or.
     +      ( dist( icmp ) .and. ( .not. kcold ) ) )then
c choose weights and abscissae for Gauss-Legengre quadrature
          r1 = rhole
          r2 = min( r, rmax )
          call GLqtab( r1, r2, npts, wt, abs )
c sum contributions
          gmassi = 0
          do i = 1, npts
            r1 = abs( i )
            gmassi = gmassi + wt( i ) * gmfuni( r1 )
          end do
        else
c use known function where possible
          gmassi = gmassd( r )
        end if
      end if
      return
      end
