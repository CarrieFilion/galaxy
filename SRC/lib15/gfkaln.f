      real*8 function gfkaln( eta )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the value of Kalnajs's function g for argument eta
c   part of the determination of the DF for a disk
c
c calling argument
      real*8 eta
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      common / kalndf / xx, m
      integer m
      real*8 xx
c
c external
      real*8 ymtau
c
c local variables
      integer n
      real*8 dfm, pa, pb, pc
c
      dfm = dfcns( 1, icmp )
      gfkaln = 0
      if( m .gt. 2 )then
c second derivative of Legendre polynomial of degree m-1
        if( eta .eq. 1. )then
          gfkaln = .125 * dble( ( m + 1 ) * m * ( m - 1 ) * ( m - 2 ) )
        else
c evaluate required Legendre polynomials using recurrence relation
          pa = 1
          pb = eta
          n = 1
          do while ( n .lt. m - 1 )
            pc = ( dble( 2 * n + 1 ) * eta * pb - dble( n ) * pa )
     +                                                  / dble( n + 1 )
            n = n + 1
            if( n .lt. m - 1 )then
              pa = pb
              pb = pc
            end if
          end do
          gfkaln = pc * ( ( dfm - 2. ) * eta * eta - dfm ) +
     +                                                  2. * eta * pb
          gfkaln = gfkaln * ( dfm - 1. ) / ( eta * eta - 1. )**2
        end if
c second factor
        gfkaln = gfkaln * ymtau( xx * eta )
      end if
      return
      end
