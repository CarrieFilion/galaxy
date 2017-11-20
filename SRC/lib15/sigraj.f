      real*8 function sigraj( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/bdpps.f'
c
c local arrays
      real*8 poly( 25, 25 )
      integer iuse( 25 )
c
c local variables
      integer i, j, j1, n, naj1
      real*8 a, b, c, rm, r12
      include 'inc/pi.f'
c
      data iuse / 25 * -1 /
c
      naj1 = naj + 1
      if( iuse( naj1 ) .ne. kaj )then
c evaluate polynomial coefficients
        c = 2 * ( naj + kaj ) + maj
        c = c + .5
        a = dble( maj + naj ) - .5
        b = naj
        n = 2 * kaj
        do i = 1, n
          a = a + 1.
          b = b + 1.
          c = c * a * b / dble( i**2 )
        end do
        c = 2. * sqrt( c / pi )
        a = -.5
        do i = 1, kaj
          a = a + 1.
          c = c * dble( i ) / a
        end do
c normalisation for azimuthal integral (ignored in Agris's paper ii)
      if( maj .eq. 0 )then
        c = c / ( 2. * pi )
      else
        c = c / ( sqrt( 2. ) * pi )
      end if
c j = 0 term
        poly( 1, naj1 ) = c
c general terms
        do j = 1, naj
          j1 = j - 1
          poly( j + 1, naj1 ) = poly( j, naj1 ) * dble( j1 - naj )
     +     * ( dble( naj + 2 * kaj + maj + j1 ) + .5 ) * dble( kaj + j )
     +     / ( ( dble( kaj + j1 ) + .5 ) * dble( j * ( 2 * kaj + j ) ) )
        end do
c reset flag
        iuse( naj1 ) = kaj
      end if
c polynomial in ( 1 - r**2 )
      sigraj = 0.
      r12 = 1. - r * r
      rm = 0.
      if( maj .eq. 0 )rm = 1
      if( r .gt. 0. )rm = r**maj * r12**( dble( kaj ) - .5 )
c sum over j from 0 to n
      do j = 1, naj1
        sigraj = sigraj + poly( j, naj1 ) * rm
        rm = rm * r12
      end do
c fudge that is not in Agris's eq ( 23 )!
      a = ( -1 )**naj
      sigraj = a * sigraj
      return
      end
