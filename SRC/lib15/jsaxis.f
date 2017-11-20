      subroutine jsaxis( aa, label, iflag )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw an axis, select and label tickmarks, and output the
c   given axis label
c   no label is written if the label flag is blank
c   if iflag=0, the tickmark values are not output, otherwise they are
c
c calling arguments
      character*1 aa
      character*(*) label
      integer iflag
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local variables
      character*1 a
      integer i, nchr, ndec, nw, nwmx
      logical integ, even
      real dx, mheight, rng, tickx, ticky, tiny, x, xa, xx, y, ya
c
      a = aa
      call lowercase( a )
c choose suitable spacing for tickmarks
      rng = sx2 - sx1
      if( a .eq. 'y' )rng = sy2 - sy1
      tiny = 1.e-3 * rng
      i = log10( .2 * rng )
      if( .2 * rng .lt. 1. )i = i - 1
      dx = 10.**i
      i = rng / dx
      if( i .gt. 10 )dx = 2. * dx
      if( i .gt. 20 )dx = 2.5 * dx
c decide format for output
      i = 2. * dx
      xx = 2. * dx - real( i )
      integ = ( xx .lt. .01 ) .and. ( rng .gt. 1. )
      if( .not. integ )ndec = .999 - log10( xx )
c ensure horizontal text for tickmark labels
      call jsrota( 0. )
c reduce size of scale notation
      mheight = height
      call jschsz( .7 * mheight )
c disable clipping
      call pgsclp( 0 )
      tickx = height * ( sx2 - sx1 ) / ( 3. * ( rx2 - rx1 ) )
      ticky = height * ( sy2 - sy1 ) / ( 3. * ( ry2 - ry1 ) )
c decide which axis
      if( a .eq. 'x' )then
c x-axis
        call jspenu( sx1, sy1 )
        call jspend( sx2, sy1 )
        call jspenu( sx1, sy2 )
        call jspend( sx2, sy2 )
c add tickmarks
        i = sx1 / dx - 1.
        x = dx * real( i )
        even = mod( i, 2 ) .eq. 0
c skip to first useful tick mark
        do while ( x + dx .lt. sx1 - tiny )
          x = x + dx
          even = .not. even
        end do
        do while ( x + dx .lt. sx2 + tiny )
          x = x + dx
          even = .not. even
          call jspenu( x, sy2 )
          ya = sy2 - ticky
          if( even )ya = sy2 - 2. * ticky
          call jspend( x, ya )
          call jspenu( x, sy1 )
          ya = sy1 + ticky
          if( even )ya = sy1 + 2. * ticky
          call jspend( x, ya )
          if( ( iflag .ne. 0 ) .and. even )then
c tickmark labels
            nw = 0
            if( x .ne. 0. )nw = log10( abs( x ) ) + .001
            nw = max( nw, 0 ) + 1
            if( x .lt. 0. )nw = nw + 1
            if( integ )then
              i = x
              call jsbldi( i, nw )
            else
              nw = nw + ndec + 1
              call jsbldf( x, nw, ndec )
            end if
c offset for tickmark labels
            xa = x - 1.8 * tickx * real( nw )
            ya = sy1 - 6. * ticky
            call jswrit( xa, ya )
          end if
        end do
c restore old character height
        call jschsz( mheight )
c add label if any
        nchr = len( label )
        if( nchr .gt. 0 )then
          xa = .5 * ( sx1 + sx2 - 1.8 * tickx * real( nchr ) )
          ya = sy1 - 10.5 * ticky
          call jsbldt( label )
          call jswrit( xa, ya )
        end if
      else
        if( .not. ( a .eq. 'y' ) )call jserr
c y - axis
        call jspenu( sx1, sy1 )
        call jspend( sx1, sy2 )
        call jspenu( sx2, sy1 )
        call jspend( sx2, sy2 )
c add tickmarks
        i = sy1 / dx - 1.
        y = dx * real( i )
        even = mod( i, 2 ) .eq. 0
c skip to first useful tick mark
        do while ( y + dx .lt. sy1 - tiny )
          y = y + dx
          even = .not. even
        end do
        nwmx = 0
        do while ( y + dx .lt. sy2 + tiny )
          y = y + dx
          even = .not. even
          call jspenu( sx2, y )
          xa = sx2 - tickx
          if( even )xa = sx2 - 2. * tickx
          call jspend( xa, y )
          call jspenu( sx1, y )
          xa = sx1 + tickx
          if( even )xa = sx1 + 2. * tickx
          call jspend( xa, y )
          if( ( iflag .ne. 0 ) .and. even )then
c tickmark labels
            nw = 0
            if( y .ne. 0. )nw = log10( abs( y ) ) + .001
            nw = max( nw, 0 ) + 1
            if( y .lt. 0. )nw = nw + 1
            if( integ )then
              i = y
              call jsbldi( i, nw )
            else
              nw = nw + ndec + 1
              call jsbldf( y, nw, ndec )
            end if
            nwmx = max( nw, nwmx )
c offset for tickmark labels
            xa = sx1 - 3. * tickx * real( nw + 1 )
            ya = y - 2.1 * ticky
            call jswrit( xa, ya )
          end if
        end do
c restore old character height
        call jschsz( mheight )
c add label if any
        nchr = len( label )
        if( nchr .gt. 0 )then
c rotate text for axis label
          call jsrota( 90. )
          xa = sx1 - 3. * tickx * real( nwmx + 2 )
          ya = .5 * ( sy1 + sy2 - 1.8 * ticky * real( nchr ) )
          call jsbldt( label )
          call jswrit( xa, ya )
c restore horizontal text
          call jsrota( 0. )
        end if
      end if
c restore clipping
      call pgsclp( 1 )
      call jsebuf
      return
      end
