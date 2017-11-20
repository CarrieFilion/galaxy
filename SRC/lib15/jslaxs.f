      subroutine jslaxs( aa, label, iflag )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c draw axis and place logarithmically spaced tickmarks
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
      character a*1, c*3
      integer i, k, nchr, nw, nwmx
      logical integ, decade
      real mheight, rng, s1, s2, tickx, ticky, tiny, x, xa, xx, y, ya
c
      a = aa
      call lowercase( a )
c choose suitable spacing for tickmarks
      if( a .eq. 'x' )then
        s1 = sx1
        s2 = sx2
      else if( a .eq. 'y' )then
        s1 = sy1
        s2 = sy2
      else
        print *, 'In JSLAXS neither x- nor y-axis!'
        call jserr
        stop
      end if
      rng = s2 - s1
      tiny = 1.e-3 * rng
c start from integer value less than lowest point of range
      k = s1
      if( s1 .lt. 0. )k = k - 1
c decide format for output
      xx = 10.**s2
      integ = xx .ge. 1.
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
        xx = k
        do while ( xx .lt. sx2 + tiny )
          do i = 1, 10
            integ = i .eq. 1
            x = xx + log10( real( i ) )
            if( ( x .gt. sx1 - tiny ) .and. ( x .lt. sx2 + tiny ) )then
              decade = integ
              call jspenu( x, sy2 )
              ya = sy2 - ticky
              if( decade )ya = sy2 - 2. * ticky
              call jspend( x, ya )
              call jspenu( x, sy1 )
              ya = sy1 + ticky
              if( decade )ya = sy1 + 2. * ticky
              call jspend( x, ya )
              if( decade .and. ( iflag .ne. 0 ) )then
c tickmark labels
                k = nint( x )
                if( x .ge. 0. )then
                  if( k .lt. 4 )then
                    nw = k + 1
                    k = 10**k
                    call jsbldi( k, nw )
                  else
                    nw = 2
                    if( k .lt. 10 )then
                      write( c( 1:1 ), '( i1 )' )k
                      call jsbldt( '10^{'//c( 1:1 )//'}' )
                    else if( k .lt. 100 )then
                      write( c( 1:2 ), '( i2 )' )k
                      call jsbldt( '10^{'//c( 1:2 )//'}' )
                    else
                      write( c( 1:3 ), '( i3 )' )k
                      call jsbldt( '10^{'//c//'}' )
                    end if
                  end if
                else
                  if( k .ge. -2 )then
                    nw = 2 - k
                    xa = 10.**x
                    call jsbldf( xa, nw, -k )
                    nw = nw - 1
                  else
                    nw = 3
                    if( k .gt. -10 )then
                      write( c( 1:1 ), '( i1 )' )-k
                      call jsbldt( '10^{-'//c( 1:1 )//'}' )
                    else if( k .gt. -100 )then
                      write( c( 1:2 ), '( i2 )' )-k
                      call jsbldt( '10^{-'//c( 1:2 )//'}' )
                    else
                      write( c( 1:3 ), '( i3 )' )-k
                      call jsbldt( '10^{-'//c//'}' )
                    end if
                  end if
                end if
c offset for tickmark labels
                xa = x - 1.8 * tickx * real( nw )
                ya = sy1 - 6. * ticky
                call jswrit( xa, ya )
              end if
            end if
          end do
          xx = xx + 1
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
      else if( a .eq. 'y' )then
c y - axis
        call jspenu( sx1, sy1 )
        call jspend( sx1, sy2 )
        call jspenu( sx2, sy1 )
        call jspend( sx2, sy2 )
c add tickmarks
        nwmx = 0
        xx = k
        do while ( xx .lt. sy2 + tiny )
          do i = 1, 10
            integ = i .eq. 1
            y = xx + log10( real( i ) )
            if( ( y .gt. sy1 - tiny ) .and. ( y .lt. sy2 + tiny ) )then
              decade = integ
              call jspenu( sx2, y )
              xa = sx2 - tickx
              if( decade )xa = sx2 - 2. * tickx
              call jspend( xa, y )
              call jspenu( sx1, y )
              xa = sx1 + tickx
              if( decade )xa = sx1 + 2. * tickx
              call jspend( xa, y )
              if( decade .and. ( iflag .ne. 0 ) )then
c tickmark labels
                k = nint( y )
                if( y .ge. 0. )then
                  if( k .lt. 4 )then
                    nw = k + 1
                    k = 10**k
                    call jsbldi( k, nw )
                  else
                    nw = 2
                    if( k .lt. 10 )then
                      write( c( 1:1 ), '( i1 )' )k
                      call jsbldt( '10^{'//c( 1:1 )//'}' )
                    else if( k .lt. 100 )then
                      write( c( 1:2 ), '( i2 )' )k
                      call jsbldt( '10^{'//c( 1:2 )//'}' )
                    else
                      write( c( 1:1 ), '( i3 )' )k
                      call jsbldt( '10^{'//c//'}' )
                    end if
                  end if
                else
                  if( k .ge. -2 )then
                    nw = 2 - k
                    ya = 10.**y
                    call jsbldf( ya, nw, -k )
                    nw = nw - 1
                  else
                    nw = 3
                    if( k .gt. -10 )then
                      write( c( 1:1 ), '( i1 )' )-k
                      call jsbldt( '10^{-'//c( 1:1 )//'}' )
                    else if( k .gt. -100 )then
                      write( c( 1:2 ), '( i2 )' )-k
                      call jsbldt( '10^{-'//c( 1:2 )//'}' )
                    else
                      write( c( 1:1 ), '( i3 )' )-k
                      call jsbldt( '10^{-'//c//'}' )
                    end if
                  end if
                end if
                nwmx = max( nwmx, nw )
c offset for tickmark labels
                xa = sx1 - 3. * tickx * real( nw + 1 )
                ya = y - 2.1 * ticky
                call jswrit( xa, ya )
              end if
            end if
          end do
          xx = xx + 1
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
        end if
c restore horizontal text
        call jsrota( 0. )
      end if
c restore clipping
      call pgsclp( 1 )
      call jsebuf
      return
      end
