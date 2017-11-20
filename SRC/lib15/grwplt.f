      subroutine grwplt( mf, istp, amp, phase, tim, istart, iend )
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c part of the analysis software
c plots the time evolution of first the amplitute and then the phase
c    of expansion coefficients measured from the simulation
c  the user can iteratively select the time range over which to fit for
c    a single overstability (a growing and rotating mode)
c
c calling arguments
      integer iend, istart, istp, mf
      real amp( istp ), phase( istp ), tim( istp )
c
c external
      real roundup
c
c local variables
      character*1 yn
      integer i, ip, n
      logical jscren
      real berr, bint, pif, pit, sl, slerr, x, x1, xmax, y
      real ymax, ymin
      include 'inc/pi.f'
c
c make phases monotonic
      pif = 0.
      pit = 0.
      if( mf .ne. 0 )pit = 2. * pi / real( mf )
      do i = 2, istp
        if( ( phase( i ) + pif ) .lt. ( phase( i - 1 ) - .5 * pit ) )
     +      pif = pif + pit
        phase( i ) = phase( i ) + pif
      end do
c set window
    2 xmax = ( int( tim( istp ) * .1 ) + 1 ) * 10
      do 3 ip = 1, 2
        x = .1 + .4 * real( ip - 1 )
        x1 = x + .3
        call jssize( x, x1, 0.1, 0.7 )
        if( ip .eq. 1 )then
          ymin = 0
          do i = 1, istp
            ymin = min( ymin, amp( i ) )
          end do
          ymin = roundup( ymin )
          ymax = 0.
        else
          ymin = -5.
          ymax = 0
          do i = 1, istp
            ymax = max( ymax, phase( i ) )
          end do
          ymax = roundup( ymax )
        end if
        call jscale( 0., xmax, ymin, ymax )
c draw axes and plot data
        call jsaxis( 'x', 'time', 4 )
        if( ip .eq. 2 )then
          call jsaxis( 'y', 'radians', 11 )
          call jsymbl( tim, phase, istp, 4 )
        else
          call jsaxis( 'y', 'ln(A)', 7 )
          call jsymbl( tim, amp, istp, 4 )
c select range for linear fit
          call jsebuf
          print *, 'Default time range for fit is', istart, iend
          call gtchar( ' Do you wish to change this? (y/n)', yn )
          if( yn .eq. 'y' )then
            call gtintg( 'Enter start time', istart )
            call gtintg( 'Enter end time', iend )
            istart = max( 1, istart )
            iend = min( iend, istp )
          end if
        end if
c determine best straight line
        n = iend - istart + 1
        if( ip .eq. 1 )then
          call linlsq( tim( istart ), amp( istart ), n, sl, bint,
     +                 slerr, berr )
        else
          call linlsq( tim( istart ), phase( istart ),n,sl, bint,
     +                 slerr, berr )
        end if
c draw best fit lines
        x = 0.
        if( bint .lt. ymin )x = ( ymin - bint ) / sl
        y = max( ymin, bint )
        call jsmove( x, y )
        x = tim( istp )
        y = bint + sl * x
        if( y .gt. ymax )x = ( ymax - bint ) / sl
        y = min( y, ymax )
        call jsline( x, y )
c output value of slope
        call jsbldt( 'slope =' )
        call jsbldf( sl, 6, 4 )
        call jsbldt( '+' )
        call jsbldf( slerr, 6, 4 )
        x = .2 * xmax
        y = .1 * real( 3 - ip )
        y = ( 1. - y ) * ymin + y * ymax
        call jswrit( x, y )
c mark range of values used
        x = tim( istart )
        y = amp( istart )
        if( ip .eq. 2 )y = phase( istart )
        call jsymbl( x, y, 1, 5 )
        x = tim( iend )
        y = amp( iend )
        if( ip .eq. 2 )y = phase( iend )
        call jsymbl( x, y, 1, 5 )
    3 continue
      if( .not. jscren( 0 ) )return
c enquire whether fit is acceptable
      y = 1.1 * ymin - .1 * ymax
      call jsmove( 0., y )
      call jsline( 0., y )
      call jsebuf
      call gtchar( 'Do you want to try another fit? ( y/n )', yn )
      if( yn .ne. 'y' )return
      call jspage
      go to 2
      end
