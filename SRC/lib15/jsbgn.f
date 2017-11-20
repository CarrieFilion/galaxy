      subroutine jsbgn
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to begin plotting
c   the plotting device to be selected, fraction of the available
c   surface that will be used, whether that should be the largest
c   inscribed square, and the number of sub-regions into which the
c   surface should be divided, are all determined from standard input
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local arrays
      integer ndev
      parameter ( ndev = 9 )
      character*25 device( ndev ), wkstn( ndev )
      logical ldev( ndev )
      equivalence ( ldev( 1 ), printx )
c
c local variables
      character yn*1
      integer i, idev
      logical square
      real dim, pheight, pwidth, xd, yd
c
c data statements
c
      data ldev / ndev * .false. /
c
      data device / 'null',
     +              'xwindow',
     +              'xwindow',
     +              'postscript (portrait)',
     +              'postscript (landscape)',
     +              'color ps (portrait)',
     +              'color ps (landscape)',
     +              'gif (portrait)',
     +              'gif (landscape)' /
c table of valid workstation names
      data wkstn /  '/null',
     +              '/xwindow',
     +              '/xwindow',
     +              '/vps',
     +              '/ps',
     +              '/vcps',
     +              '/cps',
     +              '/vgif',
     +              '/gif' /
c
c initialise plot
      idev = -1
      do while ( ( idev .lt. 1 ) .or. ( idev .gt. ndev ) )
        print *, 'Choose device number'
        do i = 1, ndev
          print '( i5, 3x, a25 )', i - 1, device( i )
        end do
        read( *, *, iostat = i )idev
        if( i .eq. 0 )then
          null = idev .eq. 0
          idev = idev + 1
        else
          idev = -1
        end if
      end do
      screen = ( idev .eq. 2 ) .or. ( idev .eq. 3 )
      ldev( idev ) = .true.
c
c initialise device
c pgplot - device dimensions in inches, my units in cm!
      call pgbegin( 0, wkstn( idev ), 1, 1 )
      call pgqvp( 2, fx1, fx2, fy1, fy2 )
!      print '( ''Device dimensions in cm'', 4f8.3 )',
!     +                    .1 * fx1, .1 * fx2, .1 * fy1, .1 * fy2
      pwidth = .1 * ( fx1 + fx2 )
      pheight = .1 * ( fy1 + fy2 )
      call pgqvp( 3, fx1, fx2, fy1, fy2 )
!      print '( ''Device dimensions in pixels'', 4f8.2 )',
!     +                    fx1, fx2, fy1, fy2
c
c pixpcm is the number of device units ( pixels / inches / cm... ) per cm
c
      pixpcm = ( fx1 + fx2 ) / pwidth
c use these as default
      xcmmax = pwidth
      ycmmax = pheight
      if( null )then
c set defaults for no graphic output
        square = .false.
        nx = 1
        ny = 1
      else
        if( .not. screen )then
          print 200, pwidth, pheight
  200     format( ' Do you wish to change the paper size from',
     +           f7.3, ' by', f7.3, ' cm?' )
          read '( a1 )', yn
          call lowercase( yn )
          if( yn .eq. 'y' )then
            print *, 'Enter new dimensions in cm'
            read *, xcmmax, ycmmax
c limit to possible paper size
            xcmmax = min( xcmmax, pwidth )
            ycmmax = min( ycmmax, pheight )
          end if
        end if
c square frames?
        print *, 'Do you want square frames?'
        read '( a1 )', yn
        call lowercase( yn )
        square = yn .eq. 'y'
c number of pictures per page
    2   print *, 'Enter no of pictures across and down page'
        read ( *, *, err = 2 ), nx, ny
        nx = max( nx, 1 )
        ny = max( ny, 1 )
      end if
      npic = nx * ny
c compute sizes of sub-frames
      xdim = xcmmax / real( nx )
      ydim = ycmmax / real( ny )
c impose square sub-frames if requested
      if( square )then
        dim = min( xdim, ydim )
        xdim = dim
        ydim = dim
      end if
c size of plotting area and offsets
      xd = xdim * real( nx )
      yd = ydim * real( ny )
      xoffjs = .5 * ( pwidth - xd )
      yoffjs = .5 * ( pheight - yd )
      ipic = 0
      if( screen )ipic = npic
c
      call pgask( .false. )
c set boundaries of plotting surface in cm
      vx1 = xoffjs
      vx2 = xoffjs + xd
      vy1 = yoffjs
      vy2 = yoffjs + yd
c pgplot works in inches
      call pgvsiz( vx1 / 2.54, vx2 / 2.54, vy1 / 2.54, vy2 / 2.54 )
      call pgswin( 0., xd, 0., yd )
c set default character size to 0.3 cm
      height = .3
      angle = 0.
c initialise character string
      nchar = 1
c set first window
      fx1 = vx1
      fx2 = fx1 + xdim
      fy1 = vy1 + ydim * real( ny - 1 )
      fy2 = fy1 + ydim
c set default window size and scaling
      call jssize( .15, .95, .15, .95 )
      call jscale( 0., 1., 0., 1. )
c initial line weight
      call jsthik( 3 )
c initial character size
      call jschsz( height )
c initial line style
      dash = .false.
      partd = 0
c initial pen position
      xlast = 0
      ylast = 0
      xprev = 0
      yprev = 0
      outbox = .false.
      return
      end
