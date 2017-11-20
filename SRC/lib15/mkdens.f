      subroutine mkdens
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to plot density contrast as a color image at selected times
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/jscmmn.f'
c
      include 'inc/model.f'
c
      integer mrng
      parameter ( mrng = 501 )
      real rad( mrng )
      common / trcdent / rad
c
c externals
      logical gtlogl, jscren
      real resyn, uofr
c
c local arrays
      integer ig( 10 )
      real, allocatable :: work( :, : )
c
c local variables
      integer i, ifail, ir, ix, iy, j, jstep, kstep, mx, my, nskip
      parameter ( mx = 201, my = mx )
      real blank, du, r, scale, t1, u, x, y
      include 'inc/pi.f'
c
c works only with danl data - so far
      danl = .true.
c select component
      icmp = 0
      do i = 1, ncmp
        if( disc( i ) .and. ( nsp( i ) .gt. 0 ) .and.
     +      ( .not. testp( i ) ) )icmp = i
      end do
      if( icmp .le. 0 )call crash( 'MKDENS', 'No population selected' )
c select grid appropriate for population
      if( hybrid )then
        i = igrd( icmp )
        call switch( i )
      end if
c select times for plotting
      call lstrid( 'DANL', nskip, kstep )
      if( nskip .gt. 0 )then
        call gtreal( 'Enter time of first frame required', r )
        jstep = nint( r / ts )
c set radial range
        rmin = 0
        radm = rgrid( jgrid ) / lscale
        scale = 2. * radm / real( mx - 1 )
c range of Fourier components
        print *, 'Active sectoral harmonics are:'
        j = 0
        do i = 1, ng
          if( .not. lg( i ) )then
            j = j + 1
            if( j .gt. 10 )then
              print '( 10i8 )', ig
              j = 1
            end if
            ig( j ) = i - 1
          end if
        end do
        print '( 10i8 )', ( ig( i ), i = 1, j )
        call gtintgs( 'Enter mml, mmu', ig, 2 )
c force density contrast only
        mml = max( ig( 1 ), 1 )
        mmu = min( ig( 2 ), ng - 1 )
c decide whether relative or absolute overdensities
        bulge = gtlogl('Plot relative overdensity (no gives absolute)' )
c allocate space
        allocate ( work( mx, my ) )
        blank = -5
        if( jscren( 0 ) )blank = 5
c restart file
        call rewnd
        istep = -1
        ifail = 0
        do while ( ( istep .lt. jstep ) .and. ( ifail .eq. 0 ) )
          call nextrec( 'DANL', ifail )
        end do
c work through file
        do while ( ifail .eq. 0 )
          if( mod( istep, nskip ) .eq. 0 )then
            if( .not. jscren( 0 )
     +                       )print *, 'Plotting density at time ', time
c work over array
            j = 0
            do iy = 1, my
              y = scale * real( iy - 1 ) - radm
              do ix = 1, mx
                x = scale * real( ix - 1 ) - radm
                r = sqrt( x * x + y * y )
                if( r .lt. radm )then
                  t1 = atan2( y, x )
                  if( t1 .lt. 0. )t1 = t1 + 2 * pi
                  t1 = t1 / alpha + 1
                  u = uofr( r * lscale )
                  ir = u
                  du = u - real( ir )
                  u = ir
                  work( ix, iy ) = ( 1. - du ) * resyn( t1, u ) +
     +                                    du   * resyn( t1, u + 1. )
                else
                  work( ix, iy ) = blank
                end if
              end do
            end do
c set scale
            call jspage
            call jssize( .15, 1., .15, 1. )
            call jsescl( -radm, radm, -radm, radm )
c make plot
            call colplt( work, mx, my )
c label plot
            if( ipic .eq. 1 )then
c              call jsbldt( 'Run' )
              call jsbldi( irun, 4 )
              call jswrit( -.95 * radm, .85 * radm )
              call jsbldt( 'r_{max}' )
              call jsbldf( radm, 4, 1 )
              call jswrit( .5 * radm, -.92 * radm )
            end if
c            call jsbldt( 't =' )
            i = nint( time )
            if( abs( time - real( i ) ) .gt. .1 )then
              call jsbldf( time, 5, 1 )
            else
              call jsbldi( i, 4 )
            end if
            call jswrit( .5 * radm, .85 * radm )
          end if
          call nextrec( 'DANL', ifail )
          if( istep .gt. kstep )ifail = 1
        end do
      end if
      return
      end
