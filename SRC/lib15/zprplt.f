      subroutine zprplt
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
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
c local variables
      integer i, ifail, iskip, istart, iend, j, k
      real rbzf, rho, t, xmax, ymin, ymax, z, zbzf
c
      call rewnd
      icmp = 1
      call nextrec( 'ZPRF', ifail )
      if( ifail .eq. 1 )then
        print *, 'No ZPRF records found'
        return
      end if
c set scale variables
      i = max( 1, mr - 2 )
      rbzf = real( i ) / rtrunc( 1 )
      zbzf = real( ma / 2 ) * lscale / zm( jgrid )
      xmax = real( ma / 2 ) / zbzf
c select options
      call gtintg( 'Enter first step numnber required', istart )
      call gtintg( 'Enter last step, or 0 for all', iend )
      if( iend .eq. 0 )iend = 10**6
      call gtintg( 'Enter time step interval, or 0 for all', iskip )
      iskip = max( iskip, 1 )
c skip to first requested step
      do while ( istep .lt. istart )
        call nextrec( 'ZPRF', ifail )
        if( ifail .ne. 0 )return
      end do
c work through file
      do while ( ifail .eq. 0 )
        call jspage
c write page heading
        call jssize( 0., 1., 0., 1. )
        call jscale( 0., 1., 0., 1. )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        call jsbldt( '  time' )
        t = real( istep ) * ts
        k = t + .5
        if( abs( t - real( k ) ) .gt. .01 )then
          call jsbldf( t, 5, 1 )
        else
          call jsbldi( k, 3 )
        end if
        if( npic .gt. 1 )then
          call jswrit( .83, .97 )
        else
          call jswrit( .3, .97 )
        end if
c set scale
        ymax = 0
        do i = 1, mr
          k = ( i - 1 ) * ma
          do j = 1, ma
            k = k + 1
            ymax = max( ymax, wres( k ) )
          end do
        end do
        if( ymax .gt. 0. )then
          j = log( ymax )
          ymax = j + 1
          ymin = ymax - 12
c set frame and labels
          call jssize( .2, .95, .2, .9 )
          call jscale( -xmax, xmax, ymin, ymax )
          call jsaxis( 'x', 'z', 1 )
          call jsaxis( 'y', 'ln \\rho', 1 )
c          r = real( i - 1 ) / rbzf
c          call jsbldt( 'r=' )
c          call jsbldf( r, 5, 2 )
c          call jswrit( -0.9 * xmax, .9 * ymax + .1 * ymin )
        end if
c plot profiles
        do i = 1, mr
          k = ( i - 1 ) * ma
          do j = 1, ma
            k = k + 1
            z = real( j - 1 - ma / 2 ) / zbzf
            rho = wres( k )
            if( rho .gt. 0. )then
               rho = max( log( rho ), ymin )
            else
               rho = ymin
            end if
            if( j .eq. 1 )call jsmove( z, rho )
            call jsline( z, rho )
          end do
        end do
c get next step
        call nextrec( 'ZPRF', ifail )
        do while ( mod( istep, iskip ) .ne. 0 )
          call nextrec( 'ZPRF', ifail )
          if( ifail .ne. 0 )return
        end do
        if( istep .gt. iend )ifail = 1
      end do
      return
      end
