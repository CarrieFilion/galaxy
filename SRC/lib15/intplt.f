      subroutine intplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Displays, with suitable annotation, the time evolution of various integral
c   quantities saved from a simulation in the runXXX.res file
c
c Called from ANALYS
c Graphics routines are from JSPLOT
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
      include 'inc/comprns.f'
c
      include 'inc/grids.f'
c
      include 'inc/jscmmn.f'
c
      include 'inc/model.f'
c
c externals
      real roundup
c
c local arrays
      integer msize
      parameter ( msize = 20001 )
      integer nspop( mcmp, msize ), nlost( msize )
      real cx( msize ), cy( msize ), cz( msize ), popm( mcmp, msize )
      real px( msize ), py( msize ), pz( msize )
      real fx( msize ), fy( msize ), fz( msize )
      real hx( msize ), hy( msize ), hz( msize )
      real hxo( msize ), hyo( msize ), hzo( msize )
      real pe( msize ), ke( msize ), te( msize ), top( msize )
      real vc( msize ), t( msize ), x( msize ), y( msize ), z( msize )
      real Lz( 3, msize ), xl2( 2 ), yl2( 2 )
c
c local variables
      character*3 opt
      integer i, ifail, ipp, j, jdisc, jpop, k, krun, n, nlive
      logical linit, offanal
      real actmas, den, dv, fac, hdmax, hnorm, moi, rk, tmax, tsc
      real x1, x2, x3, x4, xl, xl1, xmax, y1, y2, y3, y4, yl, yl1
      real ymax, ymin, ymn( 2 ), ymx( 2 ), zmax
      equivalence ( k, rk )
      parameter ( dv = .1 )
c
      linit = .false.
c plotting loop
      opt = '   '
      do while ( opt .ne. 'end' )
        krun = 1
c restart file
    5   call rewnd
c read in data
        if( .not. linit )then
          n = 0
          call nextrec( 'INTG', ifail )
c
          do while ( ifail .eq. 0 )
            n = n + 1
            if( n .gt. msize )call space( msize, n, 'arrays', 'INTPLT' )
            t( n ) = time
            jpop = ma
c total energy, virial and tOP
            te( n ) = wres( 1 )
            vc( n ) = wres( 2 )
            top( n ) = wres( 3 )
c pick out integers
            rk = wres( 4 )
            nlost( n ) = k
            j = 5
            nlive = 0
            do i = 1, jpop
              if( uqmass )then
                popm( i, n ) = wres( j )
                if( popm( i, n ) .gt. 0. )then
                  nlive = nlive + 1
                  if( disc( i ) )jdisc = i
                end if
              else
                rk = wres( j )
                nspop( i, n ) = k
                if( nspop( i, n ) .gt. 0 )then
                  nlive = nlive + 1
                  if( disc( i ) )jdisc = i
                end if
              end if
              j = j + 15
            end do
c only one pop expected for a 2D model
            if( twod )then
              cx( n ) = wres( 6 )
              cy( n ) = wres( 7 )
              px( n ) = wres( 9 )
              py( n ) = wres( 10 )
              fx( n ) = wres( 12 )
              fy( n ) = wres( 13 )
              hz( n ) = wres( 15 )
              pe( n ) = wres( 18 )
              ke( n ) = wres( 19 )
              hzo( n ) = wres( 22 )
            end if
            cx( n ) = 0
            cy( n ) = 0
            cz( n ) = 0
            px( n ) = 0
            py( n ) = 0
            pz( n ) = 0
            fx( n ) = 0
            fy( n ) = 0
            fz( n ) = 0
            hx( n ) = 0
            hy( n ) = 0
            hz( n ) = 0
            pe( n ) = 0
            ke( n ) = 0
c sum over active components
            j = 5
            den = 0
            do ipp = 1, jpop
              if( uqmass )then
                fac = popm( ipp, n )
                den = den + fac
              else
                fac = nspop( ipp, n )
                den = den + fac
              end if
              cx( n ) = cx( n ) + wres( j + 1 ) * fac
              cy( n ) = cy( n ) + wres( j + 2 ) * fac
              if( threed )cz( n ) = cz( n ) + wres( j + 3 ) * fac
              px( n ) = px( n ) + wres( j + 4 )
              py( n ) = py( n ) + wres( j + 5 )
              if( threed )pz( n ) = pz( n ) + wres( j + 6 )
              fx( n ) = fx( n ) + wres( j + 7 )
              fy( n ) = fy( n ) + wres( j + 8 )
              if( threed )fz( n ) = fz( n ) + wres( j + 9 )
              hz( n ) = hz( n ) + wres( j + 10 )
              if( threed )then
                hy( n ) = hy( n ) + wres( j + 11 )
                hx( n ) = hx( n ) + wres( j + 12 )
              end if
              if( ipp .le. 3 )then
                Lz( ipp, n ) = wres( j + 10 )
              else
                call crash( 'INTPLT', 'Too many pops' )
              end if
              pe( n ) = pe( n ) + wres( j + 13 )
              ke( n ) = ke( n ) + wres( j + 14 )
              j = j + 15
            end do
            cx( n ) = cx( n ) / den
            cy( n ) = cy( n ) / den
            if( threed )cz( n ) = cz( n ) / den
c angular momentum in particles off the grid
            hzo( n ) = 0
            hyo( n ) = 0
            hxo( n ) = 0
            j = j - 1
            do ipp = 1, ncmp
              hzo( n ) = hzo( n ) + wres( j + 1 )
              if( threed )then
                hyo( n ) = hyo( n ) + wres( j + 2 )
                hxo( n ) = hxo( n ) + wres( j + 3 )
              end if
              j = j + 3
            end do
            if( pertbn )then
c read satellite data also to get the full story - PE term should be ignored
c   because it was included already by a call to satpot for each particle
              call nextrec( 'SATE', ifail )
              if( ifail .ne. 0 )call crash( 'INTPLT',
     +                                  'Satellite data not found' )
c energy terms
              if( exsphr .or. exgnrc )then
                if( mr .eq. 11 )mptbr = wres( 11 )
                ke( n ) = ke( n ) +
     +       .5 * mptbr * ( wres( 4 )**2 + wres( 5 )**2 + wres( 6 )**2 )
                te( n ) = pe( n ) + ke( n )
c momentum etc
                actmas = pmass * den / ( lscale**3 * ts**2 )
c            fac = mptbr / actmas
c            cx( n ) = cx( n ) + fac * wres( 1 )
c            cy( n ) = cy( n ) + fac * wres( 2 )
c            cz( n ) = cz( n ) + fac * wres( 3 )
                fac = actmas + mptbr
                cx( n ) = ( cx( n ) * actmas + mptbr * wres( 1 ) ) / fac
                cy( n ) = ( cy( n ) * actmas + mptbr * wres( 2 ) ) / fac
                cz( n ) = ( cz( n ) * actmas + mptbr * wres( 3 ) ) / fac
                fac = mptbr
                px( n ) = px( n ) + fac * wres( 4 )
                py( n ) = py( n ) + fac * wres( 5 )
                pz( n ) = pz( n ) + fac * wres( 6 )
                hx( n ) = hx( n ) +
     +         mptbr * ( wres( 6 ) * wres( 2 ) - wres( 5 ) * wres( 3 ) )
                hy( n ) = hy( n ) +
     +         mptbr * ( wres( 4 ) * wres( 3 ) - wres( 6 ) * wres( 1 ) )
                hz( n ) = hz( n ) +
     +         mptbr * ( wres( 5 ) * wres( 1 ) - wres( 4 ) * wres( 2 ) )
              else if( extbar )then
                if( ( n .eq. 1 ) .and. ( eptbr .eq. 0.d0 ) )then
                  eptbr = .5
                  print *, 'Bar axis ratio 2:1 assumed'
                end if
                moi = mptbr * ( 1. + eptbr**2 ) * abar**2 / 5.
                ke( n ) = ke( n ) + .5 * moi * wres( 2 )**2
                te( n ) = pe( n ) + ke( n )
                Lz( 2, n ) = moi * wres( 2 )
                hz( n ) = Lz( 1, n ) + Lz( 2, n )
              else if( .not. exspir )then
                call crash( 'INTPLT', 'Unrecognized perturber' )
              end if
            end if
            if( n .gt. 1 )then
              hz( n ) = ( hz( n ) + hzo( n ) - hz( 1 ) ) / hnorm
              hy( n ) = ( hy( n ) + hyo( n ) - hy( 1 ) ) / hnorm
              hx( n ) = ( hx( n ) + hxo( n ) - hx( 1 ) ) / hnorm
              hdmax = max( hdmax, abs( hz( n ) ), abs( hy( n ) ),
     +             abs( hx( n ) ) )
            else
              if( twod )then
                hnorm = .01 * abs( hz( 1 ) )
              else
                hnorm =
     +                .01 * sqrt( hx( 1 )**2 + hy( 1 )**2 + hz( 1 )**2 )
              end if
              if( hnorm .lt. .0001 )hnorm = 1.
              hdmax = 0
            end if
c get next record
            call nextrec( 'INTG', ifail )
c stop reading subsequent files at the same step as the first
            if( ifail .eq. 0 )then
              if( krun .eq. 1 )then
                tmax = time
              else
                if( time .gt. tmax )ifail = 1
              end if
            end if
          end do
c report
          if( krun .eq. 1 )then
            if( n .eq. 0 )then
              print *, 'No INTG data found'
              return
            else if( n .eq. 1 )then
              print *, 'Only one INTG record, plots not yet possible'
              print *, 'Initial virial ratio', ke( 1 ) / abs( vc( 1 ) )
              return
            else
              print *, n, ' moments read in to time', time
            end if
          end if
          linit = nruns .eq. 1
        end if
c plotting options follow
        hz( 1 ) = 0
        hy( 1 ) = 0
        hx( 1 ) = 0
c CoM or net accelerations
        if( ( opt .eq. 'com' ) .or. ( opt .eq. 'acc' ) )then
          if( opt .eq. 'com' )then
            do i = 1, n
              x( i ) = cx( i )
              y( i ) = cy( i )
              z( i ) = cz( i )
            end do
          else
            do i = 1, n
              x( i ) = fx( i )
              y( i ) = fy( i )
              z( i ) = fz( i )
            end do
          end if
          if( krun .eq. 1 )then
            xmax = .01
            do i = 1, n
              xmax = max( xmax, abs( x( i ) ), abs( y( i ) ) )
              if( threed )xmax = max( xmax, abs( z( i ) ) )
            end do
            xmax = roundup( xmax )
            ymax = xmax
            if( threed )zmax = xmax
c draw and mark axes
            if( twod )then
              x1 = .1
              x2 = x1 + .8
              y1 = .1
              y2 = y1 + .8
              xl = -.7 * xmax
              yl = 1.05 * ymax
            end if
            if( threed )then
              den = min( fx2 - fx1, fy2 - fy1 )
              xl = .4 * den / ( fx2 - fx1 )
              x2 = .55
              x1 = x2 - xl
              x3 = x2
              x4 = x3 + xl
              yl = .4 * den / ( fy2 -fy1 )
              y2 = .55
              y1 = y2 - yl
              y3 = y2
              y4 = y3 + yl
              xl = .3 * xmax
              yl = 3.05 * ymax
            end if
            call jspage
          end if
          call jssize( x1, x2, y1, y2 )
          call jsescl( -xmax, xmax, -ymax, ymax )
          if( nruns .gt. 1 )then
            call pgsclp( 0 )
            xl2( 1 ) = 1.3 * xmax
            yl1 = .15 * real( krun )
            yl2( 1 ) = ( 1.1 + yl1 ) * ymax + ( 1. - yl1 ) * ymin
            xl2( 2 ) = xl2( 1 ) + .2 * xmax
            yl2( 2 ) = yl2( 1 )
            call jscoin( xl2, yl2, 2, krun )
            call jsbldi( irun, 5 )
            xl1 = xl2( 2 ) + .05 * xmax
            yl1 = yl2( 2 ) - .05 * ymax
            call jswrit( xl1, yl1 )
            call pgsclp( 1 )
          end if
          if( krun .eq. 1 )then
            call jsaxis( 'x', 'x', 1 )
            call jsaxis( 'y', 'y', 1 )
c mark header
            if( opt .eq. 'com' )then
              call jsbldt( 'Center of Mass' )
            else
              call jsbldt( 'Net acceleration' )
            end if
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 4 )
            end if
            call jswrit( xl, yl )
          end if
c plot data
          call jscoin( x, y, n, krun )
          if( threed )then
            call jssize( x3, x4, y1, y2 )
            call jsescl( -xmax, xmax, -ymax, ymax )
            if( krun .eq. 1 )then
              call jsaxis( 'x', 'z', 1 )
              call jsaxis( 'y', ' ', 0 )
            end if
            call jscoin( z, y, n, krun )
            call jssize( x1, x2, y3, y4 )
            call jsescl( -xmax, xmax, -ymax, ymax )
            if( krun .eq. 1 )then
              call jsaxis( 'x', ' ', 0 )
              call jsaxis( 'y', 'z', 1 )
            end if
            call jscoin( x, z, n, krun )
          end if
        end if
c
        if( opt .eq. 'vir' )then
c set scaling
          if( krun .eq. 1 )then
            xmax = t( n )
            ymax = .5 + dv
            ymin = .5 - dv
          end if
          do i = 1, n
c estimate virial ratio using definition of Clausius
            x( i ) = ke( i ) / abs( vc( i ) )
            if( krun .eq. 1 )then
              ymax = max( ymax, x( i ) )
              ymin = min( ymin, x( i ) )
            end if
          end do
          if( krun .eq. 1 )then
            y1 = abs( ymin - .5 + dv )
            y2 = abs( ymax - .5 - dv )
            if( y1 .gt. 1.e-3 )ymin = .5 + 1.1 * ( ymin - .5 )
            if( y2 .gt. 1.e-3 )ymax = .5 + 1.1 * ( ymax - .5 )
c draw and mark axes
            call jspage
            call jssize( .15, .95, .1, .9 )
            call jscale( 0., xmax, ymin, ymax )
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', 'Virial ratio', 1 )
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 4 )
              call jswrit( .4 * xmax, 1.05 * ymax - 0.05 * ymin )
            end if
          end if
          if( nruns .gt. 1 )then
            xl1 = .06 * xmax
            yl1 = .04 * real( krun )
            yl1 = yl1 * ymax + ( 1. - yl1 ) * ymin
            call jsbldi( irun, 5 )
            call jschsz( .25 )
            call jswrit( xl1, yl1 )
            call jschsz( .3 )
            xl2( 1 ) = xl1 - .04 * xmax
            yl2( 1 ) = yl1 + .01 * ( ymax - ymin )
            xl2( 2 ) = xl2( 1 ) + .03 * xmax
            yl2( 2 ) = yl2( 1 )
            call jscoin( xl2, yl2, 2, krun )
          end if
c plot virial ratio using definition of Clausius
          call jscoin( t, x, n, krun )
c plot virial ratio using PE
c          do i = 1, n
c            x( i ) = ke( i ) / abs( pe( i ) )
c          end do
c          call jsdash( 2, 2, 2, 2 )
c          call jsjoin( t, x, n )
c          call jsdash( 0, 0, 0, 0 )
        end if
c
        if( opt .eq. 'ene' )then
c total energy differences
          if( krun .eq. 1 )then
            ymax = 0
            ymin = 0
            do i = 1, n
              ymax = max( ymax, ke( i ), pe( i ) )
              ymin = min( ymin, ke( i ), pe( i ) )
            end do
c set scaling
            xmax = t( n )
            ymax = roundup( ymax )
            ymin = -roundup( -ymin )
            ymn( 1 ) = ymin
            ymx( 1 ) = ymax
            call jspage
          else
            ymin = ymn( 1 )
            ymax = ymx( 1 )
          end if
c draw and mark axes
          call jssize( .15, .95, .6, .95 )
          call jscale( 0., xmax, ymin, ymax )
          if( krun .eq. 1 )then
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', 'PE & KE', 1 )
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 4 )
              call jswrit( .4 * xmax, 1.05 * ymax - 0.05 * ymin )
            end if
          end if
c plot data
          if( nruns .gt. 1 )then
            xl1 = .06 * xmax
            yl1 = .3 + .06 * real( krun )
            yl1 = yl1 * ymax + ( 1. - yl1 ) * ymin
            call jsbldi( irun, 5 )
            call jschsz( .25 )
            call jswrit( xl1, yl1 )
            call jschsz( .3 )
            xl2( 1 ) = xl1 - .04 * xmax
            yl2( 1 ) = yl1 + .01 * ( ymax - ymin )
            xl2( 2 ) = xl2( 1 ) + .03 * xmax
            yl2( 2 ) = yl2( 1 )
            call jscoin( xl2, yl2, 2, krun )
          end if
          call jscoin( t, pe, n, krun )
          call jscoin( t, ke, n, krun )
c total energy differences
          if( krun .eq. 1 )ymax = .9
          do i = 1, n
            x( i ) = 100. * ( te( i ) - te( 1 ) ) / abs( te( 1 ) )
            if( krun .eq. 1 )ymax = max( ymax, abs( x( i ) ) )
          end do
c set scaling
          if( krun .eq. 1 )then
            xmax = t( n )
            ymax = roundup( ymax )
            ymx( 2 ) = ymax
          else
            ymax = ymx( 2 )
          end if
c draw and mark axes
          call jssize( .15, .95, .15, .5 )
          call jscale( 0., xmax, -ymax, ymax )
          if( krun .eq. 1 )then
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', '% energy change', 1 )
          end if
c plot data
          call jscoin( t, x, n, krun )
        end if
c
        if( ( opt .eq. 'lmo' ) .or. ( opt .eq. 'amo' ) )then
c set scaling
          if( krun .eq. 1 )then
            xmax = t( n )
            ymax = 0.1
            if( opt .eq. 'amo' )then
              ymax = roundup( hdmax )
              ymax = max( ymax, 1.e-4 )
            else
              ymax = 0
              do i = 1, n
                ymax = max( ymax, abs( px( i ) ),
     +                    abs( py( i ) ), abs( pz( i ) ) )
              end do
              ymax = roundup( ymax )
            end if
c draw and mark axes
            call jspage
          end if
          if( threed .or. ( opt .eq. 'amo' ) )then
            call jssize( .19, .99, .71, .99 )
            call jscale( 0., xmax, -ymax, ymax )
            if( krun .eq. 1 )then
              if( threed )then
                call jsaxis( 'x', ' ', 0 )
              else
                call jsaxis( 'x', 'time', 1 )
              end if
            end if
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 4 )
              call jswrit( .05 * xmax, .8 * ymax )
            end if
          end if
c
          if( opt .eq. 'lmo' )then
            if( krun .eq. 1 )then
              call jsaxis( 'y', 'z momentum', 1 )
            end if
            if( nruns .gt. 1 )then
              xl2( 1 ) = .05 * xmax
              yl1 = .05 + .07 * real( krun )
              yl2( 1 ) = yl1 * ymax - ( 1. - yl1 ) * ymax
              xl2( 2 ) = xl2( 1 ) + .05 * xmax
              yl2( 2 ) = yl2( 1 )
              call jscoin( xl2, yl2, 2, krun )
              call jsbldi( irun, 5 )
              xl1 = xl2( 2 ) + .01 * xmax
              yl1 = yl2( 2 ) - .05 * ymax
              call jswrit( xl1, yl1 )
            end if
            call jscoin( t, pz, n, krun )
          else
            if( krun .eq. 1 )then
              if( hnorm .eq. 1. )then
                call jsaxis( 'y', '\\Deltaz ang mom', 1 )
              else
                call jsaxis( 'y', '% z ang mom change', 1 )
              end if
            end if
            if( nruns .gt. 1 )then
              xl2( 1 ) = .05 * xmax
              yl1 = .05 + .07 * real( krun )
              yl2( 1 ) = yl1 * ymax - ( 1. - yl1 ) * ymax
              xl2( 2 ) = xl2 ( 1 ) + .05 * xmax
              yl2( 2 ) = yl2( 1 )
              call jscoin( xl2, yl2, 2, krun )
              call jsbldi( irun, 5 )
              xl1 = xl2( 2 ) + .01 * xmax
              yl1 = yl2( 2 ) - .05 * ymax
              call jswrit( xl1, yl1 )
            end if
            call jscoin( t, hz, n, krun )
          end if
c
          if( threed .or. ( opt .eq. 'lmo' ) )then
            call jssize( .19, .99, .41, .69 )
            call jscale( 0., xmax, -ymax, ymax )
            if( krun .eq. 1 )then
              call jsaxis( 'x', ' ', 0 )
              if( ( .not. threed ) .and. ( nruns .eq. 1 ) )then
                call jsbldt( 'Run no' )
                call jsbldi( irun, 4 )
                call jswrit( .05 * xmax, 1.1 * ymax )
              end if
            end if
            if( opt .eq. 'lmo' )then
              if( krun .eq. 1 )call jsaxis( 'y', 'y momentum', 1 )
              call jscoin( t, py, n, krun )
            else
              if( krun .eq. 1 )then
                if( hnorm .eq. 1. )then
                  call jsaxis( 'y', '\\Deltay ang mom', 1 )
                else
                  call jsaxis( 'y', '% y ang mom change', 1 )
                end if
              end if
              call jscoin( t, hy, n, krun )
            end if
c
            call jssize( .19, .99, .11, .39 )
            call jscale( 0., xmax, -ymax, ymax )
            if( krun .eq. 1 )call jsaxis( 'x', 'time', 1 )
            if( opt .eq. 'lmo' )then
              if( krun .eq. 1 )call jsaxis( 'y', 'x momentum', 1 )
              call jscoin( t, px, n, krun )
            else
              if( krun .eq. 1 )then
                if( hnorm .eq. 1. )then
                  call jsaxis( 'y', '\\Deltax ang mom', 1 )
                else
                  call jsaxis( 'y', '% x ang mom change', 1 )
                end if
              end if
              call jscoin( t, hx, n, krun )
            end if
          end if
        end if
c
        if( opt .eq. 'dht' )then
          if( ( extbar .and. jpop .eq. 1 ) .or. ( nlive .gt. 1 ) )then
            if( krun .eq. 1 )then
              xmax = t( n )
              y2 = 0
              y1 = 0
            end if
c find range of Lz and form total Lz
            do i = 1, n
              y( i ) = 0
              do j = 1, nlive
                y( i ) = y( i ) + Lz( j, i )
                if( krun .eq. 1 )then
                  y2 = max( y2, Lz( j, i ) )
                  y1 = min( y1, Lz( j, i ) )
                end if
              end do
              if( krun .eq. 1 )y2 = max( y2, y( i ) )
            end do
            if( krun .eq. 1 )then
              call jspage
              y2 = roundup( y2 )
              if( y1 .lt. 0. )y1 = -roundup( -y1 )
            end if
c set up plot
            call jssize( .15, .5, .15, .95 )
            call jscale( 0., xmax, y1, y2 )
            if( krun .eq. 1 )then
              call jsaxis( 'x', 't', 1 )
              if( nruns .eq. 1 )then
                call jsbldt( 'Run no' )
                call jsbldi( irun, 4 )
                call jswrit( .05 * xmax, 1.1 * y2 )
              end if
              call jsaxis( 'y', 'L_z', 1 )
            end if
            if( nruns .gt. 1 )then
              xl2( 1 ) = .1 * xmax
              yl1 = .2 + .03 * real( krun )
              yl2( 1 ) = yl1 * y2
              xl2( 2 ) = xl2( 1 ) + .1 * xmax
              yl2( 2 ) = yl2( 1 )
              call jscoin( xl2, yl2, 2, krun )
              call jsbldi( irun, 5 )
              xl1 = xl2( 2 ) + .02 * xmax
              yl1 = yl2( 2 ) - .01 * y2
              call jswrit( xl1, yl1 )
            end if
c plot total Lz
            call jscoin( t, y, n, krun )
c plot separate components
            if( ( nlive .gt. 1 ) .or. extbar )then
              call jsdash( 2, 2, 2, 2 )
              do j = 1, nlive
                do i = 1, n
                  y( i ) = Lz( j, i )
                end do
c                if( j .ne. jdisc .and. nlive .gt. 2 )call pgsci( j + 1 )
                call jscoin( t, y, n, krun )
                x1 = -100
                if( extbar )then
                  if( nlive .eq. 1 )then
c torque on halo from the bar
                    call linlsq( t, y, n, x1, y1, x2, y2 )
                  else
                    call crash( 'INTPLT',
     +                        'More than 1 live component + rigid bar' )
                  end if
                else
c compute torque on live bulge/ halo
                  if(
     +              j .ne. jdisc )call linlsq( t, y, n, x1, y1, x2, y2 )
                end if
c report torque(s)
                if( ( nruns .eq. 1 ) .and. ( x1 .gt. -100. ) )then
                  call jsbldt( 'Mean torque' )
                  y1 = .9
                  if( nlive .gt. 2 )then
                    call jsbldt( 'on comp' )
                    call jsbldi( j, 1 )
                    y1 = .95 - .05 * real( j )
                  end if
                  call jsblde( x1, 12, 4 )
                  call jswrit( .1 * xmax, y1 * y2 )
                end if
              end do
              call jsdash( 0, 0, 0, 0 )
            end if
c instantaneous torque
            if( krun .eq. 1 )then
              y4 = 0
              y3 = 0
              tsc = ( t( n ) - t( 1 ) ) / real( n - 1 )
            end if
c compute y range
            do j = 1, nlive
              if( j .ne. jdisc )then
                y( 1 ) = 0
                do i = 2, n
                  y( i ) = ( Lz( j, i ) - Lz( j, i - 1 ) ) / tsc
                end do
                k = 9
                if( n .gt. 180 )k = n / 20
                call bcsmth( y, n, k )
                if( krun .eq. 1 )then
                  do i = 2, n
                    y4 = max( y4, y( i ) )
                    y3 = min( y3, y( i ) )
                  end do
                end if
              end if
            end do
c set up plot
            if( krun .eq. 1 )then
              y4 = roundup( y4 )
              if( y3 .lt. 0. )y3 = -roundup( -y3 )
            end if
            if( y4 .gt. y3 )then
              call jssize( .6, .95, .15, .95 )
              call jscale( 0., xmax, y3, y4 )
              if( krun .eq. 1 )then
                call jsaxis( 'x', 't', 1 )
                if( nruns .eq. 1 )then
                  call jsbldt( 'Run no' )
                  call jsbldi( irun, 4 )
                  call jswrit( .05 * xmax, 1.02 * y4 )
                end if
                call jsaxis( 'y', '\\tau_z', 1 )
              end if
c recompute torques and plot
              do j = 1, nlive
                if( j .ne. jdisc )then
                  y( 1 ) = 0
                  do i = 2, n
                    y( i ) = ( Lz( j, i ) - Lz( j, i - 1 ) ) / tsc
                  end do
                  call bcsmth( y, n, k )
c                  if( nlive .gt. 2 )call pgsci( j + 1 )
                  call jscoin( t( k / 2 + 1 ), y, n - k + 1, krun )
                end if
              end do
            else
              print *, 'difference too small to plot'
            end if
          end if
        end if
c
        if( opt .eq. 'top' )then
c set scaling
          if( krun .eq. 1 )then
            xmax = t( n )
            ymax = .5
c draw and mark axes
            call jspage
            call jssize( .15, .95, .1, .9 )
            call jscale( 0., xmax, 0., ymax )
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', 't_{OP}', 1 )
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 4 )
              call jswrit( .4 * xmax, 1.05 * ymax )
            end if
          end if
c plot data
          if( nruns .gt. 1 )then
            xl1 = .06 * xmax
            yl1 = .97 - .05 * real( krun )
            yl1 = yl1 * ymax
            call jsbldi( irun, 5 )
            call jschsz( .25 )
            call jswrit( xl1, yl1 )
            call jschsz( .3 )
            xl2( 1 ) = xl1 - .04 * xmax
            yl2( 1 ) = yl1 + .01 * ymax
            xl2( 2 ) = xl2( 1 ) + .03 * xmax
            yl2( 2 ) = yl2( 1 )
            call jscoin( xl2, yl2, 2, krun )
          end if
          call jscoin( t, top, n, krun )
        end if
c escapers
        if( opt .eq. 'esc' )then
          if( krun .eq. 1 )ymax = 1
          do i = 1, n
            if( offanal )then
              k = 0
              if( uqmass )k = nlost( n )
            else
              k = nlost( n )
            end if
            do j = 1, jpop
              if( uqmass )then
c only initial number is available
                k = k + nsp( j )
              else
c use current number
                k = k + nspop( j, n )
              end if
            end do
            x( i ) = 100 * real( nlost( i ) ) / real( k )
            if( krun .eq. 1 )ymax = max( ymax, x( i ) )
          end do
c set scaling
          if( krun .eq. 1 )then
            xmax = t( n )
            ymax = roundup( ymax )
c draw and mark axes
            call jspage
            call jssize( .15, .95, .1, .9 )
            call jscale( 0., xmax, 0., ymax )
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', '% escapers', 1 )
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 4 )
              call jswrit( .4 * xmax, 1.05 * ymax )
            end if
          end if
c plot data
          if( nruns .gt. 1 )then
            xl1 = .06 * xmax
            yl1 = .97 - .05 * real( krun )
            yl1 = yl1 * ymax
            call jsbldi( irun, 5 )
            call jschsz( .25 )
            call jswrit( xl1, yl1 )
            call jschsz( .3 )
            xl2( 1 ) = xl1 - .04 * xmax
            yl2( 1 ) = yl1 + .01 * ymax
            xl2( 2 ) = xl2( 1 ) + .03 * xmax
            yl2( 2 ) = yl2( 1 )
            call jscoin( xl2, yl2, 2, krun )
          end if
          call jscoin( t, x, n, krun )
        end if
c switch to next run or original if done
        krun = krun + 1
        if( krun .le. nruns )then
          call cmprun( krun )
          go to 5
        else
          call cmprun( 1 )
          krun = 1
        end if
c choose option
        call gtchar( 'Select option: com, energy, lmom, amom, acc, ' //
     +               'top, virial, escapers, dht or end', opt )
      end do
      return
      end
