      subroutine posplt( x, y, z, n, first, rnew, znew )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Plots positions or velocities of the particles
c
c Called from PNTPLT & DMPLOT
c Graphics routines are from JSPLOT
c
c calling arguments
      integer n
      logical first
      real x( n ), y( n ), z( n ), rnew, znew
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c external
      real grofu
c
c local variables
      logical tkl
      real a, rmn, rmx, smx, xmx, x1, x2, x3, x4, ymx, y1, y2, y3, y4
      real zmx
c
      rmn = 0
      rmx = rnew
      tkl = .false.
      if( p2d .or. sf2d )then
        if( uoffs .ne. 0. )rmn = rinh( jgrid ) / lscale
      else
c other grid types
        if( p3a )rmn = grofu( 0. ) / lscale
        xmx = rmx
        ymx = rmx
        zmx = znew
        if( p3d )rmx = max( xmx, zmx )
      end if
c up set frame boundaries
      if( threed )then
        smx = max( xmx + zmx, ymx + zmx, 1.2 * xmx )
        smx = .85 / smx
        x1 = .1
        x2 = x1 + smx * xmx
        x3 = x2
        x4 = x3 + smx * zmx
        y1 = .1
        y2 = y1 + smx * ymx
        y3 = y2
        y4 = y3 + smx * zmx
        tkl = zmx .gt. 1.
      end if
c plot particles
      if( twod )then
        call jssize( .05, .95, 0., .9 )
        call jsescl( -rmx, rmx, -rmx, rmx )
        if( first )then
          if( c2d )then
            call jssize( .05, .95, .05, .95 )
            call jschsz( .3 )
            if( scale_set )then
              a = rmx * unit_L
              call jsescl( -a, a, -a, a )
              call jsaxis( 'x', 'x (kpc)', 1 )
              call jsaxis( 'y', 'y (kpc)', 1 )
              call jsescl( -rmx, rmx, -rmx, rmx )
            else
              call jsescl( -rmx, rmx, -rmx, rmx )
              call jsaxis( 'x', 'x', 1 )
              call jsaxis( 'y', 'y', 1 )
            end if
            call jschsz( .2 )
          else
            if( rmn .gt. 0. )call jscirc( 0., 0., rmn )
            call jscirc( 0., 0., rmx )
          end if
        end if
        call jsymbl( x, y, n, 0 )
      else
c 3-D
        call jssize( x1, x2, y1, y2 )
        call jscale( -xmx, xmx, -ymx, ymx )
        if( first )then
          if( scale_set )then
            call jscale( -xmx * unit_L, xmx * unit_L,
     +                   -ymx * unit_L, ymx * unit_L )
            call jsaxis( 'x', 'x (kpc)', 1 )
            call jsaxis( 'y', 'y (kpc)', 1 )
            call jscale( -xmx, xmx, -ymx, ymx )
          else
            call jsaxis( 'x', 'x', 1 )
            call jsaxis( 'y', 'y', 1 )
          end if
          if( p3a .or. p3d .or. s3d )call jscirc( 0., 0., xmx )
        end if
        call jsymbl( x, y, n, 0 )
c
        call jssize( x3, x4, y1, y2 )
        call jscale( -zmx, zmx, -ymx, ymx )
        if( first )then
          if( scale_set )then
            call jscale( -zmx * unit_L, zmx * unit_L,
     +                   -ymx * unit_L, ymx * unit_L )
            if( tkl )then
              call jsaxis( 'x', 'z (kpc)', 1 )
            else
              call jsaxis( 'x', 'z (kpc)', 0 )
            end if
            call jsaxis( 'y', ' ', 1 )
            call jscale( -zmx, zmx, -ymx, ymx )
          else
            if( tkl )then
              call jsaxis( 'x', 'z', 1 )
            else
              call jsaxis( 'x', 'z', 0 )
            end if
            call jsaxis( 'y', ' ', 0 )
          end if
          a = zm( 1 ) / lscale
          if( zmx .gt. a )then
            call jsdash( 2, 2, 2, 2 )
            call jsmove( a, -ymx )
            call jsline( a, ymx )
            call jsmove( -a, -ymx )
            call jsline( -a, ymx )
            call jsdash( 0, 0, 0, 0 )
          end if
          if( s3d )call jscirc( 0., 0., xmx )
        end if
        call jsymbl( z, y, n, 0 )
c
        call jssize( x1, x2, y3, y4 )
        call jscale( -xmx, xmx, -zmx, zmx )
        if( first )then
          if( scale_set )then
            call jscale( -xmx * unit_L, xmx * unit_L,
     +                   -zmx * unit_L, zmx * unit_L )
            call jsaxis( 'x', ' ', 1 )
            if( tkl )then
              call jsaxis( 'y', 'z (kpc)', 1 )
            else
              call jsaxis( 'y', 'z (kpc)', 0 )
            end if
            call jscale( -xmx, xmx, -zmx, zmx )
          else
            call jsaxis( 'x', ' ', 0 )
            if( tkl )then
              call jsaxis( 'y', 'z', 1 )
            else
              call jsaxis( 'y', 'z', 0 )
            end if
          end if
          if( zmx .gt. a )then
            call jsdash( 2, 2, 2, 2 )
            call jsmove( -xmx, a )
            call jsline( xmx, a )
            call jsmove( -xmx, -a )
            call jsline( xmx, -a )
            call jsdash( 0, 0, 0, 0 )
          end if
          if( s3d )call jscirc( 0., 0., xmx )
        end if
        call jsymbl( x, z, n, 0 )
        call jssize( x1, x2, y1, y2 )
        call jscale( -xmx, xmx, -ymx, ymx )
      end if
      return
      end
