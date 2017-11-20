      program isotropy
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c program to produce a graphic comparison of departures from Kepler of the
c   radial variation of the acceleration computed using the grid of field
c   particles on shells or rings at fixed distances from a unit mass source
c   particle.  The location of the source particle is chosen randomly a
c   given number of times.
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local arays
      integer mrad
      parameter ( mrad = 50 )
      integer ksource, msource, mfield( mrad )
      real*8 fr( mrad ), fr2( mrad ), ft2( mrad ), rad( mrad )
      real*8 ftot(mrad), frmax( mrad ), frmin( mrad )
      real*8 phi( mrad ), phi2( mrad )
      common / variance / msource, ksource, mfield, rad, fr, fr2, ft2,
     +                    ftot, frmax, frmin, phi, phi2
c
c externals
      external sforce, sphi
      real sphi
c
c local variables
      integer ir, kk
      real x, xmax, xx, y, ymax, ymin, yy, zz
c
      call setnag
c initialize
      call getset
      call boilrp( .false. )
c set parallel flag
      parallel = numprocs .gt. 1
      if( parallel )call crash( 'MAIN', 'single processor version' )
c restart .lis file
      rewind no
c read .dat file
      call msetup
      call grdset
      call setrun
      close( ni )
c reset particle mass
      pmass = 1
c create or check Greens function
      call greenm( .false. )
      lprint = .false.
c reserve space
      call setspc
      call switch( 1 )
c initialize results arrays
      do ir = 1, mrad
c distances in grid units
        if( p2d .or. p3d )then
          rad( ir ) = .1 * softl * real( ir )
        else
          rad( ir ) = .2 * real( ir )
        end if
        fr( ir ) = 0
        ftot( ir ) = 0
        fr2( ir ) = 0
        ft2( ir ) = 0
        phi( ir ) = 0
        phi2( ir ) = 0
        frmin( ir ) = 100
        frmax( ir ) = -100
        mfield( ir ) = 0
      end do
      phys = .true.
c number of positions for source mass
      call gtintg( 'Enter number of source masses', msource )
      msource = max( msource, 1 )
      msource = min( msource, 1000 )
      do kk = 1, msource
c place a particle at x,y,z position on the grid
        call place( xx, yy, zz )
        call masset
c calculate field
        call findf( .true. )
c examine field
        call examine
        if( mod( kk, 20 ) .eq. 0 )print *, 'done', kk
      end do
c compute means and dispersions
      do ir = 1, mrad
        if( mfield( ir ) .gt. 0 )then
          fr( ir ) = fr( ir ) / real( mfield( ir ) )
          ftot( ir ) = ftot( ir ) / real( mfield( ir ) )
          fr2( ir ) =
     +         max( fr2( ir ) - mfield( ir ) * fr( ir )**2, 0.d0 )
          fr2( ir ) = sqrt( fr2( ir ) / real( mfield( ir ) ) )
          ft2( ir ) = max( ft2( ir ), 0.d0 )
          ft2( ir ) = sqrt( ft2( ir ) / real( mfield( ir ) ) )
          Phi( ir ) = Phi( ir ) / real( mfield( ir ) )
          Phi2( ir ) =
     +         max( Phi2( ir ) - mfield( ir ) * Phi( ir )**2, 0.d0 )
          Phi2( ir ) = sqrt( Phi2( ir ) / real( mfield( ir ) ) )
        end if
        fr( ir ) = fr( ir ) * rad( ir )**2
        ftot( ir ) = ftot( ir ) * rad( ir )**2
        fr2( ir ) = fr2( ir ) * rad( ir )**2
        ft2( ir ) = ft2( ir ) * rad( ir )**2
        frmin( ir ) = frmin( ir ) * rad( ir )**2
        frmax( ir ) = frmax( ir ) * rad( ir )**2
      end do
c start plot
      call jsbgn
c force plot
      call jspage
      ymax = 1.5
      xmax = rad( mrad )
      call jscale( 0., xmax, -.3 * ymax, ymax )
      call jsaxis( 'x', '\\xi', 1 )
      call jsaxis( 'y', '-\\xi^{2}f', 1 )
      call jsbldt( 'Grid =' )
      if( c3d )then
        call jsbldi( ngx, 3 )
        call jsbldt( 'x' )
        call jsbldi( ngy, 3 )
        call jsbldt( 'x' )
        call jsbldi( ngz, 3 )
      else if( c2d )then
        call jsbldi( ngx, 3 )
        call jsbldt( 'x' )
        call jsbldi( ngy, 3 )
        call jsbldt( '\\epsilon =' )
      else if( p2d .or. p3d )then
        call jsbldi( nr( 1 ), 3 )
        call jsbldt( 'x' )
        call jsbldi( na, 3 )
        if( p3d )then
          call jsbldt( 'x' )
          call jsbldi( ngz, 3 )
        end if
      end if
      call jswrit( .1 * xmax, .9 * ymax )
      if( c3d )then
c compare with Monaghan kernel
        tsoft = 2
        softl = 1.8
        softl2 = softl**2
        call jsbldt( 'effective' )
      end if
      call jsbldt( '\\epsilon =' )
      call jsbldf( softl, 6, 3 )
      call jswrit( .1 * xmax, .8 * ymax )
c radial force
      do ir = 1, mrad
        x = rad( ir )
        y = fr( ir )
        if( ir .eq. 1 )call jsmove( x, y )
        call jsline( x, y )
      end do
c non-radial force
      do ir = 1, mrad
        x = rad( ir )
        y = fr( ir ) + fr2( ir )
        call jsmove( x, y )
        y = fr( ir ) - fr2( ir )
        call jsline( x, y )
        x = x + .01
        y = ft2( ir )
        call jsmove( x, y )
        y = -ft2( ir )
        call jsline( x, y )
      end do
c plot force from a softened mass
      call jsdash( 2, 2, 2, 2 )
      call jsplot( sforce )
      call jsdash( 0, 0, 0, 0 )
c potentials
      call jspage
      ymin = -1.1
      if( .not. c3d )ymin = 1.1 * sphi( 0. )
      call jscale( 0., xmax, ymin, 0. )
      call jsaxis( 'x', '\\xi', 1 )
      call jsaxis( 'y', '\\Phi(\\xi)', 1 )
      call jsbldt( 'Grid =' )
      if( c3d )then
        call jsbldi( ngx, 3 )
        call jsbldt( 'x' )
        call jsbldi( ngy, 3 )
        call jsbldt( 'x' )
        call jsbldi( ngz, 3 )
      else if( c2d )then
        call jsbldi( ngx, 3 )
        call jsbldt( 'x' )
        call jsbldi( ngy, 3 )
      else if( p2d .or. p3d )then
        call jsbldi( nr( 1 ), 3 )
        call jsbldt( 'x' )
        call jsbldi( na, 3 )
        if( p3d )then
          call jsbldt( 'x' )
          call jsbldi( ngz, 3 )
        end if
      end if
      call jswrit( .1 * xmax, .1 * ymin )
c plot potentials
      do ir = 1, mrad
        x = rad( ir )
        y = Phi( ir )
        if( ir .eq. 1 )call jsmove( x, y )
        call jsline( x, y )
      end do
      do ir = 1, mrad
        x = rad( ir )
        y = Phi( ir ) + Phi2( ir )
        call jsmove( x, y )
        y = Phi( ir ) - Phi2( ir )
        call jsline( x, y )
      end do
c plot potential function of softening kernel
      call jsdash( 2, 2, 2, 2 )
      call jsplot( sphi )
      call jsdash( 0, 0, 0, 0 )
      call jsend
      end

      real function sforce( r )
      implicit none
      real r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c external
      real softm
c
c local variable
      real rs
c
      if( tsoft .eq. 1 )then
        sforce = r**3 / ( r**2 + softl2 )**1.5
      else if( tsoft .eq. 2 )then
        rs = r / softl
        sforce = 0
        if( r .gt. 0. )sforce = softm( rs )
      else
        call crash( 'sforce', 'unrecognized softening rule' )
      end if
      return
      end

      real function sphi( r )
      implicit none
      real r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c external
      real sftpot
c
      sphi = sftpot( r * r )
      return
      end

      subroutine place( xx, yy, zz )
      implicit none
c
c calling arguments
      real xx, yy, zz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c externals
      logical offgrd
      real*8 ranuni
c
c local variable
      integer i
      real scale
c
      jgrid = 1
      call switch( 1 )
c initialize list pointers
      do i = 1, nlists
        islist( 2, i, 1 ) = -1
      end do
      scale = 2
      if( p2d .or. p3d )scale = 2 * rgrid( 1 )
c generate and save a particle at a random sub-cell position
    1 newc( 1, 1 ) = scale * ( ranuni( 0. ) - .5 )
      newc( 2, 1 ) = scale * ( ranuni( 0. ) - .5 )
      xx = newc( 1, 1 )
      yy = newc( 2, 1 )
      zz = 0
      if( threed )then
        newc( 3, 1 ) = zm( 1 ) * ( ranuni( 0. ) - .5 )
        zz = newc( 3, 1 )
      end if
      if( ( p2d .or. p3d ) .and. 
     +    ( sqrt( xx**2 + yy**2 ) .gt. .1 * rgrid( 1 ) ) )go to 1
      if( offgrd( 1 ) )go to 1
c set labels
      iz( 1 ) = 1
      iflag( 1 ) = 1
      loc( 1 ) = 0
      pwt( 1 ) = 1
      label( 1 ) = 1
      call relabl( 1 )
c save particle
      call scattr( 1 )
c update list table
      do i = 1, nlists
        islist( 1, i, 1 ) = islist( 2, i, 1 )
      end do
      return
      end 

      subroutine examine
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      integer mrad
      parameter ( mrad = 50 )
      integer msource, ksource, mfield( mrad )
      real*8 fr( mrad ), fr2( mrad ), ft2( mrad ), rad( mrad )
      real*8 ftot(mrad), frmax( mrad ), frmin( mrad )
      real*8 phi( mrad ), phi2( mrad )
      common / variance / msource, ksource, mfield, rad, fr, fr2, ft2,
     +                    ftot, frmax, frmin, phi, phi2
c
c externals
      logical offgrd, ongrd1
      real*8 ranuni
c
c local variables
      integer i, ir, j, m, mfld
      parameter ( mfld = 100 )
      real*8 dx, dy, dz, frad, ftot2, r, x, y, z, theta1, phi1, z1
      include 'inc/pi.f'
c
c get position of source mass
      j = 0
      do ilist = 1, nlists
        call interpret
        inext = islist( 1, ilist, 1 )
        if( inext .ge. 0 )then
          call gather( i )
          j = j + i
          if( j .ne. 1 )call
     +                   crash( 'EXAMINE', 'Wrong number of particles' )
          x = oldc( 1, 1 )
          y = oldc( 2, 1 )
          if( threed )z = oldc( 3, 1 )
c work over radii
          do ir = 1, mrad
            r = rad( ir )
c generate mfld particles randomly on a sphere or circle
            do i = 1, mfld
              if( threed )then
                z1 = ranuni( 0. ) - .5
                phi1 = 2.0 * pi * ranuni( 0. )
                theta1 = asin( z1 )
                newc( 1, i ) = x + r * cos( theta1 ) * cos( phi1 )
                newc( 2, i ) = y + r * cos( theta1 ) * sin( phi1 )
                newc( 3, i ) = z + r * z1
              else
                phi1 = 2.0 * pi * ranuni( 0. )
                newc( 1, i ) = x + r * cos( phi1 )
                newc( 2, i ) = y + r * sin( phi1 )
              end if
              iz( i ) = 1
              if( offgrd( i ) )iz( i ) = nzones
              iflag( i ) = 1
              if( hybrid )then
                if( ongrd1( i ) )then
                  label( i ) = 1
                else
                  label( i ) = 2
                end if
              else
                label( i ) = 1
              end if
              oldc( 1, i ) = newc( 1, i )
              oldc( 2, i ) = newc( 2, i )
              if( threed )oldc( 3, i ) = newc( 3, i )
              call relabl( i )
            end do
            iplane = -1
            call getacc( mfld )
c average round ring or sphere
            m = 0
            do i = 1, mfld
              if( nskip( i ) )then
                mfield( ir ) = mfield( ir ) + 1
                dx = x - oldc( 1, i )
                dy = y - oldc( 2, i )
                frad = dx * acc( 1, i ) + dy * acc( 2, i ) 
                ftot2 = acc( 1, i )**2 + acc( 2, i )**2
                if( threed )then
                  dz = z - oldc( 3, i )
                  frad = frad + dz * acc( 3, i )
                  ftot2 = ftot2 + acc( 3, i )**2
                end if
                frad = frad / r
                frmin( ir ) = min( frad, frmin( ir ) )
                frmax( ir ) = max( frad, frmax( ir ) )
                fr( ir ) = fr( ir ) + frad
                fr2( ir ) = fr2( ir ) + frad**2
                ft2( ir ) = ft2( ir ) + ftot2 - frad**2
                ftot( ir ) = ftot( ir ) + sqrt( ftot2 )
                phi( ir ) = phi( ir ) + gpot( i )
                phi2( ir ) = phi2( ir ) + gpot( i )**2
              end if
            end do
          end do
        end if
      end do
      return
      end

      include 'inc/pgiden.f'
