      program exptst
      use aarrays
      implicit none
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
      integer itype
      real fldr, srcr
      common / local / fldr, srcr, itype
c
c external
      external frad, fperp, phig
      real frad, fperp, phig, roundup
      real*8 ranuni
c
c local variables
      integer i, j, maxs3l
      real e, f, r, s, x, xmax, xmin, y, ymax, ymin, z
      include 'inc/pi.f'

      integer k
      real grofu, uofgr

c
      call setnag
c read run number from standard input and open main ASCII I/O files
      call getset
c tidy up if there was a previous .lis file
      rewind no
c read .dat file (ASCII input)
      call msetup
      call grdset
      call setrun
      close( ni )
      maxs3l = s3lmax
c allocate space
      call setspc
c reset particle mass
      pmass = 1
c place one particle at r=1
      x = .7
      y = 0
      z = sqrt( 1. - x* x )
c
      srcr = -1
      r = s3rad( nr( 1 ) )
      do while ( srcr .lt. 0. .or. srcr .gt. r )
        call gtreal( 'Enter source rad', srcr )
      end do
      x = srcr * x
      y = srcr * y
      z = srcr * z
      call place( x, y, z )
      istep = 0
      isfld = istep
      isubst = 1
c get radius of field point
      fldr = -1
      do while ( fldr .lt. 0. .or. fldr .gt. 1.05 * r )
        call gtreal( 'Enter fldr', fldr )
      end do
c
      call jsbgn
      itype = 1
c work over acceleration then potential
      do j = 1, 3, 2
        if( itype .eq. 1 )then
          if( j .eq. 1 )then
            f = 1. / ( srcr - fldr )**2
            if( srcr .gt. fldr )then
              ymax = roundup( f )
              ymin = -.1 * ymax
            else
              ymin = -roundup( f )
              ymax = -.1 * ymin
            end if
          else if( j .eq. 2 )then
            f = 1. / ( srcr - fldr )**2
            ymax = roundup( f )
            ymin = -.1 * ymax            
          else
            ymin = 1. / abs( srcr - fldr )
            ymin = -roundup( ymin )
            ymax = -.1 * ymin
          end if
        else
          ymax = 0
          ymin = -250
          if( j .eq. 1 )ymin = -5000
          xmin = 1.e-3 * srcr
          xmax = 1.5 * srcr
        end if
c initialize plot
        call jspage
        if( itype .eq. 1 )then
          xmax = pi
          xmin = -xmax
          call jscale( xmin, xmax, ymin, ymax )
          call jsaxis( 'x', '\\phi', 1 )
        else
          call jscale( xmin, xmax, ymin, ymax )
          call jsaxis( 'x', 'Field r', 1 )
        end if
        if( j .eq. 1 )then
          call jsaxis( 'y', 'f_r', 1 )
        else if( j .eq. 2 )then
          call jsaxis( 'y', 'f_{perp}', 1 )
        else
          call jsaxis( 'y', '\\Phi', 1 )
        end if
c note radii of source and field points
        call jsbldt( 'Src r=' )
        call jsbldf( srcr, 5, 3 )
        if( ymax .gt. -ymin )then
          y = .9 * ymax + .1 * ymin
        else
          y = .1 * ymax + .9 * ymin
        end if
        call jswrit( .95 * xmin + .05 * xmax, y )
        if( itype .eq. 1 )then
          call jsbldt( 'Fld r=' )
          call jsbldf( fldr, 5, 3 )
          call jswrit( .4 * xmin + .6 * xmax, y )
        end if
c work over lmax
        do i = maxs3l, 0, -1
          s3lmax = i
          if( s3lmax .le. s3maxl )then
            s3ntm = ( ( s3lmax + 1 ) * ( s3lmax + 2 ) ) / 2
c set expansion coeffs
            call masset
c save expansion value at phi = 0
            if( i .eq. maxs3l )then
              if( itype .eq. 1 )xmin = 0
              if( j .eq. 1 )then
                e = frad( xmin )
              else if( j .eq. 2 )then
                e = fperp( xmin )
              else
                e = phig( xmin )
              end if
            end if
            if( j .eq. 1 )then
              call jsplot( frad )
            else if( j .eq. 2 )then
              if( s3lmax .gt. 0 )call jsplot( fperp )
            else
              call jsplot( phig )
            end if
          end if
        end do
c plot exact value
        s3lmax = -1
        call pgsci( 2 )
        call jsdash( 2, 2, 2, 2 )
        if( j .eq. 1 )then
          call jsplot( frad )
        else if( j .eq. 2 )then
          call jsplot( fperp )
        else
          call jsplot( phig )
        end if
        call jsdash( 0, 0, 0, 0 )
        call pgsci( 1 )
c save exact value at phi = 0
        if( j .eq. 1 )then
          f = frad( xmin )
        else if( j .eq. 2 )then
          f = 1! fperp( xmin )
        else
          f = phig( xmin )
        end if
        if( abs( f ) .gt. 1.e-6 )then
          r = abs( ( e - f ) / f )
          if( j .eq. 1 )then
            print *, 'fractional rad force error at phi = 0', r
          else if( j .eq. 3 )then
            print *, 'fractional potential error at phi = 0', r
          end if
        end if
      end do
      call jsend
      stop
      end

      subroutine place( xx, yy, zz )
      use aarrays
      implicit none
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
c
c local variable
      integer i
c
      jgrid = 1
      call switch( 1 )
c initialize list pointers
      do i = 1, nlists
        islist( 2, i, 1 ) = -1
      end do
c set and save this particle
    1 newc( 1, 1 ) = xx
      newc( 2, 1 ) = yy
      if( threed )newc( 3, 1 ) = zz
      if( offgrd( 1 )
     +         )call crash( 'Place', 'Source particle is off the grid' )
c set labels
      iz( 1 ) = 1
      iflag( 1 ) = 1
      loc( 1 ) = 0
      pwt( 1 ) = 1
      label( 1 ) = 1
      call relabl( 1 )
c save particle
      call scattr( 1, lptcls, ptcls )
c update list table
      do i = 1, nlists
        islist( 1, i, 1 ) = islist( 2, i, 1 )
      end do
      return
      end 

      real function frad( phi )
      implicit none
c
      real phi
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
      integer itype
      real fldr, srcr
      common / local / fldr, srcr, itype
c
c externals
      logical offgrd, ongrd1
c
c local variables
      integer i, ir, j, m, mfld
      real*8 r, t, xs( 3 ), phi1, tiny
      parameter ( tiny = 1.d-8 )
      include 'inc/pi.f'
c
      frad = 0
c get position of source mass
      j = 0
      do ilist = 1, nlists
        call interpret
        inext = islist( 1, ilist, 1 )
        if( inext .ge. 0 )then
          call gather( i )
          j = j + i
          if( j .ne. 1
     +             )call crash( 'EXAMINE', 'Wrong number of particles' )
c save source position
          do i = 1, 3
            xs( i ) = oldc( i, 1 )
          end do
c create a field particle on the circle
          r = sqrt( xs( 1 )**2 + xs( 2 )**2 )
          if( itype .eq. 1 )then
            phi1 = phi + atan2( xs( 2 ), xs( 1 ) )
            oldc( 1, 1 ) = fldr * r * cos( phi1 ) / ( srcr + tiny )
            oldc( 2, 1 ) = fldr * r * sin( phi1 ) / ( srcr + tiny )
            oldc( 3, 1 ) = fldr * xs( 3 ) / ( srcr + tiny )
          else
c create a field particle on the radius vector through the source
            phi1 = atan2( xs( 2 ), xs( 1 ) ) + .1
            oldc( 1, 1 ) = phi * r * cos( phi1 ) / ( srcr + tiny )
            oldc( 2, 1 ) = phi * r * sin( phi1 ) / ( srcr + tiny )
            oldc( 3, 1 ) = phi * xs( 3 ) / ( srcr + tiny )
          end if
c need newc for these externals
          do i = 1, 3
            newc( i, 1 ) = oldc( i, 1 )
          end do
          iz( 1 ) = 1
          if( offgrd( 1 ) )iz( 1 ) = nzones
          iflag( 1 ) = 1
          if( hybrid )then
            if( ongrd1( 1 ) )then
              label( 1 ) = 1
            else
              label( 1 ) = 2
            end if
          else
            label( 1 ) = 1
          end if
          call relabl( 1 )
c restore original position of source mass
          do i = 1, 3
            newc( i, 1 ) = xs( i )
          end do
c compute radial force - oldc is the field position, newc the source
          if( s3lmax .lt. 0 )then
c exact expression
            r = 0
            do i = 1, 3
              r = r + ( newc( i, 1 ) - oldc( i, 1 ) )**2
            end do
            t = 1. / ( r + tiny )**1.5
            do i = 1, 3
              acc( i, 1 ) = t * ( newc( i, 1 ) - oldc( i, 1 ) )
            end do
          else
c get grid force components
            call getacc( 1 )
          end if
          frad = 0
          do i = 1, 3
            frad = frad + acc( i, 1 ) * oldc( i, 1 )
          end do
          frad = frad / ( fldr + tiny )
        end if
      end do
      return
      end

      real function fperp( phi )
      implicit none
c
      real phi
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
      integer itype
      real fldr, srcr
      common / local / fldr, srcr, itype
c
c externals
      logical offgrd, ongrd1
c
c local variables
      integer i, ir, j, m, mfld
      real*8 frad, ft, r, t, xs( 3 ), phi1, tiny
      parameter ( tiny = 1.d-8 )
      include 'inc/pi.f'
c
      fperp = 0
c get position of source mass
      j = 0
      do ilist = 1, nlists
        call interpret
        inext = islist( 1, ilist, 1 )
        if( inext .ge. 0 )then
          call gather( i )
          j = j + i
          if( j .ne. 1
     +             )call crash( 'EXAMINE', 'Wrong number of particles' )
c save source position
          do i = 1, 3
            xs( i ) = oldc( i, 1 )
          end do
c create a field particle on the circle
          r = sqrt( xs( 1 )**2 + xs( 2 )**2 )
          if( itype .eq. 1 )then
            phi1 = phi + atan2( xs( 2 ), xs( 1 ) )
            oldc( 1, 1 ) = fldr * r * cos( phi1 ) / ( srcr + tiny )
            oldc( 2, 1 ) = fldr * r * sin( phi1 ) / ( srcr + tiny )
            oldc( 3, 1 ) = fldr * xs( 3 ) / ( srcr + tiny )
          else
c create a field particle on the radius vector through the source
            phi1 = atan2( xs( 2 ), xs( 1 ) ) + .1
            oldc( 1, 1 ) = phi * r * cos( phi1 ) / ( srcr + tiny )
            oldc( 2, 1 ) = phi * r * sin( phi1 ) / ( srcr + tiny )
            oldc( 3, 1 ) = phi * xs( 3 ) / ( srcr + tiny )
          end if
c need newc for these externals
          do i = 1, 3
            newc( i, 1 ) = oldc( i, 1 )
          end do
          iz( 1 ) = 1
          if( offgrd( 1 ) )iz( 1 ) = nzones
          iflag( 1 ) = 1
          if( hybrid )then
            if( ongrd1( 1 ) )then
              label( 1 ) = 1
            else
              label( 1 ) = 2
            end if
          else
            label( 1 ) = 1
          end if
          call relabl( 1 )
c restore original position of source mass
          do i = 1, 3
            newc( i, 1 ) = xs( i )
          end do
c compute radial force - oldc is the field position, newc the source
          if( s3lmax .lt. 0 )then
c exact expression
            r = 0
            do i = 1, 3
              r = r + ( newc( i, 1 ) - oldc( i, 1 ) )**2
            end do
            t = 1. / ( r + tiny )**1.5
            do i = 1, 3
              acc( i, 1 ) = t * ( newc( i, 1 ) - oldc( i, 1 ) )
            end do
          else
c get grid force components
            call getacc( 1 )
          end if
          ft = 0
          do i = 1, 3
            ft = ft + acc( i, 1 )**2
          end do
          frad = 0
          do i = 1, 3
            frad = frad + acc( i, 1 ) * oldc( i, 1 )
          end do
          frad = frad / ( fldr + tiny )
          fperp = sqrt( max( ft - frad**2, 0.d0 ) )
        end if
      end do
      return
      end

      real function phig( phi )
      implicit none
c
      real phi
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
      integer itype
      real fldr, srcr
      common / local / fldr, srcr, itype
c
c externals
      logical offgrd, ongrd1
c
c local variables
      integer i, ir, j, m, mfld
      real*8 r, t, x, y, z, phi1, tiny
      parameter ( tiny = 1.d-6 )
      include 'inc/pi.f'
c
      phig = 0

      if( srcr .eq. 0. )then
        phig = -1. / fldr
        return
      end if

c get position of source mass
      j = 0
      do ilist = 1, nlists
        call interpret
        inext = islist( 1, ilist, 1 )
        if( inext .ge. 0 )then
          call gather( i )
          j = j + i
          if( j .ne. 1
     +             )call crash( 'EXAMINE', 'Wrong number of particles' )
          x = oldc( 1, 1 )
          y = oldc( 2, 1 )
          if( threed )z = oldc( 3, 1 )
          r = sqrt( x * x + y * y )
c create a field particle on the circle
          if( itype .eq. 1 )then
            phi1 = atan2( y, x ) + phi
            oldc( 1, 1 ) = fldr * r * cos( phi1 ) / ( srcr + tiny )
            oldc( 2, 1 ) = fldr * r * sin( phi1 ) / ( srcr + tiny )
            oldc( 3, 1 ) = fldr * z / ( srcr + tiny )
          else
c create a field particle on the radius vector through the source
            phi1 = atan2( y, x ) + .1
            oldc( 1, 1 ) = phi * r * cos( phi1 ) / ( srcr + tiny )
            oldc( 2, 1 ) = phi * r * sin( phi1 ) / ( srcr + tiny )
            oldc( 3, 1 ) = phi * z / ( srcr + tiny )
          end if
          do i = 1, 3
            newc( i, 1 ) = oldc( i, 1 )
          end do
          iz( 1 ) = 1
          if( offgrd( 1 ) )iz( 1 ) = nzones
          iflag( 1 ) = 1
          if( hybrid )then
            if( ongrd1( 1 ) )then
              label( 1 ) = 1
            else
              label( 1 ) = 2
            end if
          else
            label( 1 ) = 1
          end if
          call relabl( 1 )
c restore original position of source mass
          newc( 1, 1 ) = x
          newc( 2, 1 ) = y
          if( threed )newc( 3, 1 ) = z
c compute potential
          if( s3lmax .lt. 0 )then
            r = 0
            do i = 1, 3
              r = r + ( newc( i, 1 ) - oldc( i, 1 ) )**2
            end do
            r = sqrt( r ) + tiny
            phig = -1 / r
          else
            call getacc( 1 )
            phig = gpot( 1 )
          end if
        end if
      end do
      return
      end

      include 'inc/pgiden.f'

      subroutine ncheck
      implicit none
      return
      end
