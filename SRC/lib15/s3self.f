      subroutine s3self( is, s )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c Routine to determine the contribution from one particle to the total mass
c   on the spherical grid - which was assigned from its last position stored
c   in /ptcls/.  It could contribute to both the interior and exterior grid
c   points, but those needed to null its self-force enter in different ways
c   depending on whether it is now inside or outside the nearest grid point,
c   or in a completely different shell.
c
c calling arguments
      integer is
      real s( 2, 4, 1 )
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
c external
      real grofu, uofgr
c
c local arrays
      real cmp( 0:s3maxl ), smp( 0:s3maxl )
c
c local variables
      integer i, k, l, lm, m
      real arg, cp, du, rad, rc, rfac, rl, sp, u, w
      real w1, w2, w3, w4, x, y, z
c
      if( s3lmax .gt. s3maxl )call crash( 'S3SELF',
     +                                        'local arrays too small' )
c just one particle - use original position
      x = oldpos( 1, is ) - xcen( 1, 1 )
      y = oldpos( 2, is ) - xcen( 2, 1 )
      z = oldpos( 3, is ) - xcen( 3, 1 )
      rad = sqrt( x * x + y * y + z * z )
      u = uofgr( rad )
      k = u + 1.5
      if( k .eq. 1 )then
        w1 = 0
        w2 = 1
      else
        du = u - real( k - 1 )
        w1 = .5 - du
        w2 = .5 + du
        w3 = grofu( .5 * ( real( k - 1 ) + u - .5 ) )
        w4 = grofu( .5 * ( real( k - 1 ) + u + .5 ) )
      end if
c make table of cos(-m\phi) and sin(-m\phi)
      rc = sqrt( x**2 + y**2 ) + 1.e-8
      cp = x / rc
      sp = y / rc
      cmp( 0 ) = 1
      smp( 0 ) = 0
      do m = 1, s3lmax
        cmp( m ) = cmp( m - 1 ) * cp + smp( m - 1 ) * sp
        smp( m ) = smp( m - 1 ) * cp - cmp( m - 1 ) * sp
      end do
c get Plms
      arg = z / rad
      call s3tplm( dble( arg ), .false. )
      if( ( k .gt. 1 ) .and. ( k .lt. nr( jgrid ) ) )then
c fragments of particle's mass interior and exterior to current ring
        do i = 1, 4
          if( i .eq. 1 )then
            rl = 1. / s3rad( k )
            rfac = rl * w3
            w = w1
          else if( i .eq. 2 )then
            rl = 1. / s3rad( k + 1 )
            rfac = rl * w4
            w = w2
          else if( i .eq. 3 )then
            rl = 1. / w3
            rfac = rl * s3rad( k - 1 )
            w = w1
          else
            rl = 1. / w4
            rfac = rl * s3rad( k )
            w = w2
          end if
          lm = 0
          do l = 0, s3lmax
            if( l .gt. 0 )rl = rl * rfac
c surface harmonics selection
            if( lg( l + 1 ) )then
              lm = lm + l + 1
            else
              do m = 0, l
                lm = lm + 1
                s( 1, i, lm ) = w * cmp( m ) * plm( lm ) * rl
                s( 2, i, lm ) = w * smp( m ) * plm( lm ) * rl
              end do
            end if
          end do
        end do
      else if( k .eq. 1 )then
c no mass can be interior to centre
        lm = 0
        rl = w2 / s3rad( 2 )
        rfac = rad / s3rad( 2 )
        do l = 0, s3lmax
          if( l .gt. 0 )rl = rl * rfac
          do m = 0, l
            lm = lm + 1
c surface harmonics selection
            if( .not. lg( m + 1 ) )then
              s( 1, 2, lm ) = cmp( m ) * plm( lm ) * rl
              s( 2, 2, lm ) = smp( m ) * plm( lm ) * rl
            end if
          end do
        end do
c monopole term only for mass exterior to central point
        if( .not. lg( 1 ) )then
          s( 1, 1, 1 ) = 1. / rad
        end if
      else
c only a fraction of the particle is interior to outer ring
        do i = 1, 2
          if( i .eq. 1 )then
            rl = 1. / s3rad( nr( jgrid ) )
            rfac = rl * w3
          else
            rl = 1. / w3
            rfac = rl * s3rad( nr( jgrid ) - 1 )
          end if
          lm = 0
          do l = 0, s3lmax
            rl = rl * rfac
            do m = 0, l
              lm = lm + 1
c apply Fourier filter
              if( .not. lg( m + 1 ) )then
                s( 1, i, lm ) = w1 * cmp( m ) * plm( lm ) * rl
                s( 2, i, lm ) = w1 * smp( m ) * plm( lm ) * rl
              end if
            end do
          end do
        end do
      end if
c combine terms - interior mass fragments
      lm = 0
      rl = s3rad( k ) / s3rad( k + 1 )
      rfac = rl
      if( ncl( is ) .gt. k )then
        w1 = s3rad( k + 1 ) / s3rad( ncl( is ) )
        w2 = s3rad( k + 1 ) / s3rad( ncl( is ) + 1 )
        w3 = w1
        w4 = w2
      end if
      do l = 0, s3lmax
        do m = 0, l
          lm = lm + 1
          if( ncl( is ) .gt. k )then
            s( 1, 2, lm ) = s( 1, 2, lm ) + s( 1, 1, lm ) * rfac
            s( 2, 2, lm ) = s( 2, 2, lm ) + s( 2, 1, lm ) * rfac
            s( 1, 1, lm ) = s( 1, 2, lm ) * w3
            s( 2, 1, lm ) = s( 2, 2, lm ) * w3
            s( 1, 2, lm ) = s( 1, 2, lm ) * w4
            s( 2, 2, lm ) = s( 2, 2, lm ) * w4
          else if( ncl( is ) .eq. k )then
            s( 1, 2, lm ) = s( 1, 2, lm ) + s( 1, 1, lm ) * rfac
            s( 2, 2, lm ) = s( 2, 2, lm ) + s( 2, 1, lm ) * rfac
          else if( ncl( is ) .eq. k - 1 )then
            s( 1, 2, lm ) = s( 1, 1, lm )
            s( 2, 2, lm ) = s( 2, 1, lm )
            s( 1, 1, lm ) = 0
            s( 2, 1, lm ) = 0
          else
            s( 1, 1, lm ) = 0
            s( 2, 1, lm ) = 0
            s( 1, 2, lm ) = 0
            s( 2, 2, lm ) = 0
          end if
        end do
        rfac = rl * rfac
        if( ncl( is ) .gt. k )then
          w3 = w3 * w1
          w4 = w4 * w2
        end if
      end do
c combine terms - exterior mass fragments
      lm = 0
      rl = 0
      if( k .gt. 1 )rl = s3rad( k - 1 ) / s3rad( k )
      rfac = 1
      if( ncl( is ) .lt. k - 1 )then
        w1 = s3rad( ncl( is ) ) / s3rad( k - 1 )
        w2 = s3rad( ncl( is ) + 1 ) / s3rad( k - 1 )
        w3 = 1
        w4 = 1
      end if
      do l = 0, s3lmax
        do m = 0, l
          lm = lm + 1
          if( ncl( is ) .eq. k )then
            s( 1, 3, lm ) = s( 1, 4, lm )
            s( 2, 3, lm ) = s( 2, 4, lm )
            s( 1, 4, lm ) = 0
            s( 2, 4, lm ) = 0
          else if( ncl( is ) .eq. k - 1 )then
            s( 1, 3, lm ) = s( 1, 3, lm ) + s( 1, 4, lm ) * rfac
            s( 2, 3, lm ) = s( 2, 3, lm ) + s( 2, 4, lm ) * rfac
          else if( ncl( is ) .lt. k - 1 )then
            s( 1, 3, lm ) = s( 1, 3, lm ) + s( 1, 4, lm ) * rfac
            s( 2, 3, lm ) = s( 2, 3, lm ) + s( 2, 4, lm ) * rfac
            s( 1, 4, lm ) = s( 1, 3, lm ) * w4
            s( 2, 4, lm ) = s( 2, 3, lm ) * w4
            s( 1, 3, lm ) = s( 1, 3, lm ) * w3
            s( 2, 3, lm ) = s( 2, 3, lm ) * w3
          else
            s( 1, 3, lm ) = 0
            s( 2, 3, lm ) = 0
            s( 1, 4, lm ) = 0
            s( 2, 4, lm ) = 0
          end if
        end do
        rfac = rl * rfac
        if( ncl( is ) .lt. k - 1 )then
          w3 = w3 * w1
          w4 = w4 * w2
        end if
      end do
c scale for mass of particle
      do lm = 1, s3ntm
        do i = 1, 4
          s( 1, i, lm ) = pmass * s( 1, i, lm )
          s( 2, i, lm ) = pmass * s( 2, i, lm )
        end do
      end do
      return
      end
