      subroutine s3dsum( jst )
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to assign mass of particles to the spherical grid.  Each particle
c   contributes to both the interior and exterior grid points
c The shell number and weights should have been determined by a previous
c   call to weight.  ncl( is ) is the number of the largest grid shell just
c   interior to the particle.
      use aarrays
      implicit none
c
c calling arguments
      integer jst
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
c local arrays
      real cmp( 0:s3maxl ), smp( 0:s3maxl )
c
c local variables
      integer i, is, j, jr, jz, k, kg, l, lm, m, nrg
      logical grid2
      real arg, cp, rad, rc, rfac, rl, sp, splm, w, x, y, z
c
      if( s3lmax .gt. s3maxl )call crash( 'S3DSUM',
     +                                        'local arrays too small' )
      grid2 = hybrid .and. ( jlist .eq. 1 )
      nrg = nr( jgrid )
c work over particles
      do is = 1, jst
c skip particles off this grid
        if( nskip( is ) )then
          jz = max( iz( is ), 1 )
          x = newc( 1, is ) - xcpred( 1, jz, jgrid )
          y = newc( 2, is ) - xcpred( 2, jz, jgrid )
          z = newc( 3, is ) - xcpred( 3, jz, jgrid )
          rad = rr( is ) + 1.e-8
          k = ncl( is )
          kg = jgrid
          if( grid2 .and. ( label( is ) .eq. 1 ) )kg = 1
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
          if( ( k .gt. 1 ) .and. ( k .lt. nrg ) )then
c fragments of particle's mass interior and exterior to current ring
            do i = 1, 4
              if( i .eq. 1 )then
                rl = 1. / s3rad( k )
                rfac = rl * wt( 3, is )
                jr = k
                w = wt( 1, is )
                j = 1
              else if( i .eq. 2 )then
                rl = 1. / s3rad( k + 1 )
                rfac = rl * wt( 4, is )
                jr = k + 1
                w = wt( 2, is )
                j = 1
              else if( i .eq. 3 )then
                rl = 1. / wt( 3, is )
                rfac = rl * s3rad( k - 1 )
                jr = k - 1
                w = wt( 1, is )
                j = 3
              else
                rl = 1. / wt( 4, is )
                rfac = rl * s3rad( k )
                jr = k
                w = wt( 2, is )
                j = 3
              end if
              j = 4 * s3ntm * ( jr - 1 ) + j - 4
              lm = 0
              do l = 0, s3lmax
                if( l .gt. 0 )rl = rl * rfac
c surface harmonics selection
                if( lg( l + 1 ) )then
                  j = j + 4 * ( l + 1 )
                  lm = lm + l + 1
                else
                  do m = 0, l
                    j = j + 4
                    lm = lm + 1
                    splm = plm( lm )
                    s3dmss( j, jz, kg ) =
     +                s3dmss( j, jz, kg ) + w * cmp( m ) * splm * rl
                    s3dmss( j + 1, jz, kg ) =
     +            s3dmss( j + 1, jz, kg ) + w * smp( m ) * splm * rl
                  end do
                end if
              end do
            end do
          else if( k .eq. 1 )then
c no mass can be interior to centre
            j = 4 * s3ntm - 4
            lm = 0
            rl = wt( 2, is ) / s3rad( 2 )
            rfac = rad / s3rad( 2 )
            do l = 0, s3lmax
              if( l .gt. 0 )rl = rl * rfac
c surface harmonics selection
              if( lg( l + 1 ) )then
                j = j + 4 * ( l + 1 )
                lm = lm + l + 1
              else
                do m = 0, l
                  j = j + 4
                  lm = lm + 1
                  splm = plm( lm )
                  s3dmss( j + 1, jz, kg ) =
     +                s3dmss( j + 1, jz, kg ) + cmp( m ) * splm * rl
                  s3dmss( j + 2, jz, kg ) =
     +                s3dmss( j + 2, jz, kg ) + smp( m ) * splm * rl
                end do
              end if
            end do
c monopole term only for mass exterior to central point
            if( .not. lg( 1 ) )then
c              s3dmss( 1, jz, kg ) =
c     +               s3dmss( 1, jz, kg ) + wt( 1, is )     ! wt( 1, is ) = 0
              s3dmss( 3, jz, kg ) = s3dmss( 3, jz, kg ) + 1. / rad
            end if
          else
c only a fraction of the particle is interior to outer ring
            do i = 1, 2
              if( i .eq. 1 )then
                rl = 1. / s3rad( nrg )
                rfac = rl * wt( 3, is )
                jr = nrg
                j = 1
              else
                rl = 1. / wt( 3, is )
                rfac = rl * s3rad( nrg - 1 )
                jr = nrg - 1
                j = 3
              end if
              j = 4 * s3ntm * ( jr - 1 ) + j - 4
              lm = 0
              do l = 0, s3lmax
                rl = rl * rfac
c surface harmonics selection
                if( lg( l + 1 ) )then
                  j = j + 4 * ( l + 1 )
                  lm = lm + l + 1
                else
                  do m = 0, l
                    j = j + 4
                    lm = lm + 1
                    splm = plm( lm )
                    s3dmss( j, jz, kg ) = s3dmss( j, jz, kg ) +
     +                            wt( 1, is ) * cmp( m ) * splm * rl
                    s3dmss( j + 1, jz, kg ) = s3dmss( j + 1, jz, kg ) +
     +                            wt( 1, is ) * smp( m ) * splm * rl
                  end do
                end if
              end do
            end do
          end if
        end if
      end do
      return
      end
