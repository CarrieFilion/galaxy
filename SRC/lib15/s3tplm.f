      subroutine s3tplm( x, derivs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to tabulate all the P_l^m values, and derivatives if required,
c   for a single value of the calling argument
c algorithm and notation are lifted from Numerical Recipes pp 182-183
c   modified to include the sqrt( ( l - m )! / ( l + m )! ) factor
c   see notes in herschel:/home/sellwood/docs/codes/s3d/forces.tex
c
c calling arguments
      logical derivs
      real*8 x
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local array
      integer mfac
      parameter ( mfac = 3 * ( ( s3maxl + 1 ) * s3maxl ) / 2 )
      real*8 rfac( mfac )
c
c local variables
      integer i, iuse, j, k, l, m
      real*8 fact, f2, pll, pm, pmm, pmp1, somx2
      save rfac, iuse
c
      data iuse / 0 /
c
      plm( 1 ) = 1
      dplm( 1 ) = 0
c finished if monopoles only
      if( s3lmax .gt. 0 )then
c pretabulate constant factors
        if( iuse .ne. s3lmax )then
c check space
          j = 3 * ( ( s3lmax + 1 ) * s3lmax ) / 2
          if( j .gt. mfac )call crash( 'S3TPLM',
     +                                         'Local array too small' )
          i = 0
          do m = 0, s3lmax
            if( m .gt. 0 )then
              i = i + 1
              rfac( i ) = fact / sqrt( dble( 2 * m * ( 2 * m - 1 ) ) )
              fact = fact + 2
            else
              fact = 1
            end if
            if( m .lt. s3lmax )then
              i = i + 1
              rfac( i ) = sqrt( dble( 2 * m + 1 ) )
              if( m + 1 .lt. s3lmax )then
                do l = m + 2, s3lmax
                  i = i + 1
                  f2 = dble( l - 1 - m ) / dble( l - 1 + m )
                  rfac( i ) = sqrt( f2 ) * ( l + m - 1 )
                  i = i + 1
                  f2 = ( l - m ) * ( l + m )
                  rfac( i ) = 1. / sqrt( f2 )
                end do
              end if
            end if
          end do
c derivs
          do l = 1, s3lmax
            do m = 0, l
              if( m .lt. l )then
                i = i + 1
                f2 = dble( l - m ) / dble( l + m )
                rfac( i ) = sqrt( f2 ) * ( l + m )
              end if
            end do
          end do
          iuse = s3lmax
          if( i .ne. j )call crash( 'S3TPLM', 'Algorithm failed' )
        end if
c compute P_m^m
        i = 0
        do m = 0, s3lmax
          if( m .gt. 0 )then
            i = i + 1
            pmm = -pmm * somx2 * rfac( i )
            j = k + m + 1
            plm( j ) = pmm
            k = j
          else
            pmm = 1
            somx2 = sqrt( ( 1 - x ) * ( 1 + x ) )
            j = 1
            k = 1
          end if
          if( m .lt. s3lmax )then
            pm = pmm
c compute P_{m+1}^m
            i = i + 1
            pmp1 = x * rfac( i ) * pm
            j = j + m + 1
            plm( j ) = pmp1
            if( m + 1 .lt. s3lmax )then
c compute P_l^m when l > m + 1
              do l = m + 2, s3lmax
                pll = ( x * ( 2 * l - 1 ) * pmp1 -
     +              rfac( i + 1 ) * pm ) * rfac( i + 2 )
                i = i + 2
                pm = pmp1
                pmp1 = pll
                j = j + l
                plm( j ) = pll
              end do
            end if
          end if
        end do
        if( derivs )then
c compute and tabulate derivatives - formula (8.5.4) from Abramowitz & Stegun
          fact = x**2 - 1
c OK to return zeros if x=1, since derivatives don't contribute to forces
          if( fact .ne. 0. )fact = 1 / fact
c P00 term already set
          k = 1
          do l = 1, s3lmax
            do m = 0, l
              k = k + 1
              if( m .lt. l )then
                pm = plm( k - l )
                i = i + 1
                f2 = rfac( i )
              else
                pm = 0
                f2 = 0
              end if
              dplm( k ) = ( l * x * plm( k ) - f2 * pm ) * fact
            end do
          end do
        end if
      end if
      return
      end
