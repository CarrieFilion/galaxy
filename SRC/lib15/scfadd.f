      subroutine scfadd( jst )
      use aarrays
      implicit none
c routine to increment the coefficients of the basis expansion by adding the
c   contribution of the current group of particles.  Code adapted from that
c   originally written by G. Quinlan and L. Hernquist
c Modified by Jerry Sellwood and included in this package with permission
c
c Called by: MASSGN
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local arrays
      real*8 cmp( 0:s3maxl ), smp( 0:s3maxl )
c
c local variables
      integer is, j, jz, k, l, m, n
      real*8 costh, cp, r, rc, sp, t1, t2, t3, x, xi, y, z
c
      if( .not. ( ( basset .eq. 'hern' ) .or. ( basset .eq. 'plum' ) ) )
     +                      call crash( 'SCFADD', 'Unrecognized basis' )
      do is = 1, jst
        jz = max( iz( is ), 1 )
        x = newc( 1, is ) - xcpred( 1, jz, jgrid )
        y = newc( 2, is ) - xcpred( 2, jz, jgrid )
        z = newc( 3, is ) - xcpred( 3, jz, jgrid )
        r = sqrt( x * x + y * y + z * z )
c ultraspherical harmonics
        if( basset .eq. 'hern' )then
          xi = ( r - 1. ) / ( r + 1. )
        else
          xi = ( r**2 - 1. ) / ( r**2 + 1. )
        end if
        do l = 0, s3lmax
          ultrasp( 0, l ) = 1.0
          ultrasp( 1, l ) = twoalpha( l ) * xi
          t1 = ultrasp( 1, l )
          t2 = 1.0
          do n = 1, nbas - 1
            ultrasp( n + 1, l ) = ( c1( n, l ) * xi * t1 -
     +                              c2( n, l ) * t2 ) * c3( n )
            t2 = t1
            t1 = ultrasp( n + 1, l )
          end do
          do n = 0, nbas
            ultraspt( n, l ) = ultrasp( n, l ) * anltilde( n, l )
          end do
        end do
c make a table of cos(m\phi) and sin(m\phi)
        rc = sqrt( x * x + y * y + 1.d-20 )
        cp = x / rc
        sp = y / rc
        cmp( 0 ) = 1
        smp( 0 ) = 0
        do m = 1, s3lmax
          cmp( m ) = cmp( m - 1 ) * cp - smp( m - 1 ) * sp
          smp( m ) = smp( m - 1 ) * cp + cmp( m - 1 ) * sp
        end do
c get Plms
        costh = z / r
        call s3tplm( costh, .false. )
c
        k = 0
        j = 0
        do l = 0, s3lmax
c apply filter
          if( lg( l + 1 ) )then
            j = j + 2 * ( nbas + 1 ) * ( l + 1 )
            k = k + l + 1
          else
            t1 = r**l * pwt( is )
            if( basset .eq. 'hern' )then
              t1 = t1 / ( ( 1. + r )**( 2 * l + 1 ) )
            else
              t1 = t1 / ( ( 1. + r**2 )**l * sqrt( 1. + r**2 ) )
            end if
            do m = 0, l
              k = k + 1
              t2 = t1 * plm( k ) * coeflm( l, m ) * cmp( m )
              t3 = t1 * plm( k ) * coeflm( l, m ) * smp( m )
              do n = 0, nbas
                sfpmss( j + 1, jz ) =
     +                       sfpmss( j + 1, jz ) + t2 * ultraspt( n, l )
                sfpmss( j + 2, jz ) =
     +                       sfpmss( j + 2, jz ) + t3 * ultraspt( n, l )
                j = j + 2
              end do
            end do
          end if
        end do
      end do
      ncontrib = ncontrib + jst
      return
      end
