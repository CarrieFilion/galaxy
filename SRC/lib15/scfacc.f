      subroutine scfacc( jst )
      use aarrays
      implicit none
c Computes the acceleration components for the current group of particles
c    using the basis function expansion coeffs only.  Code adapted from
c    that originally written by G. Quinlan and L. Hernquist
c Modified by Jerry Sellwood and included in this package with permission
c
c This routine is the SCF equivalent of GRIDF and is called by GETACC
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
      integer is, j, k, l, m, n
      real*8 aphi, ar, ath, clm, cp, dlm, elm, flm, phinltil, potdum, sp
      real*8 r, rc, t1, t2, t3, t4, x, xi, y, z
c
      if( .not. ( ( basset .eq. 'hern' ) .or. ( basset .eq. 'plum' ) ) )
     +                      call crash( 'SCFACC', 'Unrecognized basis' )
      do is = 1, jst
        x = oldc( 1, is ) - xcen( 1, 1 )
        y = oldc( 2, is ) - xcen( 2, 1 )
        z = oldc( 3, is ) - xcen( 3, 1 )
        r = sqrt( x * x + y * y + z * z )
        rr( is ) = r
        nskip( is ) = .true.
c ultraspherical harmonics
        if( basset .eq. 'hern' )then
          xi = ( r - 1. ) / ( r + 1. )
        else
          xi = ( r**2 - 1. ) / ( r**2 + 1. )
        end if
        do l = 0, s3lmax
          ultrasp( 0, l ) = 1.0
          ultrasp( 1, l ) = twoalpha( l ) * xi
          ultrasp1( 0, l ) = 0.0
          ultrasp1( 1, l ) = 1.0
          t1 = ultrasp( 1, l )
          t2 = 1.0
          do n = 1, nbas - 1
            ultrasp( n + 1, l ) = ( c1( n, l ) * xi * t1 -
     +                              c2( n, l ) * t2 ) * c3( n )
            t2 = t1
            t1 = ultrasp( n + 1, l )
            ultrasp1( n + 1, l ) = ( ( twoalpha( l ) +
     +                 ( n + 1 ) - 1. ) * t2 - ( n + 1 ) * xi *
     +      ultrasp( n + 1, l ) ) / ( twoalpha( l ) * ( 1. - xi * xi ) )
          end do
        end do
c make a table of cos(m\phi) and sin(m\phi)
        rc = sqrt( x**2 + y**2 + 1.d-20 )
        cp = x / rc
        sp = y / rc
        cmp( 0 ) = 1
        smp( 0 ) = 0
        do m = 1, s3lmax
          cmp( m ) = cmp( m - 1 ) * cp - smp( m - 1 ) * sp
          smp( m ) = smp( m - 1 ) * cp + cmp( m - 1 ) * sp
        end do
c get Plms and derivatives
        t1 = z / r
        call s3tplm( t1, .true. )
c
        potdum = 0.0
        ar = 0.0
        ath = 0.0
        aphi = 0.0
        k = 0
        j = 0
        do l = 0, s3lmax
c apply filter
          if( lg( l + 1 ) )then
            j = j + 2 * ( l + 1 ) * ( nbas + 1 )
            k = k + l + 1
          else
            t1 = 0.0
            t2 = 0.0
            t3 = 0.0
            t4 = 0.0
            do m = 0, l
              clm = 0.0
              dlm = 0.0
              elm = 0.0
              flm = 0.0
              do n = 0, nbas
                clm = clm + ultrasp( n, l ) * sfpfld( j + 1, 1 )
                dlm = dlm + ultrasp( n, l ) * sfpfld( j + 2, 1 )
                elm = elm + ultrasp1( n, l ) * sfpfld( j + 1, 1 )
                flm = flm + ultrasp1( n, l ) * sfpfld( j + 2, 1 )
                j = j + 2
              end do
              k = k + 1
              t1 = t1 + plm( k ) * ( clm * cmp( m ) + dlm * smp( m ) )
              t2 = t2 - plm( k ) * ( elm * cmp( m ) + flm * smp( m ) )
              t3 = t3 - dplm( k) * ( clm * cmp( m ) + dlm * smp( m ) )
              t4 = t4 - m * plm( k ) * ( dlm * cmp( m ) -
     +                                   clm * smp( m ) )
            end do
            if( basset .eq. 'hern' )then
              phinltil = r**l / ( ( 1. + r )**( 2 * l + 1 ) )
              ar = ar + phinltil * ( -t1 * ( l / r - ( 2 * l + 1 ) /
     +       ( 1. + r ) ) + t2 * 4. * ( 2. * l + 1.5 ) / ( 1. + r )**2 )
            else
              phinltil = r**l / ( ( 1. + r**2 )**l * sqrt( 1. + r**2 ) )
              ar = ar + phinltil * ( -t1 * ( l / r - ( 2 * l + 1 ) * r /
     +    ( 1. + r**2 ) ) + t2 * 8. * r * ( l + 1 ) / ( 1. + r**2 )**2 )
            end if
            ath = ath + t3 * phinltil
            aphi = aphi + t4 * phinltil
            potdum = potdum + t1 * phinltil
          end if
        end do
        ath = -rc * ath / r**2
        aphi = aphi / rc
c resolve into Cartesians and save
        acc( 1, is ) = acc( 1, is ) +
     +                     ( ar + z * ath / rc ) * x / r - y * aphi / rc
        acc( 2, is ) = acc( 2, is ) +
     +                     ( ar + z * ath / rc ) * y / r + x * aphi / rc
        acc( 3, is ) = acc( 3, is ) + ( z * ar - rc * ath ) / r
        gpot( is ) = gpot( is ) + potdum
      end do
      return
      end
