      subroutine sfpacc( jst )
c  Copyright (C) 2015, Jerry Sellwood
c
c Computes the acceleration components for the current group of particles
c    using the basis function expansion coeffs only - i.e., compute
c
c                                    dPhi    1 dPhi    dPhi
c     a  =  ( a , a   , a )   =  ( - ----, - - ----, - ---- ) .
c              R   phi   z            dR     R dphi     dz
c
c This routine is the SFP equivalent of GRIDF and is called by GETACC
c    The cylindrical radius, plane number and weights are pre-calculated
c    by a call to WEIGHT
c The SFP force is defined in cylindrical polar components and has to be
c    converted to Cartesians when needed.
      use aarrays
      implicit none
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
      real cmph( mostm, mbuff ), smph( mostm, mbuff )
c
c local variables
      integer ifun, is, j, l, m, n
      real al, ar, at, az, c, ci, cr, d, e, f, k, r, s, x, y
c
c compute plane numbers and weights for interpolation
      call weight( jst, .true. )
c evaluate required cosines and sines
      if( maxm .gt. 0 )then
c work through all particles
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
            x = oldc( 1, is ) - xcen( 1, 1 )
            y = oldc( 2, is ) - xcen( 2, 1 )
            r = max( sqrt( x * x + y * y ), 1.e-5 )
c use cos = x/r & sin = y/r
            cmph( 1, is ) = x / r
            smph( 1, is ) = y / r
c recurrence relations
            if( maxm .gt. 1 )then
              do m = 2, maxm
                cmph( m, is ) = cmph( m - 1, is ) * cmph( 1, is ) -
     +                          smph( m - 1, is ) * smph( 1, is )
                smph( m, is ) = cmph( m - 1, is ) * smph( 1, is ) +
     +                          smph( m - 1, is ) * cmph( 1, is )
              end do
            end if
          end if
        end do
      end if
c
c 'standard' 2-D basis sets
c
      if( ( basset .eq. 'ablj' ) .or.
     +    ( ( basset .eq. 'bess' ) .and. sf2d ) )then
c work through all particles
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
c look up functions and derivatives for this r
            call sfpvld( dble( rr( is ) ) )
            rr( is ) = max( rr( is ), 1.e-5 )
c sum over functions
            do ifun = 1, lastf
              m = msel( ifun )
              n = nsel( ifun )
              f = funvals( ifun )
              d = dervals( ifun )
              l = n * ngx + maxm + 1
              cr = sfpfld( l + m, 1 )
              if( m .eq. 0 )then
                acc( 1, is ) = acc( 1, is ) + d * cr
                if( phys )gpot( is ) = gpot( is ) - cr * f
              else
                c = cmph( m, is )
                s = smph( m, is )
                ci = sfpfld( l - m, 1 )
                acc( 1, is ) = acc( 1, is ) + d * ( ci * s + cr * c )
                acc( 2, is ) = acc( 2, is ) +
     +                               real( m ) * f * ( ci * c - cr * s )
                if( phys )gpot( is ) = gpot( is ) -
     +                                           f * ( ci * s + cr * c )
              end if
            end do
c radius factor in tangential force - rr( is ) is the cylindrical radius
            if( basset .eq. 'ablj' )then
              acc( 2, is ) = acc( 2, is ) * maxr / rr( is )
            else
              acc( 2, is ) = acc( 2, is ) / rr( is )
            end if
c end skip over off-grid particles
          end if
        end do
c
c only 3-D basis set
c
      else if( ( basset .eq. 'bess' ) .and. sf3d )then
c work through all particles
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
c look up functions and derivatives for this r
            call sfpvld( dble( rr( is ) ) )
            rr( is ) = max( rr( is ), 1.e-5 )
c particle within range of grid planes
            if( ( ncl( is ) .gt. 0 ) .and.
     +          ( ncl( is ) .lt. ngz ) )then
              j = ncl( is )
c sum over functions
              do ifun = 1, lastf
                m = msel( ifun )
                n = nsel( ifun )
                f = funvals( ifun )
                d = dervals( ifun )
c interpolate between planes
                l = ( j - 1 ) * ngxy + n * ngx + maxm + 1
                if( m .eq. 0 )then
                  ar = sfpfld( l, 1 ) * wt( 1, is ) +
     +                 sfpfld( l + ngxy, 1 ) * wt( 2, is )
                  az = sfpfld( l, 2 ) * wt( 1, is ) +
     +                 sfpfld( l + ngxy, 2 ) * wt( 2, is )
                  acc( 1, is ) = acc( 1, is ) + d * ar
                  acc( 3, is ) = acc( 3, is ) + f * az
                  if( phys )gpot( is ) = gpot( is ) - f * ar
                else
                  c = cmph( m, is )
                  s = smph( m, is )
                  ar = ( sfpfld( l - m, 1 ) * s +
     +                   sfpfld( l + m, 1 ) * c ) * wt( 1, is ) +
     +                 ( sfpfld( l - m + ngxy, 1 ) * s +
     +                   sfpfld( l + m + ngxy, 1 ) * c ) * wt( 2, is )
                  at = ( sfpfld( l - m, 1 ) * c -
     +                   sfpfld( l + m, 1 ) * s ) * wt( 1, is ) +
     +                 ( sfpfld( l - m + ngxy, 1 ) * c -
     +                   sfpfld( l + m + ngxy, 1 ) * s ) * wt( 2, is )
                  az = ( sfpfld( l - m, 2 ) * s +
     +                   sfpfld( l + m, 2 ) * c ) * wt( 1, is ) +
     +                 ( sfpfld( l - m + ngxy, 2 ) * s +
     +                   sfpfld( l - m + ngxy, 2 ) * c ) * wt( 2, is )
                  acc( 1, is ) = acc( 1, is ) + d * ar
                  acc( 2, is ) = acc( 2, is ) + real( m ) * f * at
                  acc( 3, is ) = acc( 3, is ) + f * az
                  if( phys )gpot( is ) = gpot( is ) - f * ar
                end if
              end do
            else
c particle above or below grid planes
              j = max( ncl( is ), 1 )
              j = min( j, ngz )
c sum over functions
              do ifun = 1, lastf
                m = msel( ifun )
                n = nsel( ifun )
                k = deltak * real( n )
                e = exp( -k * wt( 1, is ) )
                f = funvals( ifun )
                d = dervals( ifun )
c exponential dilution factor from nearest plane
                l = ( j - 1 ) * ngxy + n * ngx + maxm + 1
                if( m .eq. 0 )then
                  ar = sfpfld( l, 1 ) * e
                  az = ar * k * wt( 2, is )
                  acc( 1, is ) = acc( 1, is ) + d * ar
                  acc( 3, is ) = acc( 3, is ) + f * az
                  if( phys )gpot( is ) = gpot( is ) - f * ar
                else
                  c = cmph( m, is )
                  s = smph( m, is )
                  ar =
     +           ( sfpfld( l - m, 1 ) * s + sfpfld( l + m, 1 ) * c ) * e
                  at =
     +           ( sfpfld( l - m, 1 ) * c - sfpfld( l + m, 1 ) * s ) * e
                  az = ar * k * wt( 2, is )
                  acc( 1, is ) = acc( 1, is ) + d * ar
                  acc( 2, is ) = acc( 2, is ) + real( m ) * f * at
                  acc( 3, is ) = acc( 3, is ) + f * az
                  if( phys )gpot( is ) = gpot( is ) - f * ar
                end if
              end do
            end if
c radius factor in tangential force - rr( is ) is the cylindrical radius
            acc( 2, is ) = acc( 2, is ) / rr( is )
c end skip over off-grid particles
          end if
        end do
c
c complex 2-D basis sets
c
      else if( basset .eq. 'lgsp' )then
c work through all particles
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
c look up functions and derivatives for this r
            call sfpvld( dble( rr( is ) ) )
            rr( is ) = max( rr( is ), 1.e-5 )
c sum over pitch angles
            do ifun = 1, lastf
              m = msel( ifun )
              n = nsel( ifun )
              al = real( n - maxn / 2 ) * deltak
              f = funvals( ifun )
              d = dervals( ifun )
              l = n * ngx + maxm + 1
              cr = sfpfld( l + m, 1 )
              ci = sfpfld( l - m, 1 )
              if( m .eq. 0 )then
                call crash( 'SFPACC', 'Option not programmed' )
                acc( 1, is ) = acc( 1, is ) +
     +               cr * ( .5 * f + al * d ) - ci * ( .5 * d - al * f )
                if( phys )gpot( is ) = gpot( is ) + cr * f - ci * d
              else
                c = cmph( m, is )
                s = smph( m, is )
                acc( 1, is ) = acc( 1, is ) -
     +                       ( .5 * cr + al * ci ) * ( c * f - s * d ) +
     +                       ( .5 * ci - al * cr ) * ( s * f + c * d )
                acc( 2, is ) = acc( 2, is ) - real( m ) *
     +               ( cr * ( s * f + c * d ) + ci * ( c * f - s * d ) )
                if( phys )gpot( is ) = gpot( is ) -
     +                   cr * ( c * f - s * d ) + ci * ( s * f + c * d )
              end if
            end do
c radius factor
            acc( 1, is ) = acc( 1, is ) / rr( is )
            acc( 2, is ) = acc( 2, is ) / rr( is )
c end skip over off-grid particles
          end if
        end do
      else
        call crash( 'SFPACC', 'Unrecognized basis' )
      end if
c resolve forces into Cartesian components
      do is = 1, jst
c skip particles off grid
        if( nskip( is ) )then
          ar = acc( 1, is )
          at = acc( 2, is )
          acc( 1, is ) = ar * cmph( 1, is ) - at * smph( 1, is )
          acc( 2, is ) = ar * smph( 1, is ) + at * cmph( 1, is )
c end skip over off-grid particles
        end if
      end do
      return
      end
