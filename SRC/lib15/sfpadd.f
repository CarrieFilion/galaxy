      subroutine sfpadd( jst )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to increment the coefficients of the basis expansion by adding the
c   contribution of the current group of particles
c The cylindrical radius, plane number and weights were pre-calculated in
c   a call to WEIGHT
      use aarrays
      implicit none
c
c Called by: MASSGN or SFPSUM
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
      real*8, allocatable :: cmph(:,:), smph(:,:)
c
c local variables
      integer ifun, is, j, jz, k, m, n
      real*8 c, d, f, s, x, y
c
      allocate ( cmph( mostm, mbuff ) )
      allocate ( smph( mostm, mbuff ) )
c evaluate needed cosines and sines
      if( maxm .gt. 0 )then
        do is = 1, jst
          if( nskip( is ) )then
            j = max( iz( is ), 1 )
            x = newc( 1, is ) - xcpred( 1, j, jgrid )
            y = newc( 2, is ) - xcpred( 2, j, jgrid )
c use cos = x/r & sin = y/r
            cmph( 1, is ) = x / rr( is )
            smph( 1, is ) = y / rr( is )
c recurrence relations
            if( maxm .gt. 1 )then
              do m = 2, maxm
                cmph( m, is ) = cmph( m - 1, is ) * cmph( 1, is ) -
     +                          smph( m - 1, is ) * smph( 1, is )
                smph( m, is ) = cmph( m - 1, is ) * smph( 1, is ) +
     +                          smph( m - 1, is ) * cmph( 1, is )
              end do
            end if
c end skip over particles outside active region
          end if
c end loop over particles
        end do
      end if
c
c 'standard' 2-D basis sets
c
      if( ( basset .eq. 'ablj' ) .or.
     +    ( ( basset .eq. 'bess' ) .and. sf2d ) )then
        do is = 1, jst
          if( nskip( is ) )then
c count particles contributing to coefficients
            ncontrib = ncontrib + 1
c look up all function values for this radius
            call sfpval( dble( rr( is ) ) )
            jz = max( iz( is ), 1 )
c loop through all selected functions
            do ifun = 1, lastf
              m = msel( ifun )
              n = nsel( ifun )
              f = funvals( ifun )
c update the coefficients
              k = n * ngx + maxm + 1
              if( m .eq. 0 )then
                sfpmss( k, jz ) = sfpmss( k, jz ) + f
              else
                c = cmph( m, is )
                s = smph( m, is )
c positive m refers to cos terms and negative m to sin terms
                sfpmss( k + m, jz ) = sfpmss( k + m, jz ) + f * c
                sfpmss( k - m, jz ) = sfpmss( k - m, jz ) + f * s
              end if
            end do
c end skip over particles outside active region
          end if
c end loop over particles
        end do
c
c the only 3-D basis
c
      else if( ( basset .eq. 'bess' ) .and. sf3d )then
        do is = 1, jst
          if( nskip( is ) )then
c count particles contributing to coefficients
            ncontrib = ncontrib + 1
c look up all function values for this radius
            call sfpval( dble( rr( is ) ) )
            jz = max( iz( is ), 1 )
c plane number
            j = ncl( is )
c loop through all selected functions
            do ifun = 1, lastf
              m = msel( ifun )
              n = nsel( ifun )
              f = funvals( ifun )
              k = ( j - 1 ) * ngxy + n * ngx + maxm + 1
              if( m .eq. 0 )then
                sfpmss( k, jz ) = sfpmss( k, jz ) + wt( 1, is ) * f
                k = k + ngxy
                sfpmss( k, jz ) = sfpmss( k, jz ) + wt( 2, is ) * f
              else
                c = cmph( m, is )
                s = smph( m, is )
c positive m refers to cos terms and negative m to sin terms
                sfpmss( k + m, jz ) =
     +                         sfpmss( k + m, jz ) + wt( 1, is ) * f * c
                sfpmss( k - m, jz ) =
     +                         sfpmss( k - m, jz ) + wt( 1, is ) * f * s
                k = k + ngxy
                sfpmss( k + m, jz ) =
     +                         sfpmss( k + m, jz ) + wt( 2, is ) * f * c
                sfpmss( k - m, jz ) =
     +                         sfpmss( k - m, jz ) + wt( 2, is ) * f * s
              end if
            end do
c end skip over particles outside active region
          end if
c end loop over particles
        end do
c
c a complex 2-D basis
c
      else if( basset .eq. 'lgsp' )then
        do is = 1, jst
          if( nskip( is ) )then
c count particles contributing to coefficients
            ncontrib = ncontrib + 1
c look up all function values for this radius
            call sfpvld( dble( rr( is ) ) )
            jz = max( iz( is ), 1 )
c loop through all selected functions
            do ifun = 1, lastf
              m = msel( ifun )
              n = nsel( ifun )
              f = funvals( ifun )
              d = dervals( ifun )
c update the coefficients
              k = n * ngx + maxm + 1
              if( m .eq. 0 )then
                call crash( 'SFPADD', 'Option not programmed' )
                sfpmss( k, jz ) = sfpmss( k, jz ) + f
              else
                c = cmph( m, is )
                s = smph( m, is )
c positive m refers to real parts and negative m to imaginary parts
                sfpmss( k + m, jz ) =
     +                               sfpmss( k + m, jz ) + c * f - s * d
                sfpmss( k - m, jz ) =
     +                               sfpmss( k - m, jz ) - s * f - c * d
              end if
            end do
c end skip over particles outside active region
          end if
c end loop over particles
        end do
      else
        call crash( 'SFPADD', 'Unrecognized basis' )
      end if
      deallocate ( cmph )
      deallocate ( smph )
      return
      end
