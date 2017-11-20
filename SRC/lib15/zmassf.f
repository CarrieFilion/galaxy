      real*8 function zmassf( z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the mass fraction integrated from -infinity to z / z0
c   uses QUADPACK routines DQAGS and DQAGI
c
c calling argument
      real*8 z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      external rhozs
      real*8 algrng2, errfn, quad_gnr, quad_inf
c
c local arrays
      integer m
      parameter ( m = 101 )
      real*8, allocatable :: mt( : ), zt( : )
      save mt, zt
c
c local variables
      integer ier, iuse, iz
      real*8 epsa, epsr, ex, res
      save iuse
c
      data iuse / 0 /
c
      if( ( iztyp( icmp ) .eq. 1 ) .or. ( iztyp( icmp ) .eq. 3 ) )then
c function for Spitzer sheet is ( 1 + tanh( z ) ) / 2
        ex = exp( .5 * z )
        ex = ( ex - 1. / ex ) / ( ex + 1. / ex )
        zmassf = .5 * ( 1. + ex )
      else if( iztyp( icmp ) .eq. 2 )then
c Gaussian rule - error function, ier is ignored
        ex = z / sqrt( 2. )
        zmassf = .5 + .5 * errfn( ex, ier )
      else if( iztyp( icmp ) .eq. 5 )then
c exponential function
        zmassf = 1. - exp( -abs( z ) )
        zmassf = .5 * ( 1. + sign( zmassf, z ) )
      else
c integral of some arbitrary density function
        if( iuse .ne. icmp )then
c (re-)compute normalization
          znorm( icmp ) = 1
          epsr = 1.d-8
          epsa = epsr
          res = quad_inf( rhozs, 0.d0, 1, epsa, epsr, ier )
          if( ier .gt. 0 )then
            print *, 'ier =', ier, ' from QUAD_INF'
            call crash( 'ZMASSF', 'QUADPACK error' )
          end if
          znorm( icmp ) = 2. * res
          allocate ( mt( m ) )
          allocate ( zt( m ) )
c compute table
          do iz = 1, m
            zt( iz ) = .1 * real( iz - 1 )
c compute mass fraction to input height
            epsr = 1.d-8
            epsa = epsr
            ier = 0
            res = quad_gnr( rhozs, 0.d0, zt( iz ), epsa, epsr, ier )
            if( ier .gt. 0 )then
              print *, 'ier =', ier, ' from QUAD_GNR'
              call crash( 'ZMASSF', 'QUADPACK error' )
            end if
            mt( iz ) = res
          end do
          iuse = icmp
        end if
        if( abs( z ) .lt. zt( m ) )then
          zmassf = algrng2( zt, mt, m, abs( z ) )
          if( z .gt. 0. )then
            zmassf = zmassf + .5
          else
            zmassf = .5 - zmassf
          end if
        else
          zmassf = 0
          if( z .gt. 0.d0 )zmassf = 1
        end if
      end if
      return
      end
