      real*8 function Phispl( r, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Returns a value from a spline fit to an azimuthally or spherically averaged
c   potential.  The calling arguments and returned value are in natural units.
c   A dummy potential is returned when analyt is set to .true. as needed
c   for the start of an iterative solutions.
c   uses NAG routines E02BBF and E02DEF
c
c calling arguments
      real*8 r, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c allocatable arrays
      real*8, allocatable, save :: Bc1(:), Bc2(:), tr(:), tz(:), w(:)
      real*8, allocatable, save :: pt(:), rt(:)
c
c externals
      real*8 db2val, isopsi, splint2
c
c local arrays
      integer jdftyp( mcmp ), jmpar( mcmp ), jmtyp( mcmp )
      real rtc( mcmp ), sdfcn( 3, mcmp )
c
c local variables
      integer icode, ifail, ip, ir, iu, iuse, j, jcmp, k, lcmp, lr, lz
      integer nrspl, nzspl
      real a, stype
      real*8 r2, z2
      save iuse, lr, lz, nrspl, nzspl
c
      data iuse / 0 /
c
c arbitrary analytic potential needed for start-up of DFITER
      if( analyt )then
        r2 = sqrt( r * r + z * z )
c Plummer form
        if( ( impar( icmp ) .eq. 0 ) .or. ( impar( icmp ) .eq. 1 ) )then
          Phispl = -cmpmas( icmp ) / sqrt( 1. + r2 * r2 )
c King model
        else if( impar( icmp ) .eq. 2 )then
          Phispl = -cmpmas( icmp ) / rtidal - isopsi( r2 )
c Jaffe model form
        else if( impar( icmp ) .eq. 3 )then
          r2 = max( r2, 1.d-12 )
          Phispl = cmpmas( icmp ) * log( r2 / ( 1. + r2 ) )
c Hernquist model
        else if( impar( icmp ) .eq. 4 )then
          Phispl = -cmpmas( icmp ) / ( 1. + r2 )
        else
          call crash( 'PHISPL', 'Unknown analytic potential' )
        end if
      else
c read spline fit if required
        if( initls )then
          iu = -1
          call opnfil( iu, 'pft', 'unformatted', 'old', 'seq', ifail )
          if( ifail .ne. 0 )call crash( 'PHISPL',
     +                                       'old .pft file not found' )
          read( iu )nrspl, nzspl, lcmp, jcmp, icode
c read and check headers
          do ip = 1, lcmp
            read( iu )jmtyp( ip ), jmpar( ip ), jdftyp( ip ),
     +                ( sdfcn( j, ip ), j = 1, 3 ), cmpmas( ip ),
     +                fmass( ip ), rscale( ip ), rtc( ip ),
     +                disc( ip ), stype
            write( ctype( ip ), '( a4 )' )stype
          end do
c check header values
          k = abs( ncmp - lcmp ) + abs( ncode - icode )
          a = 0
          do ip = 1, lcmp
            k = k + abs( jmtyp( ip ) - imtyp( ip ) )
            k = k + abs( jmpar( ip ) - impar( ip ) )
            k = k + abs( jdftyp( ip ) - idftyp( ip ) )
            a = a + abs( sdfcn( 1, ip ) - dfcns( 1, ip ) )
            do j = 2, 3
              dfcns( j, ip ) = sdfcn( j, ip )
            end do
          end do
          if( ( k .gt. 0 ) .or. ( a .gt. 0. ) )then
            if( master )then
              write( no, * )
     +               'Number components expected', ncmp, ' found', lcmp
              write( no, * )'ncode expected', ncode, ' found', icode
              do ip = 1, lcmp
                write( no, * )'For component', ip
                write( no, '( ''  Expected values'', 3i4, 3f10.5 )' )
     +                           imtyp( ip ), impar( ip ), idftyp( ip ),
     +                                      ( dfcns( j, ip ), j = 1, 3 )
                write( no, '( '' .pft file values'', 3i4, 3f10.5 )' )
     +                           jmtyp( ip ), jmpar( ip ), jdftyp( ip ),
     +                                      ( sdfcn( j, ip ), j = 1, 3 )
              end do
            end if
            call crash( 'PHISPL', 'Wrong .pft file' )
          end if
c check truncation radii
          do ip = 1, ncmp
            if( abs( rtc( ip ) - rtrunc( ip ) ) .gt. 1.e-4 )then
              if( master )then
                print 201, ip, rtrunc( ip ), rtc( ip )
                write( no, 201 )ip, rtrunc( ip ), rtc( ip )
              end if
  201         format( 'Warning: rtrunc for compnt', i2, ' was', f7.2,
     +                    '   Value reset to', f6.2, ' from .pft file' )
              rtrunc( ip ) = rtc( ip )
            end if
          end do
          if( icmp .ne. jcmp )call crash( 'PHISPL', 'Wrong component' )
          lr = nrspl - 4
          lz = nzspl - 4
c read data
          if( sphrod( icmp ) )then
            if( iuse .eq. 0 )then
              allocate ( tr( nrspl ) )
              allocate ( tz( nzspl ) )
              allocate ( Bc2( lr * lz ) )
              allocate ( w( lz ) )
            end if
            read( iu )( tr( ir ), ir = 1, nrspl ),
     +                ( tz( ir ), ir = 1, nzspl ),
     +                ( Bc2( ir ), ir = 1, lr * lz )
            iuse = 1
          else
            allocate ( rt( lr ) )
            allocate ( pt( lr ) )
            allocate ( tr( nrspl ) )
            allocate ( bc1( nrspl ) )
            read( iu )( tr( ir ), ir = 1, nrspl ),
     +                ( Bc1( ir ), ir = 1, nrspl )
c assume abscissae are the same as the knots
            do ir = 1, lr
              rt( ir ) = tr( ir + 2 )
            end do
          end if
          close( iu )
          initls = .false.
        end if
        if( sphrod( icmp ) )then
c spheroidal potential - ensure point is within box
          r2 = min( abs( r ), 0.999d0 * tr( nrspl ) )
          z2 = min( abs( z ), 0.999d0 * tz( nzspl ) )
          if( .not. allocated( w ) )allocate( w( 4 * lr ) )
          Phispl = db2val( z2, r2, 0, 0, tz, tr, lz, lr, 4, 4, Bc2, w )
c approximation for values outside defined range
          if( ( abs( r ) .gt. r2 ) .or.
     +        ( abs( z ) .gt. z2 ) )Phispl = Phispl *
     +                       sqrt( r2**2 + z2**2 ) / sqrt( r**2 + z**2 )
        else
c spherical potential
          ifail = 0
          if( z .ne. 0. )call crash( 'PHISPL', 'Non-zero z' )
c ensure point is within range
          r2 = min( abs( r ), rt( lr ) )
          Phispl = splint2( rt, pt, lr, r2, tr, Bc1, .false. )
c approximation for values outside defined range
          if( abs( r ) .gt. rt( lr ) )Phispl =
     +                                      Phispl * rt( lr ) / abs( r )
        end if
      end if
      return
      end
