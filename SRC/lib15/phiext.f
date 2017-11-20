      real*8 function Phiext( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns potential due to a rigid mass component for compressed halos
c   may use QUADPACK routine DQAGI
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      external frext
      real zthick
      real*8 Phihal, splint2, quad_inf
c
c local arrays
      integer nt
      parameter ( nt = 201 )
      real*8, allocatable :: c( : ), lam( : ), pt( : ), rt( : )
c
c local variables
      integer i, ier, iuse
      real*8 epsa, epsr, p, r2, x
      include 'inc/pi.f'
      save iuse, c, lam, pt, rt
      data iuse / 0 /
c
      Phiext = 1
      if( .not. disc( icmp ) )then
c spherical masses are straightforward
        Phiext = Phihal( r )
      else
        x = r / rscale( icmp )
        if( zthick( icmp ) .gt. 0. )then
          if( iuse .eq. 0 )then
            allocate ( c( nt + 4 ) )
            allocate ( lam( nt + 4 ) )
            allocate ( pt( nt ) )
            allocate ( rt( nt ) )
            if( master )print *, 'Building a table in Phiext'
            do i = 1, nt
c choose abscissae in table
              r2 = real( i - 1 ) / real( nt - 1 )
              r2 = 50. * sin( .5 * pi * r2 )**2
c integrate radial force from enclosed mass
              epsa = 1.e-8
              epsr = epsa
              p = quad_inf( frext, r2, 1, epsa, epsr, ier )
c ier = 2 means requested accuracy was not achieved
              if( ier .ne. 0 .and. ier .ne. 2 )then
                print *, 'ier =', ier, ' from QUAD_INF'
                call crash( 'PHIEXT', 'QUADPACK error' )
              end if
              rt( i ) = r2
              pt( i ) = p
            end do
c initialize spline
            r2 = .5 * r2
            p = splint2( rt, pt, nt, r2, lam, c, .true. )
            iuse = 1
            if( master )print *, 'Table ready'
          end if
c interpolate in table
          if( r .lt. rt( nt ) )then
            Phiext = splint2( rt, pt, nt, r, lam, c, .false. )
          else
c assume all mass is enclosed at last measured point
            Phiext = pt( nt ) * rt( nt ) / r
          end if
        else
c this function needs the monopole term of the disk only - ie a spherical mass!
          if( ctype( icmp ) .eq. 'KT  ' )then
            if( x .gt. 1.d-4 )then
              Phiext = -( 1. - ( sqrt( 1. + x * x ) - 1. ) / x )
            else
              Phiext = -( 1. - .5 * x + .125 * x**3 - .0625 * x**5 )
            end if
          else if( ctype( icmp ) .eq. 'EXP ' )then
c exponential sphere
            if( r .gt. 1.d-4 )then
              Phiext = -( 1. - exp( -x ) ) / x
            else
c expansion for small arguments
              Phiext = -( 1. - .5 * x + x**2 / 6. - x**3 / 24. )
            end if
c tabulated M(r)
          else if( ctype( icmp ) .eq. 'MTAB' )then
            call crash( 'PHIEXT', 'zero thickness option unavailable' )
          else
            call crash( 'PHIEXT', 'Unrecognzed disk type' )
          end if
          Phiext = cmpmas( icmp ) * Phiext / rscale( icmp )
        end if
      end if
      if( Phiext .gt. 0.d0 )call crash( 'Phiext', 'Value not set' )
      return
      end
