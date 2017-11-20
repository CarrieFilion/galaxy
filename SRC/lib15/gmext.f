      real*8 function gmext( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c for adiabatic compression: returns externally imposed mass interior to r
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      integer nt, iusegm
      parameter ( nt = 1001 )
      real*8, allocatable :: radtab( : ), gmtab( : ), lam( : ), c( : )
      save c, gmtab, iusegm, lam, radtab
c
c externals
      external gdsphf
      real zthick
      real*8 gmassh, gmtabd, quad_PAT, splint2
c
c local variables
      integer i, ifail, n
      real*8 epsr, r2
      data iusegm / 0 /
c
      gmext = -1
      if( disc( icmp ) )then
        r2 = r / rscale( icmp )
c finitely thick disks
        if( zthick( r ) .gt. 0. )then
          if( iusegm .eq. 0 )then
            if( master )print *, 'GMEXT: Building a table'
            allocate ( radtab( nt ) )
            allocate ( gmtab( nt ) )
            radtab( 1 ) = 0
            gmtab( 1 ) = 0
            epsr = 1.d-6
            do i = 2, nt
              r2 = rtrunc( icmp ) * real( i - 1 ) / real( nt - 1 )
              radtab( i ) = r2
              n = 0
              ifail = 0
              gmtab( i ) = quad_Pat( gdsphf, 0.d0, r2, epsr, ifail )
            end do
            if( master )print *, 'GMEXT: Table ready'
            iusegm = 1
c set up spline
            allocate ( lam( nt + 4) )
            allocate ( c( nt + 4 ) )
            gmext = splint2(
     +                radtab, gmtab, nt, .5 * r2, lam, c, .true. )
            r2 = r / rscale( icmp )
          end if
c interpolate
          gmext = 0
          if( ( r2 .gt. 0.d0 ) .and. ( r2 .lt. radtab( nt ) ) )then
            gmext =
     +           splint2( radtab, gmtab, nt, r2, lam, c, .false. )
            gmext = max( gmext, 0.d0 )
          else if( r2 .ge. radtab( nt ) )then
            gmext = gmtab( nt )
          end if
c razor thin disks
        else
c Kuz'min/Toomre sphere
          if( ctype( icmp ) .eq. 'KT  ' )then
            r2 = r2 * r2
            if( r2 .lt. 1.d-5 )then
              gmext = r2 * ( .5 - 7. * r2 / 24. )
            else
              gmext = 1. - 1. / sqrt( 1. + r2 )
            end if
            gmext = gmext * cmpmas( icmp )
          else if( ctype( icmp ) .eq. 'EXP ' )then
c exponential sphere
            if( r2 .gt. 1.d-4 )then
              gmext = 1. - ( 1. + r2 ) * exp( -r2 )
            else
c expansion for small arguments
              gmext = r2**2 * ( .5 - r2 / 3. + r2**2 / 8. )
            end if
            gmext = gmext * cmpmas( icmp )
          else if( ctype( icmp ) .eq. 'MTAB' )then
c tabulated M(R)
            gmext = gmtabd( r2 )
          end if
        end if
c not a disk
      else
        gmext = gmassh( r )
      end if
      if( gmext .lt. 0. )call crash( 'GMEXT', 'Unrecognized component' )
      return
      end
