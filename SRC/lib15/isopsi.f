      real*8 function isopsi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine for numerical evaluation of psi(r) for spherical isotropic models
c   as in chapter 4 of Binney & Tremaine
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
c externals
      external isopsf
      real*8 aitken2
c
c local array
      real*8 y( 2 )
c
c local variables
      integer i, ifail
      real*8 r0, r1, tol
c
c inside maximum radius
      if( r .lt. rtidal )then
c initialise table on first call
        if( iusedf .eq. 0 )then
          if( master )print *, 'ISOPSI: Building a table'
          do i = 1, nrdd
            rtab( i ) = rtidal * real( i - 1 ) / real( nrdd - 1 )
            rtab( i ) = min( rtab( i ), rtidal )
            r0 = 0
            tol = 1.d-8
            r1 = rtab( i )
            if( ctype( icmp ) .eq. 'POLY' )then
              y( 1 ) = 1
            else if( ctype( icmp ) .eq. 'KING' )then
              y( 1 ) = dfcns( 3, icmp )
            else
              call crash( 'ISOPSI', 'Unrecognized halo type' )
            end if
            y( 2 ) = 0
            ifail = 1
c integrate 2 ODEs to find Psi
            call ode_int( r0, r1, 2, y, tol, isopsf, ifail )
            psitab( i ) = y( 1 )
          end do
          iusedf = 1
          if( master )print *, 'ISOPSI: Table ready'
        end if
c look up value
        isopsi = aitken2( rtab, psitab, nrdd, r )
      else
c outside tidal radius
        isopsi = cmpmas( icmp ) * ( 1. / r - 1. / rtidal )
      end if
      return
      end
