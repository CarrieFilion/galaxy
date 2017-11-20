      subroutine p2fndf( jzone, kzone )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to calculate the force field and potential of a single zone of
c   the 2-D polar grid
      use aarrays
      implicit none
c
c calling arguments
      integer jzone, kzone
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local variables
      integer ibase, ipm, ispac, itype, jrf, jrs, krf, krs, ntypes
c
c determine radial range which contains all sources
      jrs = nr( jgrid )
      krs = 0
      do ipm = 1, nzones
        jrs = min( jrs, jrad( ipm ) )
        krs = max( krs, krad( ipm ) )
        if( ipm .gt. 1 )krad( ipm ) = max( krad( ipm ), krs )
      end do
      if( ( jrs .gt. krs ) .or. ( krs .eq. 0 ) )then
        if( master )write( no, * )' No mass assigned to grid?', jrs, krs
        call crash( 'P2FNDF', 'Invalid radial range for source masses' )
      end if
      ispac = ( krs - jrs + 1 ) * na
      ibase = na * ( jrs - 1 ) + 1
c find Fourier coefficients of mass distribution already combined into zone 1
      call p2manl( 1, jrs, krs )
c assume forces are required over the whole grid unless told otherwise
      jrf = 1
      if( ( jzone .gt. 0 ) .and.
     +    ( jrad( jzone ) .lt. nr( jgrid ) ) )jrf = jrad( jzone )
      krf = nr( jgrid )
      if( ( kzone .le. nzones ) .and.
     +    ( .not. wholeg ) )then
        krf = krad( kzone )
        do while ( krf .lt. jrf )
          kzone = kzone + 1
          if(
     +       kzone .gt. nzones )call crash( 'P2FNDF', 'kzone > nzones' )
          krf = krad( kzone )
        end do
      end if
c work over types
      ntypes = 2
      if( potl )ntypes = 3
      do itype = 1, ntypes
c copy mass coeffs to target area
        call blkcpy( grdmss( ibase, 1 ), grdfld( ibase, itype ), ispac )
c convolve with Green function
        call p2conv( itype, jrs, krs, jrf, krf )
c re-synthesize
        call p2fsyn( itype, jrf, krf )
      end do
      return
      end
