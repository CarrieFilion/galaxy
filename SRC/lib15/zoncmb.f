      subroutine zoncmb
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to sum masses in different zones on the grid prior to force
c   determination
      use aarrays
      implicit none
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
      integer i, ig
      integer nzone
      real fract !, e
c
c update the step numbers of new arrays noslfg only - usually done in masscl
      if( noslfg )then
        if( mzone .gt. 1 )then
          do nzone = 2, mzone
            if( lstep( 2, nzone ) .lt. mstep( nzone ) )then
              lstep( 1, nzone ) = lstep( 2, nzone )
              lstep( 2, nzone ) = mstep( nzone )
            end if
          end do
        end if
      end if
c number of zones for next step
      if( nzones .gt. 1 )then
        do i = 1, nzones
          if( mod( istep, nstep( i ) ) .eq. 0 )mzone = i
        end do
      else
        mzone = 1
      end if
c
c combine masses in different zones by linear interpolation between times
c
      if( ( .not. noslfg ) .and. ( nzones .gt. 1 ) )then
        do ig = 1, ngrid
          call switch( ig )
          do nzone = 2, nzones
c check most recent mass array (usually needed)
            i = lstep( 2, nzone ) - istep
            if( ( i .lt. 0 ) .or. ( i .gt. nstep( nzone ) ) )then
              if( master )write( no, * )istep, lstep( 2, nzone )
              call crash( 'ZONCMB', 'New mass array at a useless time' )
            end if
            if( i .ne. nstep( nzone ) )then
              fract = 1. - real( i ) / real( nstep( nzone ) )
              if( master .and. lprint .and. ( jgrid .eq. 1 ) )
     +                   write( no, 200 )fract, nzone, lstep( 2, nzone )
c add fraction of most recent masses to current masses
              if( lsfp )then
                do i = 1, mesh( jgrid )
                  sfpmss( i, 1 ) =
     +                       sfpmss( i, 1 ) + fract * sfpmss( i, nzone )
                end do
              else if( s3d )then
                do i = 1, mesh( jgrid )
                  s3dmss( i, 1, 1 ) =
     +                 s3dmss( i, 1, 1 ) + fract * s3dmss( i, nzone, 1 )
                end do
                if( hybrid .and. ( jgrid .eq. 2 ) )then
                  do i = 1, mesh( jgrid )
                    s3dmss( i, 1, 2 ) =
     +                 s3dmss( i, 1, 2 ) + fract * s3dmss( i, nzone, 2 )
                  end do
                end if
              else
                do i = 1, mesh( jgrid )
                  grdmss( i, 1 ) =
     +                       grdmss( i, 1 ) + fract * grdmss( i, nzone )
                end do
c end if block for methods
              end if
c end if for later fragments
            end if
c check older mass array if needed
            i = lstep( 2, nzone ) - istep
            if( i .ne. 0 )then
              i = istep - lstep( 1, nzone )
              if( ( lstep( 1, nzone ) .lt. 0 ) .or.
     +            ( i .lt. 0 ) .or. ( i .gt. nstep( nzone ) ) )then
                if( master )write( no, * )istep, lstep( 1, nzone )
                call crash( 'ZONCMB',
     +                              'Old mass array at a useless time' )
              end if
              fract = 1. - real( i ) / real( nstep( nzone ) )
              if( master .and. lprint .and. ( jgrid .eq. 1 ) )
     +                   write( no, 200 )fract, nzone, lstep( 1, nzone )
c add fraction of older masses to current masses
              if( lsfp )then
                do i = 1, mesh( jgrid )
                  sfpmss( i, 1 ) =
     +                       sfpmss( i, 1 ) + fract * sfpmoz( i, nzone )
                end do
              else if( s3d )then
                do i = 1, mesh( jgrid )
                  s3dmss( i, 1, 1 ) =
     +                 s3dmss( i, 1, 1 ) + fract * s3dmoz( i, nzone, 1 )
                end do
                if( hybrid .and. ( jgrid .eq. 2 ) )then
                  do i = 1, mesh( jgrid )
                    s3dmss( i, 1, 2 ) =
     +                 s3dmss( i, 1, 2 ) + fract * s3dmoz( i, nzone, 2 )
                  end do
                end if
              else
                do i = 1, mesh( jgrid )
                  grdmss( i, 1 ) =
     +                       grdmss( i, 1 ) + fract * grdmoz( i, nzone )
                end do
c end if block for methods
              end if
c end if for earlier fragments
            end if
c end loop over zones
          end do
c end loop over grids
        end do
      end if
      call switch( 0 )
      return
  200 format( ' Adding fraction', f5.2, ' of masses from zone no', i2,
     +        ' step', i7 )
      end
