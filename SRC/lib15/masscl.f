      subroutine masscl( normal )
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to zero out the mass array(s) or the expansion coefficients as
c   appropriate.
c The calling argument is important only for mulitple time step zones and
c   is used to determine the step number for that zone.  Calls from STEP
c   are normal, calls from MASSET are not
      use aarrays
      implicit none
c
c calling argument
      logical normal
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
c local variables
      integer i, ig, j, k, nzone
      real dt
c
c normal is `true' when called from STEP and `false' otherwise
      do nzone = 1, mzone
        if( normal )then
c the step number for mass assignment is the next step
          mstep( nzone ) = istep + nstep( nzone )
        else
c the step number for mass assignment is the current step
          mstep( nzone ) = 0
          if( istep .gt. 0 )mstep( nzone ) =
     +           ( ( istep - 1 ) / nstep( nzone ) + 1 ) * nstep( nzone )
        end if
c compute predicted position of grid center at the new time
        if( centrd )then
          dt = real( mstep( nzone ) - ipast( 1 ) ) / real( icstp )
          do ig = 1, ngrid
            do j = 1, 3
              xcpred( j, nzone, ig ) = cenfit( 1, j, ig ) * dt +
     +                                 cenfit( 2, j, ig ) * ( 1. - dt )
            end do
          end do
        end if
      end do
      if( normal )then
c clear acceleration components for off-grid particles if needed
        if( offfor .and. ( noff .gt. 0 ) )then
          do i = 1, noff
            do j = 1, 3
              accoff( j, i ) = 0
            end do
            potoff( i ) = 0
c destroy flags of particles that left the grid during the last step
            loco( i ) = abs( loco( i ) )
          end do
        end if
      end if
c save old mass arrays if zones > 1
      if( mzone .gt. 1 )then
        do nzone = 2, mzone
          do ig = 1, ngrid
            call switch( ig )
            k = mesh( jgrid )
            if( lsfp )then
              call blkcpy2( sfpmss( 1, nzone ), sfpmoz( 1, nzone ), k )
            else if( s3d )then
              call blkcpy(
     +                 s3dmss( 1, nzone, 1 ), s3dmoz( 1, nzone, 1 ), k )
            else if( lgrd )then
              call blkcpy( grdmss( 1, nzone ), grdmoz( 1, nzone ), k )
            else
              if( .not.( noslfg .or. bht ) )
     +                     call crash( 'MASSCL', 'Unrecognized method' )
            end if
            if( hybrid .and. ( jgrid .eq. 2 ) )then
              if( s3d )then
                call blkcpy(
     +                 s3dmss( 1, nzone, 2 ), s3dmoz( 1, nzone, 2 ), k )
              else
                call crash( 'MASSCL', '2nd grid not s3d' )
              end if
            end if
c update the step numbers of these arrays - do this only once per call
            if( jgrid .eq. ngrid )then
              lstep( 1, nzone ) = lstep( 2, nzone )
              lstep( 2, nzone ) = mstep( nzone )
            end if
          end do
        end do
      end if
c clear acceleration terms for perturber if present
      if( pertbn )call ptbscl( .false. )
c set all particle weights to unity, if they are equal (values could be
c   overwritten before each step)
      if( .not. uqmass )then
        do i = 1, mbuff
          pwt( i ) = 1
        end do
      end if
c clear mass arrays
      do ig = 1, ngrid
        call switch( ig )
c SFP methods
        if( lsfp )then
          if( sf2d .or. sf3d )then
            maxr = max( newmaxr, minmaxr )
            newmaxr = minmaxr
          end if
          ncontrib = 0
          do j = 1, mzone
            do i = 1, mesh( jgrid )
              sfpmss( i, j ) = 0
            end do
          end do
c PM+SH grid
        else if( s3d )then
          do j = 1, mzone
            do i = 1, mesh( jgrid )
              s3dmss( i, j, jgrid ) = 0
            end do
            if( hybrid )then
              if(
     +        jgrid .ne. 2 )call crash( 'MASSCL', 'hybrid grid mix up' )
              do i = 1, mesh( jgrid )
                s3dmss( i, j, 1 ) = 0
              end do
            end if
          end do
c all other grid methods
        else if( lgrd )then
          do j = 1, mzone
            do i = 1, mesh( jgrid )
              grdmss( i, j ) = 0
            end do
            if( p2d .or. p3d )then
              jrad( j ) = nr( jgrid )
              krad( j ) = 0
            end if
          end do
c direct-N method
        else if( dr3d )then
          do i = 1, ndrct
            do j = 5, 7
              drpt( j, i ) = 0
            end do
          end do
        else
          if( .not.( noslfg .or. bht ) )
     +                     call crash( 'MASSCL', 'Unrecognized method' )
        end if
      end do
      call switch( 0 )
      return
      end
