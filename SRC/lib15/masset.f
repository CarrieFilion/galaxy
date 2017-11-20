      subroutine masset
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to work through all the particles and to assign their masses to the
c   grid or to compute the expansion coefficients.  Because mass assignment
c   is part of a normal time step, this routine is not called at every step.
c
c For SFP methods when the radial scale is flexible, the routine first
c   determines the radius of the smallest circle that just encloses all the
c   particles and then calls SFPSUM to evaluate the coefficients.  Otherwise,
c   the SFP coefficients can be determined in one pass.
c
c Note that the coordinates retrieved by GATHER are put into the oldc
c   array, while MASSGN uses the newc array.  A block copy from oldc to
c   newc is therefore required between calls to GATHER and MASSGN
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local variables
      integer jst, n
c
c set mass arrays to zero
      mzone = nzones
      call masscl( .false. )
c assign mass of all lists except that of particles off the grid
      do ilist = 1, nlists - 1
c set izone & jgrid
        call interpret
c work through groups of particles
        n = 1
        if( parallel )n = myid + 1
        inext = islist( 1, ilist, n )
        call switch( jlist )
        do while ( inext .ge. 0 )
          call gather( jst )
c copy data to new coordinates array
          call blkcpy( oldc( 1, 1 ), newc( 1, 1 ), 6 * mbuff )
c assign mass of particles or find maxr for SFP methods
          if( .not. noslfg )then
            if( heavies )then
              call switch( jlist )
              mhyb = jlist
              call massgn( jst )
            else
              mhyb = jgrid
              call massgn( jst )
              if( hybrid .and. ( jgrid .eq. 1 ) )then
                call switch( 2 )
                mhyb = 1
                call massgn( jst )
                call switch( 1 )
              end if
            end if
          end if
c compute accelerations on perturber if needed
          if( pertbn )call pertbz( jst )
        end do
      end do
      call switch( 0 )
c may need to work through particles a second time for some SFP methods
      if( sf2d .or. sf3d )then
        call sfpsum
      end if
c combine masses
      call mascmb
c check particle count
      call ncheck
c rescale density array
      call scaled
      return
      end
