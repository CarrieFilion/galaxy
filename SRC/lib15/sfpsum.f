      subroutine sfpsum
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to evaluate the coefficients of the SFP basis expansion
c
c It is needed for only those SFP bases for which the radial scale can be
c   varied - primarily the Abel-Jacobi functions - for which previous calls
c   to MASSGN will have determined the new value of maxr.  In the usual case
c   in which the radial scale is kept fixed, the coefficients will have been
c   evaluated already by calls to MASSGN and this routine does nothing.
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
c local variable
      integer jst, n
c
      if( basset .eq. 'ablj' )then
c revise maxr for next step
        maxr = min( newmaxr, maxmaxr )
        n = 1
        if( parallel )n = myid + 1
c work through all lists of particles except those off the grid
        do ilist = 1, nlists - 1
c find zone and grid code for this list
          call interpret
          if( izone .le. mzone )then
            inext = islist( 1, izone, n )
            do while ( inext .ge. 0 )
              call gather( jst )
c copy data to new coordinates array
              call blkcpy( oldc( 1, 1 ), newc( 1, 1 ), 6 * mbuff )
c set radii, cell numbers and weights
              call weight( jst, .false. )
c add in these particles
              call sfpadd( jst )
            end do
          end if
        end do
      else
c nothing to do for these bases
        if( .not. ( ( basset .eq. 'bess' ) .or.
     +              ( basset .eq. 'lgsp' ) ) )call crash( 'SFPSUM',
     +                                        'Unrecognized basis set' )
      end if
      return
      end
