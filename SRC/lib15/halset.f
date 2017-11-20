      subroutine halset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to determine the proper potential and self-energy of a rigid halo
c   from a line integral of the radial force.  The density of an extended
c   halo is assumed to be cut off sharply at the grid edge, which allows
c   the zero of the potential to be set at infinity.
c The tabulated halo potential values are saved in internal program units
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
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
c external
      real*8 frhalo
c
c local variables
      integer i, ir, jcmp
      real am, frcfac, r, r1, r2, tfrad, u
c
c set rigidh
      rigidh = .false.
      do i = 1, ncmp
        rigidh = rigidh .or. ( rigidp( i ) .and. ( .not. disc( i ) ) )
      end do
c nothing else to do if no rigid halo or if total radial force is analytic
      if( ( .not. rigidh ) .or. fixrad )then
        selfe = 0
      else
c choose radii for tabulated values of supplementary force
c spacing independent of the grid
        r = rgrid( jgrid )
        nrsup = mrsup - 1
        alp2 = log( r + 1. ) / real( nrsup - 1 )
        do ir = 1, nrsup
          u = ir - 1
          htab( 1, ir ) = exp( alp2 * u ) - 1
        end do
c work over all populations
        jcmp = icmp
        frcfac = ts / gvfac
        do ir = 1, nrsup
          r = htab( 1, ir ) / lscale
c theoretical radial force
          tfrad = 0
          do icmp = 1, ncmp
            if( rigidp( icmp ) )then
              if( cmprssd )then
                icomp = 0
                do i = 1, ncomp
                  if( iccmp( i ) .eq. icmp )icomp = i
                end do
              end if
              tfrad = tfrad + frhalo( dble( r ) )
            end if
          end do
          htab( 2, ir ) = frcfac * tfrad
        end do
c restore original icmp
        icmp = jcmp
c assume effective density drops to 0 at outer edge
        r1 = rgrid( jgrid )
        u = nrsup
        r2 = exp( alp2 * u ) - 1
        htab( 2, nrsup + 1 ) = htab( 2, nrsup ) * ( r1 / r2 )**2
c evaluate potential at and one space beyond grid edge
        hpot( nrsup ) = r1 * htab( 2, nrsup )
        hpot( nrsup + 1 ) = hpot( nrsup ) * r1 / r2
c evaluate halo potential at interior points from work required
        rigidh = .true.
        do ir = 2, nrsup
          i = nrsup + 1 - ir
          r2 = htab( 1, i )
          hpot( i ) = hpot( i + 1 )
     +          - .5 * ( r2 - r1 ) * ( htab( 2, i ) + htab( 2, i + 1 ) )
          r1 = r2
        end do
c compute self energy of halo mass distribution
        selfe = 0.
        r1 = 0.
        do ir = 2, nrsup
          r = htab( 1, ir )
          r2 = r * r
          am = r1 * htab( 2, ir - 1 ) - r2 * htab( 2, ir )
          selfe = selfe + am * .5 * ( hpot( ir - 1 ) + hpot( ir ) )
          r1 = r2
        end do
c convert to external units
        selfe = .5 * selfe / ( lscale**5 * ts**4 )
        if( master )write( no,
     +                 '( ''Self energy of the halo ='', f10.6 )' )selfe
      end if
      return
      end
