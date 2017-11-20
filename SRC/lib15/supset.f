      subroutine supset
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to determine and store the differences between the actual central
c   attraction of the model, as calculated on the grid, and the theoretically
c   expected function in which the DF would be in equilibrium.  These are
c   different for a number of reasons:
c      (a) softening and/or grid resolution,
c      (b) truncations or tapers,
c      (c) (in the case of discs) finite thickness, etc.
c The difference table is used to supplement the radial force acting on
c   every particle at every step.
c For disc components, the corrective forces are determined in the mid-plane
c   only, whereas spherical symmetry is assumed for halo components.
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
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
c externals
      logical gtlogl
      real meanrf, meansf
      real*8 frdisc, frhalo, frtot
c
c local variables
      integer i, ip, ir
      logical llocal( mcmp ), lovrd
      real actmas, am, corr, frad, frcfac, r, r1, r2, se, tfrad, u
c
c skip if nothing to do
      if( suppl )then
c determine which population is active
        ir = 0
        do i = 1, ncmp
          if( ( .not. testp( i ) ) .and. ( .not. rigidp( i ) ) )then
            ip = i
            ir = ir + 1
          end if
        end do
        if( ir .eq. 0 )call crash( 'SUPSET',
     +                                 'No active population detected' )
        if( ir .gt. 1 )then
          if( master )then
            print *, ir, ' populations active'
            print *, ' This routine can be used only if the' //
     +                                ' potentials of each are the same'
          end if
          if( ir .gt. 2 )call crash( 'SUPSET', 'Option not programmed' )
          lovrd = gtlogl( 'Do you want to continue' )
          if( .not. lovrd )call crash( 'SUPSET',
     +                    'Multiple populations option not programmed' )
        end if
        icmp = ip
c turn off any rigid masses and suppl for this part to avoid recursion
        do i = 1, ncmp
          llocal( i ) = rigidp( i )
          rigidp( i ) = .false.
        end do
        rigidh = .false.
        suppl = .false.
        if( lovrd )then
          icmp = 1
          actmas = cmpmas( 1 )
          cmpmas( 1 ) = 1
          cmpmas( 2 ) = 0
        end if
c choose radii for tabulated values of supplementary force
        r = rgrid( jgrid )
        if( c3d )r = .99 * r
        nrsup = mrsup - 1
        alp2 = log( r + 1. ) / real( nrsup - 1 )
        do ir = 1, nrsup + 1
          u = ir - 1
          htab( 1, ir ) = exp( alp2 * u ) - 1
        end do
c correct disc force for softening and truncations
        if( master .and. lprint )then
          write( no, * )
          write( no, * )'   radius th disc frc  actual frc  correction',
     +                  ' th halo frc  total force'
        end if
        frcfac = ts / gvfac
        do ir = 1, nrsup
          r = htab( 1, ir )
          if( disc( ip ) )then
            frad = meanrf( r ) / frcfac
          else
            frad = meansf( r ) / frcfac
          end if
          r = r / lscale
c theoretical radial force
          if( disc( ip ) )then
            tfrad = frdisc( dble( r ) )
          else
            tfrad = frhalo( dble( r ) )
          end if
c difference from theoretical force
          corr = tfrad - frad
c total radial force
          if( master .and. lprint )write( no, '( 3x, f6.3, 5f12.6 )' )r,
     +             frdisc( dble( r ) ), frad, corr, frhalo( dble( r ) ),
     +             frtot( dble( r ) )
          htab( 2, ir ) = frcfac * corr
        end do
c assume effective density drops to 0 at outer edge
        r1 = htab( 1, nrsup )
        r2 = htab( 1, nrsup + 1 )
        htab( 2, nrsup + 1 ) = htab( 2, nrsup ) * ( r1 / r2 )**2
c supplementary potential
        htab( 3, nrsup ) = r1 * htab( 2, nrsup )
        htab( 3, nrsup + 1 ) = htab( 3, nrsup ) * r1 / r2
c evaluate supplementary potential at interior points from work required
        do ir = 2, nrsup
          i = nrsup + 1 - ir
          r2 = htab( 1, i )
          htab( 3, i ) = htab( 3, i + 1 )
     +          - .5 * ( r2 - r1 ) * ( htab( 2, i ) + htab( 2, i + 1 ) )
          r1 = r2
        end do
c compute self energy of supplementary mass distribution
        se = 0.
        r1 = 0.
        do ir = 2, nrsup
          r = htab( 1, ir )
          r2 = r * r
          am = r1 * htab( 2, ir - 1 ) - r2 * htab( 2, ir )
          se = se + am * .5 * ( htab( 3, ir - 1 ) + htab( 3, ir ) )
          r1 = r2
        end do
c convert to external units and add to value that may have been set in halset
        selfe = selfe + .5 * se / ( lscale**5 * ts**4 )
c reset rigidh and supplementary force flag to apply these corrections
        rigidh = .false.
        do i = 1, ncmp
          rigidp( i ) = llocal( i )
          rigidh = rigidh .or. ( rigidp( i ) .and. ( .not. disc( i ) ) )
        end do
        suppl = .true.
c set flag to indicate table is ready
        lsupst = .true.
c restore original mass fractions
        if( lovrd )then
          cmpmas( 1 ) = actmas
          cmpmas( 2 ) = 1 - actmas
        end if
      end if
      return
      end
