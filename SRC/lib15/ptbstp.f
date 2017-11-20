      subroutine ptbstp
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c advances the motion of an external perturber using the previously summed
c   reaction forces from the active particles.  Saves information to the
c   .lis and .res files as appropriate
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
c local variables
      character*4 bstr
      integer i
      real ahol, pmfac
c
      data bstr / 'SATE' /
c
c check whether a perturber is present
      if( pertbn )then
        call ptbscl( .true. )
        pmfac = pmass / ( lscale**3 * ts**2 )
        if( cvers .lt. 9.115 )pmfac = pmfac * real( nsect )
c set whether to shift grid
        gshift = .false.
c motion of perturber makes sense only if there is no imposed rotnl symm
        if( nsect .eq. 1 )then
c compute grid shift to adjust for CoM motion of galaxy relative to perturber
          if( gshift )then
            do i = 1, 3
              vgal( i ) = -mptbr * ( vptrb( i ) + .5 * accptb( i ) ) /
     +                                          ( pmfac * real( nbod ) )
              xshift( i ) = -vgal( i ) * ts * lscale
            end do
          else
            do i = 1, 3
              xshift( i ) = 0
            end do
          end if
        end if
c save data if requested
        if( phys .and. master )then
          read( bstr, '( a4 )' )ahol
          if( exsphr .or. exgnrc )then
            peptrb = peptrb * pmfac
            i = 10
            write( nphys )irun, ahol, istep, i, nptrb
            if( lfrv )then
              write( nphys )xptrb,
     +              ( vptrb( i ) + .5 * accptb( i ) * ts, i = 1, 3 ),
     +              accptb, peptrb
            else
              write( nphys )( xptrb( i ) - .5 * vptrb( i ) * ts,
     +                                 i = 1, 3 ), vptrb, accptb, peptrb
            end if
          else if( extbar )then
            i = 3
            write( nphys )irun, ahol, istep, i, nptrb
            if( lfrv )then
              write( nphys )bphase, omegab + .5 * ts * accptb( 1 ),
     +                      accptb( 1 )
            else
              write( nphys )bphase - .5 * ts * omegab, omegab,
     +                      accptb( 1 )
            end if
          else if( exspir )then
            i = 2
            write( nphys )irun, ahol, istep, i, nptrb
            if( lfrv )then
              write( nphys )sphase, omegsp + .5 * ts * accptb( 1 )
            else
              write( nphys )sphase - .5 * ts * omegsp, omegsp
            end if
          else
            call crash( 'PTBSTP', 'Unknown perturber' )
          end if
        end if
c leap-frog for perturber motion in external units
        if( exsphr .or. exgnrc )then
c motion of perturber makes sense only if there is no imposed rotnl symm
          if( nsect .eq. 1 )then
            do i = 1, 3
              vptrb( i ) = vptrb( i ) + accptb( i ) * ts
              xptrb( i ) = xptrb( i ) + vptrb( i ) * ts +
     +                                              xshift( i ) / lscale
            end do
          end if
        else if( extbar )then
c bar center fixed at origin, but it is torqued
          omegab = omegab + ts * accptb( 1 )
          bphase = bphase + ts * omegab
        else if( exspir )then
c spiral simply rotates at constant speed
          sphase = sphase + ts * omegsp
        else
          call crash( 'PTBSTP', 'Unknown perturber 2' )
        end if
      end if
      return
      end
