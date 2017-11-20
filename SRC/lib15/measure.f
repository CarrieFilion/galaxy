      subroutine measure( after,
     +                    nact, lprop, prop, nbplot, dispy, lglen,
     +                    alspi, jlen, besc, nLbin, nLa, lzman, zman,
     +                    lzpr, zpbin, lmoni, ehmon, nwring, wring,
     +                    llz, Lzval )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c performs various initialization and completion tasks for analysis of
c   the paricles during one step of the simulation.  It is called from
c   STEP both before and after the particles have been advanced.
c
c the number of active populations is set by ncmp, which is ncmp - 1 only
c   when rings are active
c
c calling arguments
      integer jlen, lglen, llz, lmoni, lprop, lzman, lzpr, nact, nbplot
      integer nLbin, nwring
      logical after
      integer nLa( nact, nLbin )
      real alspi( nact, 2, lglen ), besc( nact, 2, jlen )
      real dispy( nbplot ), ehmon( 6, lmoni ), Lzval( llz )
      real prop( nact, lprop )
      real wring( nwring ), zman( nact, 2, lzman ), zpbin( nact, lzpr )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local variables
      character*4 astr
      integer i, im, ip, j, k, n
      real ahol, amsfac
c
      if( parallel )call crash( 'MEASURE', 'Parallel version needed' )
      if( after )then
c finish up
        if( plot )then
          call disply( -1, i, nbplot, dispy )
        end if
c
        if( lval .and. master )then
          astr = 'LVAL'
          read( astr, '( a4 )' )ahol
          i = 1
          write( nphys )irun, ahol, istep, nbod, i
          write( nphys )( Lzval( i ), i = 1, nbod )
        end if
c
        call phyprp( nact, lprop, prop, lglen, alspi, jlen, besc,
     +               nLbin, nLa, lzman, zman, lzpr, zpbin,
     +               lmoni, ehmon, nwring, wring )
      else
c get radial density profile
        if( rhor )then
          call rhoprf
        end if
        if( lval )then
          do i = 1, llz
            Lzval( i ) = 0
          end do
        end if
c initialise energy check variables
        do j = 1, ncmp
          do i = 1, 3
            cm( i, j ) = 0
            pp( i, j ) = 0
            for( i, j ) = 0
            ang( i, j ) = 0
          end do
          nspop( j ) = 0
          popm( j ) = 0
          pe( j ) = 0
          ke( j ) = 0
        end do
        claus = 0
c convert angular momentum lost since last analysis to natural units and add
        amsfac = pmass / ( lscale**5 * ts**3 )
        do i = 1, 2
          do j = 1, 3
            angoff( j, i ) = angoff( j, i ) + amsfac * amoff( j, i )
            amoff( j, i ) = 0
          end do
        end do
c initialise angular momentum array
        if( angm )then
          k = maxLz + 2
          do i = 1, k
            do j = 1, ncmp
              nLa( j, i ) = 0
            end do
          end do
        end if
c initialise velocity field arrays
        if( vfld .or. vflh )then
          do i = 1, lprop
            do j = 1, nact
              prop( j, i ) = 0
            end do
          end do
        end if
c initialise logarithmic spiral arrays
        if( lgsp )then
          j = 0
          do ip = 1, np
            do im = 1, nm
              j = j + 1
              do i = 1, nact
                alspi( i, 1, j ) = 0.
                alspi( i, 2, j ) = 0.
              end do
            end do
          end do
        end if
c initialise integral monitoring arrays
        if( moni )then
          n = 6
          if( threed )n = 4
          do i = 1, nmonit
            do j = 1, n
              ehmon( j, i ) = 0
            end do
          end do
        end if
c initialise zmean analysis arrays
        if( zanl )then
          j = 0
          do ip = 1, nrz
            do im = 1, nmz + 1
              j = j + 1
              do i = 1, ncmp
                zman( i, 1, j ) = 0.
                zman( i, 2, j ) = 0.
              end do
            end do
          end do
        end if
c initialise spherical Bessel function arrays
        if( sphb )then
          do j = 1, jlen
            do i = 1, nact
              besc( i, 1, j ) = 0.
              besc( i, 2, j ) = 0.
            end do
          end do
          do i = 1, ncmp
            npbess( i ) = 0
          end do
        end if
c initialise z density profile array
        if( zprf )then
          k = nbzr * nbzz
          do j = 1, k
            do i = 1, ncmp
              zpbin( i, j ) = 0.
            end do
          end do
        end if
c clear rings array
        if( rngs )then
          do j = 1, nwring
            wring( j ) = 0
          end do
        end if
c this is normally done in masscl, but that routine is now called later
        if( sf2d .or. sf3d )maxr = max( newmaxr, minmaxr )
c compute and save moment of inertia tensor
        if( moit )then
          call getmoi
        end if
c compute and save frequencies
        if( frqs )then
          call freqs
        end if
c initialise plot file creation if required
        if( plot )then
          call disply( 0, i, nbplot, dispy )
        end if
      end if
      return
      end
