      program dfiter
c  Copyright (C) 2014, Jerry Sellwood
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'inc/setup.f'
c
      real aold, potold, zold
      logical recalc
      common / potsav / recalc, aold, zold, potold
c
c external
      logical gtlogl
c
c local variables
      integer i
      real crrnt, target
      parameter ( target = 1. )
      include 'inc/pi.f'
c
      data ni, no / 7, 8 /
c
      call setnag
c open data file and determine grid type
      call getset
      call boilrp( .true. )
c decide whether this is a cold or a warm start
      print *, 'You will usually want to select a cold start'
      print *, 'A warm start will try to refine an existing .pft file'
      call gtintg( 'Enter 1 for a cold start, 0 for a warm start',
     +                                                          iusepf )
c restart .lis file for cold starts only
      if( iusepf .eq. 1 )rewind no
      iusepf = -iusepf
c
c setting up section
      call msetup
      call grdset
c restrict to axial symmetry
      call restrict
      call setrun
      close( ni )
c multiple zones not allowed in dfiter
      nzones = 1
      isfld = 0
c identify halo population
      icmp = 0
      do i = 1, ncmp
        if( ( .not. disc( i ) ) .and.
     +      ( cdft( i ) .eq. 'DFIT' ) )then
          if( icmp .eq. 0 )then
            icmp = i
          else
            print *, 'more than one halo component', icmp, i
            call gtintg( 'Choose which for this DFITER run', icmp )
          end if
        end if
      end do
      if( icmp .eq. 0 )call crash( 'DFITER', 'No DFIT pop found' )
      sphrod( icmp ) = .not. gtlogl( 'Is this a spherical halo' )
c set default normalization
      dfcns( 2, icmp ) = 0
      if( impar( icmp ) .eq. 1 )then
        dfcns( 2, icmp ) = 3. / ( 7. * pi**3 )
      else if( impar( icmp ) .eq. 2 )then
        dfcns( 2, icmp ) = dfnorm( icmp )
      else if( impar( icmp ) .eq. 4 )then
        dfcns( 2, icmp ) = 1
      end if
      if( dfcns( 2, icmp ) .eq. 0. )call crash(
     +                             'DFITER', 'Unrecognized DF formula' )
c SFP methods
      if( sf2d .or. sf3d )then
        lsfpr = nrtab
        lsfpl = lsfpr * maxn * maxn
        allocate ( sfprad( lsfpr ) )
        allocate ( sfplst( 2, lsfpl ) )
        call sfptab
      else if( lgrd )then
c create or check Greens function file for grid methods
        call greenm( .false. )
      end if
c
c first guess at potential is analytic
      if( analyt )fmass( icmp ) = 1
c allocate space
      call setspc
c begin iteration
      crrnt = 0
      call switch( igrd( icmp ) )
      do while ( abs( crrnt - target ) .gt. 0.00001 )
c assign mass to largest grid
        call denset( crrnt )
c solve for new gravitational field - not necessary for s3d in hybrid methods
        if( .not. hybrid )call findf( .true. )
c reset flags for next iteration
        analyt = .false.
        recalc = .true.
c find spline fit and save
        iusepf = 0
        aold = -100
        call splfit
        print '( ''Current values'', f8.3, 3f12.8 )',
     +                              ( dfcns( i, icmp ), i = 1, 3 ), rmax
      end do
c converged
      print '( ''Final values'', f8.3, 2f12.8 )',
     +                                    ( dfcns( i, icmp ), i = 1, 3 )
      end

      subroutine restrict
      implicit none
c routine to override instructions in .DAT file by imposing axial symmetry
c   for the purposes of dfiter only
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
c local variables
      integer i, m, n
c
c set active functions
      if( sf2d .or. sf3d )then
        lastf = 0
        if( minm .gt. 0 )call
     +     crash( 'restrict', 'fixed axisymmetric forces!' )
        maxm = 0
        m = 0
        do n = minn, maxn
          lastf = lastf + 1
          msel( lastf ) = m
          nsel( lastf ) = n
        end do
      else if( .not. s3d )then
c set Fourier filter to skip all m>0 (Note: m = i - 1)
        do i = 2, ng
          lg( i ) = .true.
        end do
      end if
      return
      end
