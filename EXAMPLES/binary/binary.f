      program binary
c  Copyright (C) 2015, Jerry Sellwood
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
      include 'inc/lunits.f'
c
      include 'inc/setup.f'
c
c local variables
      character*80 line
      integer i, ilast
      logical smooth

c$$$      integer j, k, ind( 16 )
c$$$      common / local / ind
c$$$      include 'inc/grids.f'
c$$$      real gridac, gridms

c
      data smooth / .false. /
c
c read run number from standard input and open main ASCII I/O files
      call getset
c rewind .lis file - needed only if there was a previous false start
      rewind no
c read .dat file (ASCII input)
      call msetup
      call grdset
      call setrun
      call getline( ni, line )
      read( line, '( 10x, i10 )' )ilast
      close( ni )
c check or create data files or tables for field method
      call greenm( .false. )
c place particles in their initial positions
      phys = .true.
      lprint = .true.
c reserve space
      call setspc
c place particles in initial positions and (possibly) set some vel components
      call psetup
c open results file
      call opnfil( nphys, 'res', 'unformatted', 'unknown', 'seq', i )
      call opnfil( 33, 'frc', 'unformatted', 'unknown', 'seq', i )
c rewind and write header record on results file
      call hedrec( nphys, .false. )
c determine field of disc particles
      phys = .false.
      call masset
      phys = .true.
c determine gravitational field of a smooth (or the initial) mass distribution
      call findf( phys )
      call adjvls( .true. )
c time-centre velocities
      call tmcent( .true. )
      call masset
c set up flag file
      nflg = 3
      call opnfil( nflg, 'flg', 'formatted', 'unknown', 'seq', i )
      write( nflg, * )ilast
      close( unit = nflg )

c$$$      call jsbgn
c$$$      call jspage
c$$$c      call jscale( 0., 5., -.0001, .0001 )
c$$$      call jscale( 0., 5., -.1, .1 )
c$$$      call jsaxis( 'x', 'z', 1 )
c$$$      call jsaxis( 'y', '\\Deltaa_z', 1 )

c main time-step sequence
    1 lprint = .false.
c    1 lprint = .true.

      print *, 'starting step no', istep
c$$$      if( p3d )then
c$$$        j = 0
c$$$        do i = 1, mesh( 1 )
c$$$          if( gridms( i ) .ne. 0. )then
c$$$            j = j + 1
c$$$            print *, i, j, gridms( i )
c$$$            ind( j ) = i
c$$$          end if
c$$$        end do
c$$$      end if

c calculate forces
      call findf( phys )

c$$$      if( p3d )then
c$$$        print *, 'radial forces'
c$$$        do j = 1, 16
c$$$          i = ind( j )
c$$$          print *, i, j, gridac( i, 1 )
c$$$        end do
c$$$        read *, i
c$$$        print *, 'azimuthal forces'
c$$$        do j = 1, 16
c$$$          i = ind( j )
c$$$          print *, i, j, gridac( i, 2 )
c$$$        end do
c$$$        read *, i
c$$$        print *, 'vertical forces'
c$$$        do j = 1, 16
c$$$          i = ind( j )
c$$$          print *, i, j, gridac( i, 3 )
c$$$        end do
c$$$        read *, i
c$$$        print *, 'potentials'
c$$$        do j = 1, 16
c$$$          i = ind( j )
c$$$          print *, i, j, gridac( i, 4 )
c$$$        end do
c$$$        read *, i
c$$$      end if

c move particles
      call step
      if( phys .and. master )then
        close( unit = nphys )
        call opnfil( nphys, 'res', 'unformatted', 'old', 'append', i )
      end if
      phys = mod( istep, iphys ) .eq. 0
c check whether to continue
      if( mod( istep, 10 ) .eq. 1 )then
        call opnfil( nflg, 'flg', 'formatted', 'old', 'seq', i )
        if( i .eq. 1 )stop
        read( nflg, * )ilast
        close( unit = nflg )
      end if
c decide whether next step is a main or a sub-step
      if( mod( istep, nstep( nzones ) ) .eq. 0 )then
        write( no, * )'starting step no ', istep
      end if
c start next step
      if( istep .le. ilast )go to 1

c$$$      call jsend

      stop
      end

c$$$      real function gridms( i )
c$$$      use aarrays
c$$$      implicit none
c$$$c
c$$$      integer i
c$$$c
c$$$      include 'inc/params.f'
c$$$c
c$$$      include 'inc/admin.f'
c$$$c
c$$$      include 'inc/grids.f'
c$$$c
c$$$      if( s3d )then
c$$$        gridms = s3dmss( i, 1, 1 )
c$$$      else if( lgrd )then
c$$$        gridms = grdmss( i, 1 )
c$$$      else
c$$$        call crash( 'GRIDMS', 'Unrecognized grid' )
c$$$      end if
c$$$      return
c$$$      end
c$$$
c$$$      real function gridac( i, j )
c$$$      use aarrays
c$$$      implicit none
c$$$c
c$$$      integer i, j
c$$$c
c$$$      include 'inc/params.f'
c$$$c
c$$$      include 'inc/admin.f'
c$$$c
c$$$      include 'inc/grids.f'
c$$$c
c$$$      if( s3d )then
c$$$        gridac = s3dfld( i, j )
c$$$      else if( lgrd )then
c$$$        gridac = grdfld( i, j )
c$$$      else
c$$$        call crash( 'GRIDAC', 'Unrecognized grid' )
c$$$      end if
c$$$      return
c$$$      end

      include 'getacc.f'
      include 'grdacc.f'
      include 'icheck.f'
      include 'psetup.f'
      include 'p3dset.f'
