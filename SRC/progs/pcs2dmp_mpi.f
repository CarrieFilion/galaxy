      program pcs2dmp
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c MPI version of the program to read a .pcs file containing particle
c   positions, velocities, and optionally masses.  It scales the input
c   coordinates to internal (program) units, divides particles into their
c   various time step zones, time centers the initial velocities, and
c   initializes the linked lists.  The final step is to output all the
c   information needed for GALAXY into the restart file (.dmp)
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
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'mpif.h'
c
c local variables
      integer i, j
      logical gtlogl, lcrash
c
      call setnag
c read run number from standard input and open main ASCII I/O files
      call getset
      if( bht .or. dr3d )call crash( 'PCS2DMP_MPI',
     +                         'This method has not been parallelized' )
      call boilrp( .true. )
c verify no resfil exists and offer the opportunity to abort if it does
      j = -1
      call opnfil( j, 'res', 'unformatted', 'old', 'seq', i )
      if( i .eq. 0 )then
        lcrash = .not.
     +             gtlogl( '.res file exists, do you want to continue' )
        if( lcrash
     +          )call crash( 'PCS2DMP_MPI', '.res file already exists' )
      end if
      close( j )
c read .dat file (ASCII input)
      call msetup
      call grdset
      call setrun
      close( ni )
c SFP methods
      if( sf2d .or. sf3d )then
        lsfpr = 1001
        lsfpl = lsfpr * maxn * maxn
        allocate ( sfprad( lsfpr ) )
        allocate ( sfplst( 2, lsfpl ) )
        call sfptab
      else if( lgrd )then
c check or create data files or tables only if needed
        if( lfrv )call greenm( .false. )
      end if
c allocate space
      call setspc
c read in the particles from the file
      call loadup
c read supplementary forces table if needed
      if( suppl )call suptab( .true. )
      if( master )then
        print *, 'Run', irun, ' started at step no', istep
c open old results file
        call opnfil( nphys, 'res', 'unformatted', 'old', 'append', i )
c open and write header on new results file if previous open failed
        if( i .ne. 0 )then
          call opnfil( nphys, 'res', 'unformatted', 'new', 'seq', i )
          call hedrec( nphys, .false. )
        end if
      end if
c find centroid if requested
      if( centrd )call newxcn( .true. )
c allocate space for direct methods
      if( bht )then
        bhmt = 5 * nbod
        bhnlev = 30
        allocate ( bhbox( bhnlev ) )
        allocate ( bhcom( 4, bhmt ) )
        allocate ( itup( bhmt ) )
        allocate ( itdown( bhmt / 8 ) )
      else if( dr3d )then
c find total of all particles to be computed by dr3d
        ndrct = 0
        do i = 1, ncmp
          j = igrd( i )
          if( igrid( j ) .eq. 10 )ndrct = ndrct + nsp( i )
        end do
        allocate ( drpt( 15, ndrct ) )
        call masset
      end if
      lprint = master
      phys = mod( istep, iphys ) .eq. 0
c get forces if needed
      if( lfrv )then
        call masset
        call findf( .false. )
      end if
c time-centre coordinates
      call tmcent( .true. )
c assign mass to grids - findf may have used the array
      call masset
c save initial set up
      call dumpp
c finish up
      call mpi_finalize( i )
      stop
      end

      subroutine fiddle
c  Copyright (C) 2014, Jerry Sellwood
c
c this routine, which is never executed, is to force the linker to
c   to load mpi versions of these routines
c
      call chklst
      call centre
      call hedrec( 1, .true. )
      call gtintg( 'a', i )
      call mascmb
      call ncheck
      call p2fndf( 1, 1 )
      call p3fndf( 1, 1 )
      call ptbscl
      call qudini
      call ranuni
      call second( 0. )
      return
      end
