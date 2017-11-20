      program galaxy
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c MPI version of main driving program for simulation of an isolated galaxy
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
c A variety of methods are available to determine the gravitational field
c   from the particles, including the possibility of composite combinations
c   of some methods
c
c The simulation must first be set up from a .pcs file by running pcs2dmp_mpi,
c   after which this program can be executed as often as desired
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
      include 'mpif.h'
c
c external
      character elapsd*24
c
c local variables
      character*80 line
      integer i, ierr, ilast, jrun
      real t, t2( 2 ), etime
c
      call setnag
c read run number from standard input and open main ASCII I/O files
      call getset
      if( bht .or. dr3d )call crash( 'GALAXY_MPI',
     +                         'This method has not been parallelized' )
      call boilrp( .true. )
      t = etime( t2 )
c read .dat file (ASCII input)
      call msetup
      call grdset
      call setrun
      call getline( ni, line )
      read( line( 11:40 ), * )ilast
c ensure last step is for a top zone
      i = nstep( nzones )
      ilast = i * ( ( ilast - 1 ) / i + 1 )
      close( ni )
      call mpi_barrier( mpi_comm_world, ierr )
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
c allocate space
      call setspc
c read the .dmp file
      jrun = irun
      call restrt
      if( irun .ne. jrun )then
        if( master )print *, '.dmp file from run', irun,
     +                                    ' but .dat file for run', jrun
        call crash( 'GALAXY', 'Files muddled' )
      end if
      if( master )then
        print *, 'Run', irun, ' restarted at step no', istep
        print *, 'Run will stop after step ', ilast
c open old results file
        call opnfil( nphys, 'res', 'unformatted', 'old', 'append', i )
c open and write header on new results file if previous open failed
        if( i .ne. 0 )then
          call opnfil( nphys, 'res', 'unformatted', 'new', 'seq', i )
          call hedrec( nphys, .false. )
        end if
c set up flag file
        call opnfil( nflg, 'flg', 'formatted', 'unknown', 'seq', i )
        write( nflg, * )ilast
        close( nflg )
      end if
      call mpi_barrier( mpi_comm_world, ierr )
      lprint = master
      phys = mod( istep, iphys ) .eq. 0
c main time-step sequence
      do while ( istep .le. ilast )
        lprint = .false.
        lprint = lprint .and. master
        t = etime( t2 )
        if( master )print *,
     +                    'Starting step', istep, ' after ', elapsd( t )
c save density information before it is used
        if( phys )call densty
c calculate gravitational field
        call findf( phys )
c move particles
        call step
        if( phys .and. master )then
          close( unit = nphys )
          call opnfil( nphys, 'res', 'unformatted', 'old', 'append', i )
        end if
c check whether to continue (deleting runxxx.flg causes a graceful stop)
        i = 10
        if( nstep( nzones ) .gt. 1 )i = nstep( nzones )
        if( mod( istep, i ) .eq. 0 )then
          if( master )then
            call opnfil( nflg, 'flg', 'formatted', 'old', 'seq', i )
            if( i .eq. 0 )then
              read( nflg, *, iostat = i )ilast
              close( unit = nflg )
              if( i .ne. 0 )ilast = istep - 1
            else
              ilast = istep - 1
            end if
c ensure last step is for a top zone
            i = nstep( nzones )
            ilast = i * ( ( ilast - 1 ) / i + 1 )
          end if
          call mpi_bcast( ilast, 1, mpi_integer, 0, mpi_comm_world, i )
        end if
c decide whether next step is a main or a sub-step
        phys = mod( istep, iphys ) .eq. 0
        if( ( mod( istep, nstep( nzones ) ) .eq. 0 ) .and.
     +        master )write( no, * )'starting step no ', istep
c update dump file every 25 analysis steps or if this is the last chance
        if( ( mod( istep, 25 * iphys ) .eq. 0 ) .or.
     +      ( istep .eq. ilast ) )call dumpp
      end do
      if( master )print *, 'Run stopped after step no', istep - 1
c finish up
      t = etime( t2 )
      if( master )print *, 'Total cpu time was ', elapsd( t )
      call mpi_finalize( ierr )
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
      call disply
      call getmoi
      call gtintg
      call hedrec
      call mascmb
      call measure
      call ncheck
      call p2fndf
      call p3fndf
      call phyprp
      call ptbscl
      call rhoprf
      return
      end
