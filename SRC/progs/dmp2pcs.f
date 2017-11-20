      program dmp2pcs
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, j, jrun
c
      call setnag
c read run number from standard input and open main ASCII I/O files
      call getset
      call boilrp( .true. )
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
c read the .dmp file
      jrun = irun
      call restrt
      if( irun .ne. jrun )then
        if( master )then
          print *, 'Files muddled'
          print *, '.dmp file from run', irun, ' but .dat file for run',
     +             jrun
        end if
        stop
      end if
      if( master )print *, 'Run', irun, ' restarted at step no', istep
      call chklst
c find gravitational field only if non-standard time stepping is used
      if( lfrv )then
        lprint = .true.
        phys = mod( istep, iphys ) .eq. 0
        if( bht )then
          bhmt = 5 * nbod
          bhnlev = 30
          allocate ( bhbox( bhnlev ) )
          allocate ( bhcom( 4, bhmt ) )
          allocate ( itup( bhmt ) )
          allocate ( itdown( bhmt / 8 ) )
          call bhtree
        else if( ldrct )then
c find total of all particles to be computed by dr3d
          ndrct = 0
          do i = 1, ncmp
            j = igrd( i )
            if( igrid( j ) .eq. 10 )ndrct = ndrct + nsp( i )
          end do
          allocate ( drpt( 15, ndrct ) )
          call masset
          call drctN2
        else
c assign mass to grids etc
          call masset
c determine gravitational field of a smooth (or the initial) mass distribution
          call findf( .true. )
        end if
      end if
c select output unit and create file
      ndistf = 10
      call unload
c save supplementary forces table if needed
      if( suppl )call suptab( .false. )
      end
