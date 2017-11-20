      program mkpcs
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c program to create a .pcs file for use in the simulation code GALAXY
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
c This program creates the initial positions & velocities of the particles
c   in the appropriate format and output a file ready to be read by BEGIN.
c
c It does not make sense to run this in parallel on multiple processors
c
c It reads data from a short input file (runXXX.dat) which sets the grid
c   size, softening length, number of particles, time step size, etc
c   and also determines which quantities are to be measured as the model
c   evolves.  These parameters are set in subroutines GRDSET and SETRUN
c The main files associated with the run are:
c   runXXX.dat - ASCII input: grid parameters etc & analysis instructions
c   runXXX.lis - ASCII output: a brief summary of progress
c   runXXX.pcs - binary output: the product that will be read by BEGIN 
c   runXXX.dfn - optional binary input: output from smooth
c   runXXX.grt - binary: a large Green's fn table needed for some force
c                determination methods
c
c The initial positions of the particles are set in subroutine PSETUP
c A header record is written to the .res file by HEDREC
c EITHER the masses of the particles are assigned to the grid using MASSET
c     OR the smooth analytic mass distbn is assigned to the grid using SMMASS
c The gravitational field is determined by a call to FINDF
c Any remaining velocity components not yet initialized are set by ORBSET
c Particles might be centered on the origin and net momemtum removed
c Different mass components might be shifted spatially and given bulk velocities
c Finally the .pcs file is created
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, j
      logical smooth
c
      call setnag
c read run number from standard input and open main ASCII I/O files
      call getset
c rewind .lis file - tidies up if there was a previous false start
      if( master )rewind no
      call boilrp( .true. )
c read .dat file (ASCII input)
      call msetup
      call grdset
      call setrun
      close( ni )
c we get a better setup if we use a smooth density rather than the noisy mass
c   arrays from the particles, but that part of the code needs NAG and has
c   problems with multiple active components
      j = 0
      do i = 1, ncmp
        if( ( .not. testp( i ) ) .and. ( .not. rigidp( i ) ) )j = j + 1
      end do
      smooth = lnag .and. ( j .eq. 1 )
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
c place particles in initial positions and (possibly) set some vel components
      call psetup
c assign mass to grids etc
      phys = .false.
      lprint = .true.
      if( smooth .and. lgrd )then
        call smmass
      else
        call masset
      end if
c determine gravitational field of a smooth (or the initial) mass distribution
      call findf( .true. )
c determine any supplementary radial forces
      if( suppl )call supset
c set initial velocities of disk particles only
      call orbset
c remove net offset in CoM and net momentum of the initial model if requested
      call arrest
c add any requested bulk shifts and/or motion of separate components
      call xvshft
c save initial set up
      call unload
c save supplementary forces table, if it exists
      if( suppl )call suptab( .false. )
      stop
      end
