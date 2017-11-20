      program coulomb_test
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c program to assign a few unit masses to the grid, solve for the field,
c   and check the result against a solution by direct convolution
c Requires mpi, but can also be used for a single cpu
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
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'mpif.h'
c
c local variables
      character*80 line
      integer i, n
      real t1, t2, ar( 2 ), etime
c
      call setnag
c initialize
      call getset
      call boilrp( .false. )
      if( numprocs .gt. mnodes )call crash( 'MAIN', 'mnodes too small' )
      rewind no
c set parallel flag
      parallel = numprocs .gt. 1
c find grid keyword in input .dat file
      if( master )then
        rewind ni
        read( ni, '( a )' )line
        n = 1
        do while ( line( 1:9 ) .ne. 'grid size' )
          read( ni, '( a )' )line
          call lowercase( line( 1:9 ) )
          n = n + 1
        end do
c reposition file for grdset
        rewind ni
        do i = 1, n - 1
          read( ni, '( a )' )line
        end do
      end if
c read grid parameters and create Green function
      call grdset
      call greenm( .true. )
c      call gcheck
c      read *, istep
c set number of zones
      call getline( ni, line )
      if( line( 1:6 ) .ne. 'NZONES' )call crash( 'PTEST',
     +                                           'nzones card missing' )
      read( line( 11:40 ), * )nzones
      if( master )then
        print *, 'nzones =', nzones
        write( no, '( '' nzones ='', i3 )' )nzones
      end if
c determine space needed
      call setptr
c allocate space
      call setspc
c position test masses
      call chckst
      if( parallel )call mpi_barrier( mpi_comm_world, i )
c evaluate the field
      iphys = 10
      phys = .true.
      t1 = etime( ar )
      call findf( phys )
      t2 = etime( ar )
      if( master )print *, 'cpu time was ', t2 - t1, ' sec', t2, t1
      call checkf
      if( parallel )call mpi_finalize( i )
      stop
      end

      subroutine polcat
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      if( master )print *, 'Skipping conversion to Cartesian components'
      return
      end

      subroutine fiddle
      call gtintg
      call mascmb
      call p2fndf
      call p3fndf
      end
