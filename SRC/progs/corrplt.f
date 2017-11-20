      program corrplt
c  Copyright (C) 2015, Jerry Sellwood
c
c program to plot the azimuthally averaged central attraction of the particles
c   (red dotted curve) and a comparison curve (black solid) from the predicted
c   circular speed.  The difference between the two curves is the supplemental
c   central attraction added to the forces acting on every particle at every
c   step in the subsequent simulation.
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
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      integer lnblnk
c
c local array
      integer m
      parameter ( m = 1000 )
      real data( 6, m )
c
c local variables
      character*120 line
      integer i, in, j, k, l, n
      logical gtlogl
      real fmin, fmax, rm, x, y
c
c set defaults
      call setnag
      master = .true.
      call boilrp( .false. )
c get run number
      l = 0
      call gtintg( 'Enter run number', irun )
      if( irun .gt. 0 )then
c open .lis file
        in = -1
        call opnfil( in, 'lis', 'formatted', 'old', 'seq', i )
c skip down to header created by supset
        do while ( line( 1:10 ) .ne. '    radius' )
          read( in, '( a )', end = 1 )line
        end do
        k = 5
        if( gtlogl( 'Is it a disk' ) )k = 2
c read data
        i = 0
        n = 0
        fmin = 1.e10
        fmax = 0.001
        do while ( i .eq. 0 )
          read( in, '( a )', end = 2 )line
          n = n + 1
          if( n .gt. m )call crash( 'MAIN', 'Too many lines' )
          read( line, 
     + '( 3x, f6.3, 5f12.6 )', iostat = i )( data( j, n ), j = 1, 6 )
          fmin = min( fmin, data( k, n ) )
          fmax = max( fmax, data( k, n ) )
        end do
        n = n - 2
    2   close( in )
        print *, n, ' lines read'
        fmin = 1.1 * fmin
        fmax = 1.1 * fmax
        fmax = fmax + .02 * ( fmax - fmin )
c set up plot
        l = l + 1
        if( l .eq. 1 )then
          call jsbgn
          call jspage
          rm = 1.1 * data( 1, n )
          call jscale( 0., rm, fmin, fmax )
          call jsaxis( 'x', 'R', 1 )
          call jsaxis( 'y', 'a_r', 1 )
          do i = 1, n
            x = data( 1, i )
            y = data( k, i )
            if( i .eq. 1 )call jsmove( x, y )
            call jsline( x, y )
          end do
        end if
        call pgsci( l + 1 )
        call jsdash( 0, 1, 0, 1 )
        do i = 1, n
          x = data( 1, i )
          y = data( 3, i )
          if( i .eq. 1 )call jsmove( x, y )
          call jsline( x, y )
        end do
      end if
      call jsend
      stop
    1 print *, 'end of file encountered without finding the keyword'
      print *, 'last line read was'
      i = lnblnk( line )
      print '( a )', line( 1:i )
      end

      include 'inc/pgiden.f'
