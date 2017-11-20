      program smooth
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Main driving program to make an optimal selection of particles from a
c   two-integral DF defined either for a razor-thin disk, or a sphere or
c   for an axisymmetric spheroid.  The idea is to obtain a set of particles
c   that represent the chosen DF as closely as possible - Poisson noise
c   associated with random particle selection is eliminated.
c
c The ouput from this program is a file (name.DFN) containing the requested
c   number of particles - by which is meant a set of 3 or 4 coordinate values.
c   3 coordinates (r, v_r & v_t) are created for disks & spheres, 4 coordinates
c   (r, v_r, v_t & z) are created for spheroids.  This file is in the format
c   expected for input to my program START.  The file can also be examined
c   using the program DFLOOK.
c
c The coordinates saved are the minimum number needed to define the integrals.
c   the remaining (ignorable) coordinates can be chosen trivially.  (For
c   spheroids, the v_r coordinate is the magnitude of the velocity in the
c   meridional plane - see section 3.2 of BT - which has to be resolved into
c   radial and vertical components.)
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      character*4 dftype, dtype, htype
      logical gtlogl
      real*8 distfn
c
c local variables
      character name*15
      integer i, j, k, n
      real a
c
c set defaults
      call setnag
      call gtintg( 'Enter 0 for terminal input', ni )
      if( ni .ne. 0 )call getset        
      master = .true.
      call boilrp( .true. )
c set up model
      call msetup
c check number of active DFs
      n = 0
      do j = 1, ncmp
        if( dist( j ) )n = n + 1
      end do
      if( n .eq. 0 )call crash( 'SMOOTH', 'No DF selected' )
      if( n .gt. 1 )then
        print *, n, ' DFs, which do you want?'
        do i = 1, ncmp
          if( dist( i ) )then
            if( cmprssd .and. ( .not. disc( i ) ) )then
              do j = 1, ncomp
                if( iccmp( j ) .eq. i )icomp = j
              end do
              j = ihalo0( icomp )
              k = idfn0( icomp )
              print '( '' compressed halo component'', i2, ' //
     +              '   '' has original halo type '', a4, ' //
     +           '   '' and DFtype '', a4 )', i, htype( j ), dftype( k )
            else
              j = imtyp( i )
              k = idftyp( i )
              if( disc( i ) )then
                print '( '' disk component'', i2, ' //
     +                '   '' has type '', a4, ' //
     +           '   '' and DFtype '', a4 )', i, dtype( j ), dftype( k )
              else
                print '( '' halo component'', i2, ' //
     +                '   '' has type '', a4, ' //
     +           '   '' and DFtype '', a4 )', i, htype( j ), dftype( k )
              end if
            end if
          end if
        end do
        j = 0
        do while ( j .lt. 1 .or. j .gt. ncmp )
          call gtintg( 'Select component number', j )
        end do
        do i = 1, ncmp
          if( i .ne. j )dist( i ) = .false.
        end do
      else
        do i = 1, ncmp
          if( dist( i ) )j = i
        end do
      end if
c set up numerical potential for a Shu DF in a tapered and/or thickened disk
      if( cdft( j ) .eq. 'SHUE' )then
        call smfield( j )
        close( ni )
c set flags for numerical evaluation
        numcal = .true.
        lgfld = .true.
        dist( j ) = .true.
        call cutoff( rtrunc( j ) )
c initialize the DF
        icmp = j
        a = distfn( Emine, 0. )
      end if
c number of particles
      do icmp = 1, ncmp
        if( dist( icmp ) )then
          call gtintgs( 'Enter ne, nh and nhe', jdfcs( 1, icmp ), 3 )
          jdfcs( 1, icmp ) = max( 1, jdfcs( 1, icmp ) )
          jdfcs( 2, icmp ) = max( 1, jdfcs( 2, icmp ) )
          jdfcs( 3, icmp ) = max( 1, jdfcs( 3, icmp ) )
          n = 1
          do j = 1, 3
            n = n * jdfcs( j, icmp )
          end do
          print *, n, ' particles will be generated'
          uqmass = .not. gtlogl( 'Equal mass particles' )
c name and open .dfn file
          name( 1:4 ) = cdft( icmp )
          call lowercase( name( 1:4 ) )
          n = n / 1000
          a = 0
          if( n .gt. 0 )a = log10( real( n ) ) + .001
          j = int( a ) + 1
          write( name( 5:15 ), '( i11 )' )n
          name( 5:j+4 ) = name( 16-j:15 )
          j = j + 5
          name( j:j+4 ) = 'k.dfn'
          ndistf = 1
          open( ndistf, file = name( 1:j+4 ), status = 'unknown',
     +          form = 'unformatted' )
c select random seed
          if( lnag )then
            call gtintg(
     +    'Enter random seed, or 0 for default seed', jdfcs( 4, icmp ) )
            if( jdfcs( 4, icmp ) .eq. 0 )then
              print *, 'Using default seed'
            else
              print *, 'Random seed is', jdfcs( 4, icmp )
              call rvseed( jdfcs( 4, icmp ) )
            end if
          end if
c choose particles
          call dfset
        end if
      end do
      stop
      end

c      include 'dfwght.f'
