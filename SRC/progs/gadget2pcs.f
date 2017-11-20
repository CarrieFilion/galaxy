      program gad2pcs
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c program to convert a GADGET-2 particles file to a .pcs file from which
c   GALAXY can be run
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
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c externals
      integer lnblnk
      logical gtlogl
c
c input buffers
      integer, allocatable :: id(:)
      real, allocatable :: mass(:)
      real, allocatable :: outb(:)
      real, allocatable :: v(:,:)
      real, allocatable :: x(:,:)
c
c gadget arrays
      integer npart( 6 ), nall( 6 ), nallhw( 6 ), dummy( 15 )
      real*8 massar( 6 )
c gadget variables
      integer flagage, flagcooling, flagfeedback, flagmetals, flagsfr,
     +        flag_entr_ics, numfiles
      real*8 boxsize, Hubbleparameter, Omega0, OmegaLambda, redshift,
     +       time
c
c local variables
      character fname*80, tflag*10, form*6
      logical lact( 6 ), lm
      integer i, j, k, l, lo, m, n
      real pm, t, rm( 6 ), vm( 6 ), zm( 6 )
      parameter ( l = 5000 )
      data nsp / mcmp * 0 /
      pertbn = .false.
c
      call setnag
      master = .true.
      call boilrp( .false. )
c open GADGET-2 file and read header
      call gtchar( 'Enter name of GADGET-2 file', fname )
      j = lnblnk( fname )
      open( 1, file = fname( 1:j ),
     +      form = 'unformatted', status = 'old' )
      read( 1 )npart, massar, time, redshift, flagsfr, flagfeedback,
     +         nall, flagcooling, numfiles, boxsize, omega0,
     +         omegalambda, hubbleparameter, flagage, flagmetals,
     +         nallhw, flag_entr_ics, dummy
c
      print 100, ' npart:', npart
  100 format( a, 6i12 )
      print 101, 'massar:', massar
  101 format( a, 6f12.8 )
c      print 101, 'time, redshift', time, redshift
c      print 100, 'flagsfr, flagfeedback', flagsfr, flagfeedback
c      print 100, '  nall:', nall
c      print 100, 'flagcooling, numfiles', flagcooling, numfiles
c      print 101, 'boxsize, omega0, omegalambda, hubbleparameter',
c     +          boxsize, omega0, omegalambda, hubbleparameter
c      print 100, 'flagage, flagmetals', flagage, flagmetals
c      print 100, 'nallhw:', nallhw
c      print 100, 'flag_entr_ics', flag_entr_ics
c determine array sizes
      lm = .false.
      nbod = 0
      j = 0
      do k = 1, 6
        nbod = nbod + npart( k )
        lact( k ) = npart( k ) .gt. 0
        if( lact( k ) )then
          j = j + 1
          if( j .le. 3 )then
            nsp( j ) = npart( k )
          else
            call crash( 'GAD2PCS', 'Too many populations of particles' )
          end if
        end if
        if( lact( k ) )lm = lm .or. ( massar( k ) .eq. 0.d0 )
      end do
      if( lact( 1 ) .or. lact( 5 ) .or. lact( 6 ) )call crash(
     +       'GAD2PCS', 'No gas, star or bndary particles are allowed' )
      allocate ( id( nbod ) )
      allocate ( x( 3, nbod ) )
      allocate ( v( 3, nbod ) )
      if( lm )allocate ( mass( nbod ) )
c read positions
      read( 1 )x
c read velocities
      read( 1 )v
c read ids
      read( 1 )id
c read masses if present
      if( lm )read( 1 )mass
      close( 1 )
c
c scaling
c
      print *, 'This program assumes that the input file gives' //
     +         ' distances in kpc,'
      print *, 'velocities in km/s, and masses in units of 10^10' //
     +         ' solar masses'
      if( .not. gtlogl( 'Is this true' ) )call crash(
     +             'GAD2PCS', 'Non-standard scaling requires recoding' )
c GADGET-2 manual recommends
c   UnitVelocity in cm per s  = 1e5 = 1 km/s
c   UnitLength in cm  = 3.085678e21 = 1 kpc
c   UnitMass in g        = 1.989e43 = 1e10 solar masses
c scale to GALAXY units assuming
c                       length unit = 3 kpc
c                     velocity unit = 293 km/s
c                         mass unit = 6e10 solar masses
c                        (time unit = 10 Myr)
      do j = 1, nbod
        do i = 1, 3
          x( i, j ) = x( i, j ) / 3.
          v( i, j ) = v( i, j ) / 293.
        end do
        if( lm )mass( j ) = mass( j ) / 6.
      end do
      if( .not. lm )then
        do k = 2, 4
          massar( k ) = massar( k ) / 6.
        end do
      end if
c
c end scaling
c
c determine whether masses of all particles are equal
      uqmass = lm
      if( .not. lm )then
        do k = 2, 3
          if( lact( k ) )then
            do j = k + 1, 4
              if( lact( j ) )uqmass = uqmass .or.
     +                                  ( massar( j ) .ne. massar( k ) )
            end do
          end if
        end do
      end if
      print *, 'uqmass set to', uqmass
c open output file
      call gtchar( 'Enter output file name', fname )
      j = lnblnk( fname )
      open( 10, file = fname( 1:j ),
     +      form = 'unformatted', status = 'new', iostat = i )
      if( i .ne. 0 )then
        print *,'A file with this name already exists'
        if( gtlogl( 'Do you want to overwrite it' ) ) then
           open( 10, file = fname( 1:j ),
     +           form = 'unformatted', status = 'old' )
        else
          call crash( 'GAD2PCS', 'Declined to overwite output file' )
        end if
      end if
c write a short header
      t = time
      nwpp = 6
      if( uqmass )nwpp = nwpp + 1
      pm = 1
      print *, 'Creating output file at time', t
      write( 10 )( nsp( i ), i = 1, 3 ), nwpp, l, t, pm, pertbn
c allocate an output buffer
      lo = l * nwpp
      allocate ( outb( lo ) )
c set counters
      j = 0
      n = 0
c work over active populations
      do k = 2, 4
        if( lact( k ) )then
          do m = n + 1, n + npart( k )
c output the buffer when full
            if( j .ge. lo )then
              write( 10 )outb
              j = 0
            end if
c copy the coordinates
            do i = 1, 3
              outb( j + i ) = x( i, m )
              outb( j + i + 3 ) = v( i, m )
            end do
            if( uqmass )then
              if( lm )then
                outb( j + 7 ) = mass( m )
              else
                outb( j + 7 ) = massar( k )
              end if
            end if
c increment buffer pointer
            j = j + nwpp
          end do
          n = n + npart( k )
        end if
      end do
c output any remaining particles
      if( j .gt. 0 )then
        write( 10 )( outb( i ), i = 1, j )
      end if
      close( 10 )
c report success
      j = lnblnk( fname )
      print *, 'File '// fname( 1:j ) //' created'
      n = 1
      if( t .gt. 1. )n = 1. + log10( t )
      if( n .le. 6 )then
        form = '( i1 )'
        write( form( 4:4 ), '( i1 )' )n
        tflag = '.pcs'
        write( tflag( 5:10 ), form )nint( t )
        k = lnblnk( tflag )
        if( fname( j-k+1:j ) .ne. tflag( 1:k ) )print *,
     +                    'It should be renamed runX' // tflag( 1:k ) //
     +                                             ' for use by pcs2dmp'
      end if
      end
