      subroutine decode( type, itype, nw )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine from the analysis software that interprets the data type in
c   the supplied character string and calculates the expected length of
c   data records of that type.
c The itype value is used to determine the logical unit nummber if the
c   the different data types have been split into separate files
c
c calling arguments
      character*4 type
      integer itype, nw
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
c external
      character*4 datnam
c
c local variables
      character*1 yn
      integer i
c
c initialise resfil
      if( ntype .eq. 0 )then
        byswap = .false.
        yn = datnam( 1 )
      end if
c decode name
      nw = 0
    1 itype = 1
      do while ( type .ne. datnam( itype ) )
        itype = itype + 1
        if( itype .gt. ntype )then
c check for and skip an extra header record
          itype = 0
          read( type, '( a4 )' )i
          if( p2d .and. ( i .eq. na ) )return
          if( c3d .and. ( i .eq. ngy ) )return
c unidentified string
          print '( 10x, ''Unrecognised key word'', 2x, a10 )', type
          if( .not. byswap )then
            call gtchar( 'Do you want to try byte swapping (y/n)', yn )
            if( yn .eq. 'y' )then
              byswap = .true.
              call swap( type )
              go to 1
            end if
          end if
          call crash( 'DECODE', 'Unidentified key word' )
        end if
      end do
c
c set number of words for this record type
c
      nw = -1
c picture data
      if( itype .eq. 1 )nw = mr
c density on polar grid
      if( itype .eq. 2 )nw = mr * ma
c potential distribution on polar grid
      if( itype .eq. 3 )nw = mr * ma
c frequencies
      if( itype .eq. 4 )nw = ma * mr
c logarithmic spiral coefficients
      if( itype .eq. 5 )nw = 2 * mr * ma
c angular momentum distribution
      if( itype .eq. 6 )then
        if( mr .eq. ma )then
          nw = mr
         else
          nw = mr * ma
         end if
      end if
c disc particle velocity field
      if( itype .eq. 7 )nw = mr * ma * nprop
c integrals
      if( itype .eq. 8 )nw = mr
c orbit data
      if( itype .eq. 9 )nw = mr
c grey scale data
      if( itype .eq. 10 )nw = mr * ma
c Fourier transformed density distribution
      if( itype .eq. 11 )nw = 2 * mr * ma
c integral monitoring data
      if( itype .eq. 12 )nw = mr * ma
c SFP coefficients
      if( itype .eq. 13 )nw = 2 * ma + 2
c z mean analysis coefficients
      if( itype .eq. 14 )nw = 2 * mr * ( ma + 1 )
c halo particle velocity field
      if( itype .eq. 15 )nw = mr * ma * nprop
c spherical Bessel coefficients
      if( itype .eq. 16 )nw = ma * ( mr + 1 ) * ( mr + 2 )
c z-density profile data
      if( itype .eq. 17 )nw = mr * ma
c particle rings
      if( itype .eq. 18 )nw = ( ncoor / 2 ) * mr * ma
c scf coeffs for 3-D
      if( itype .eq. 19 )nw = ma + 1
c radial volume density profile data
      if( itype .eq. 20 )nw = mr
c satellite data
      if( itype .eq. 21 )nw = mr
c logarithmic spiral coefficients of velocity
      if( itype .eq. 22 )nw = 4 * mr * ma
c moment of inertia tensors
      if( itype .eq. 23 )nw = 10 * mr
c radial surface density profile data
      if( itype .eq. 24 )nw = mr
c coeffs from s3d grid
      if( itype .eq. 25 )nw = 2 * mr * ( ma + 1 ) * ( ma + 2 )
c grid centroids
      if( itype .eq. 26 )nw = mr * ma
c mode data
c      if( itype .eq. 22 )nw = 2 * mr
c gas angular momentum distribution
c      if( itype .eq. 23 )nw = mr
c gas particle velocity field
c      if( itype .eq. 24 )nw = mr * ma * nprop
c Fourier transformed potential distribution
c      if( itype .eq. 25 )nw = 2 * mr * ma
c "critical" spiral analysis
c      if( itype .eq. 26 )nw = 2 * mr * ma
c particle angular momenta
      if( itype .eq. 27 )nw = mr
c projected denstty
      if( itype .eq. 28 )nw = ma
c
c check that nw was set
      if( nw .lt. 0 )call crash( 'DECODE', 'Unrecognized type' )
      return
      end
