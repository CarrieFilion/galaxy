      program modefit
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c Main driving program for modefitting procedure.  It performs a non-linear
c   least-squares minimization to fit normal modes to data from N-body
c   simulations.  It attempts to estimate eigenfrequencies of overstabilities
c   (or instabilities or neutral modes), together the corresponding
c   eigenfunctions.  A number of distinct modes can be fitted simultaneously
c   although it is rare for the data to be of good enough quality for find
c   more than the dominant one or two.
c
c The program reads the data to be fitted from the RUNxxx.RES file.  The
c   following different data types can be fitted:
c      danl - sectorial harmonics of the mass distbn on grid rings
c      lgsp - logarithmic spiral coeffs
c      sfpc - coeffs saved from the SFP expansion
c      sphb - spherical Bessel functions coeffs
c      zanl - sectorial harmonics of the mid-plane displacements
c
c Any fit obtained can be saved in the direct access file ../MBASE (ie, in
c   a directory one level higher in your tree, which is handy if a different
c   sub-directory is used for each run).  A summary listing of all the fits
c   from a given run can be extracted from this file using program MODELIST
c   and best estimate eigenfrequencies with error bars can be obtained using
c   program ESTFREQ
c
c For some mathematical details of the method, see
c      /home/sellwood/docs/modefit.tex
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
      integer mdata
      parameter ( mdata = 50000 )
      integer ldata
      real apodc, data( 3, mdata )
      common / ldata / ldata, apodc, data
c
      integer wpar
      parameter ( wpar = 60000 )
      integer lwork, maxb
      real w( wpar )
      common / buffer / lwork, maxb, w
c
c local variables
      character type*4, yn*1
      integer ifail, m, nmodes
      logical first, found
      real sumsq
c store size parameters in the common blocks
      ldata = mdata
      lwork = wpar
c
c set defaults
      no = 6
      call setnag
      master = .true.
      call boilrp( .true. )
c read header record and set constants
      call header( .false. )
      call switch( 1 )
c initialise plotting
      call jsbgn
      first = .true.
      yn = 'y'
      do while ( yn .eq. 'y' )
        call datype( type )
c check modes file
        if( first )then
          call nextrec( type, ifail )
          call mbase( m, nmodes, sumsq, .true., found )
          first = .false.
        end if
c call fitting software
        call modfit( type )
c
        call gtchar( 'Try another data type (y/n)', yn )
      end do
      call jsend
      stop
      end

      real*8 function dfwght( E, Lz )
c function to weight DF when unequal particle masses are requested
c  value of this function is the inverse of the desired particle weights
      implicit none
c
c calling arguments
      real*8 E, Lz
c
c dummy version for analysis only
      dfwght = 1.
      return
      end
