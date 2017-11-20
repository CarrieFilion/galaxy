      program analys
c  Copyright (C) 2015, Jerry Sellwood
c
c program for post-run analysis of the results from a single or multple
c  simulations
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
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/lunits.f'
c
c external
      logical gtlogl
c
c local array
      integer nopt
      parameter ( nopt = 24 )
      character*4 ab( nopt )
c
c local variable
      character aa*4
c
      data ab / 'lgan', 'prop', 'pnts', 'danl', 'epsi', 'cntd',
     +          'spct', 'vfld', 'ampl', 'lgav', 'rngs', 'intg',
     +          'moni', 'zprf', 'satl', 'kecp', 'frqs', 'Xplt',
     +          'mrpl', 'rhor', 'moit', 'slow', 'prjd', 'end ' /
c
c set defaults
      call setnag
      no = 6
      master = .true.
      call boilrp( .true. )
c read header record and set constants
      call header( .true. )
      call switch( 1 )
      if( gtlogl( 'Do you want physical scaling' ) )call set_scale
c initialise plot
      call jsbgn
      aa = '    '
      do while ( aa .ne. ab( nopt ) )
        print *, 'Choices are:'
        print '( 7( 6x, a4 ) )', ab
        call gtchar( 'Enter option', aa )
        if( aa .eq. ab( 1 ) )call loganl
        if( aa .eq. ab( 2 ) )call prpplt
        if( aa .eq. ab( 3 ) )call pntplt
        if( aa .eq. ab( 4 ) )call dnlplt
c        if( aa .eq. ab( 5 ) )call epsiln
c        if( aa .eq. ab( 6 ) )call contdn
        if( aa .eq. ab( 6 ) )call mkdens
        if( aa .eq. ab( 7 ) )call spctrm
        if( aa .eq. ab( 8 ) )call velfld
        if( aa .eq. ab( 9 ) )call amplot
        if( aa .eq. ab( 10 ) )call lgspav
        if( aa .eq. ab( 11 ) )call rngplt
        if( aa .eq. ab( 12 ) )call intplt
        if( aa .eq. ab( 13 ) )call monplt
        if( aa .eq. ab( 14 ) )call zprplt
        if( aa .eq. ab( 15 ) )call satanl
        if( aa .eq. ab( 16 ) )call kecmps
        if( aa .eq. ab( 17 ) )call frqplt
        if( aa .eq. ab( 18 ) )call xtplot
        if( aa .eq. ab( 19 ) )call mrplot
        if( aa .eq. ab( 20 ) )call rhoplt
        if( aa .eq. ab( 21 ) )call moiplt
        if( aa .eq. ab( 22 ) )call slowng
        if( aa .eq. ab( 23 ) )call prjplt
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
