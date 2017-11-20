      program plotpft
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c program to plot the radial and vertical variations of the potential
c   output from the program DFITER into a .pft file.  It also computes
c   the central attraction and plots the rotation curve
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
c program to read a .pft file and plot the potential it yields
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
      include 'inc/setup.f'
c
c allocatable arrays
      real*8, allocatable :: Bc1(:), Bc2(:), tr(:), tz(:)
c
c externals
      external mofr, potlf, vc, vdisc, vh
      real roundup
      real*8 potlf, vc
c
c local variables
      integer i, in, ip, j, jcmp, k, l, lr, lz, nrspl, nzspl
      real rat, ymin, ymax, xmin, xmax
c
c set defaults
      call setnag
      master = .true.
      call boilrp( .true. )
c read headers in the file
      call gtintg( 'Enter run number', irun )
      in = -1
      call opnfil( in, 'pft', 'unformatted', 'old', 'seq', i )
      if( i .ne. 0 )call crash( 'PLOTPFT', 'old .pft file not found' )
      read( in )nrspl, nzspl, ncmp, icmp, ncode
      do ip = 1, ncmp
        read( in )imtyp( ip ), impar( ip ), idftyp( ip ),
     +            ( dfcns( j, ip ), j = 1, 3 ), cmpmas( ip ),
     +            fmass( ip ), rscale( ip ), rtrunc( ip ),
     +            disc( ip ), ctype( ip )
      end do
      lr = nrspl - 4
      lz = nzspl - 4
      sphrod( icmp ) = lz .gt. 1
c read data
      allocate ( tr( nrspl ) )
      if( .not. sphrod( icmp ) )then
        allocate ( tz( nzspl ) )
        allocate ( Bc2( lr * lz ) )
        read( in, iostat = j )( tr( i ), i = 1, nrspl ),
     +                        ( tz( i ), i = 1, nzspl ),
     +                        ( Bc2( i ), i = 1, lr * lz )
      else
        allocate ( Bc1( nrspl ) )
        read( in, iostat = j )( tr( i ), i = 1, nrspl ),
     +                        ( Bc1( i ), i = 1, nrspl )
        if( j .ne. 0 )call crash( 'PLOTPFT', 'failed to read data' )
      end if
c data file is readable
      close( in )
      initls = .true.
c set frame scale
      xmin = 0
      xmax = .95 * tr( nrspl )
      rad = 0
c find central potential - rad is used just as the 2nd argument in phispl
      rad = 0
      iu = 1
      ymin = potlf( 0.d0 )
      ymin = -roundup( -ymin )
      ymax = 0
c plot potential function in two projections
      k = 1
      if( sphrod( icmp ) )k = 2
      call jsbgn
      do iu = 1, k
        call jspage
        call jscale( xmin, xmax, ymin, ymax )
        if( iu .eq. 1 )call jsaxis( 'x', 'R', 1 )
        if( iu .eq. 2 )call jsaxis( 'x', 'z', 1 )
        call jsaxis( 'y', '\\phi', 1 )
        l = 1
        if( sphrod( icmp ) )l = 7
        do i = 1, l
          rad = i - 1
          if( iu .eq. 1 )then
            rad = rad * rtrunc( icmp ) / real( l )
          end if
          call jsplt2( potlf )
        end do
      end do
c plot rotation curves
      ymin = 0
      rmax = .5 * rtrunc( icmp )
      ymax = 2. * vc( rmax )
      call jspage
      call jscale( xmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'R', 1 )
      call jsaxis( 'y', 'V_c', 1 )
      call jsplt2( vc )
      call jsdash( 2, 2, 2, 2 )
      call jsplt2( vh )
      call jsdash( 0, 1, 0, 1 )
      jcmp = icmp
      icmp = 1
      rat = cmpmas( jcmp ) / cmpmas( icmp )
      call jsplt2( vdisc )
      icmp = jcmp
      call jsbldt( 'disc dotted, halo dashed' )
      call jswrit( .05 * xmax, .95 * ymax )
      call jsbldt( 'halo/disk mass' )
      call jsbldf( rat, 4, 1 )
      if( impar( icmp ) .eq. 1 )then
        call jsbldt( 'Polytrope index=' )
        call jsbldf( dfcns( 1, icmp ), 3, 1 )
      else if( impar( icmp ) .eq. 2 )then
        call jsbldt( 'King model W_{0}=' )
        call jsbldf( dfcns( 1, icmp ), 3, 1 )
      end if
      call jswrit( .05 * xmax, .85 * ymax )
c integrated mass of (spherical) halo from central attraction in mid-plane
      call jspage
      call jsdash( 0, 0, 0, 0 )
      xmax = max( rtrunc( 1 ), rtrunc( 2 ) )
      ymax = 1.2 * cmpmas( 2 )
      call jscale( xmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'r', 1 )
      call jsaxis( 'y', 'm_{halo}(r)', 1 )
      call jsplt2( mofr )
c
      call jsend
      stop
      end

      real*8 function potlf( r )
      implicit none
c calling argument
      real*8 r
c
c common block
c
c rad is used just as the 2nd argument
      include 'inc/moment.f'
c
c external
      real*8 phispl
c
c radial variation at fixed vertical height
      if( iu .eq. 1 )then
        potlf = phispl( r, rad )
      else
c vertical variation at fixed radius
        potlf = phispl( rad, r )
      end if
      return
      end

      real*8 function vc( r )
      implicit none
c calling argument
      real*8 r
c
c common block
c
      include 'inc/moment.f'
c
c externals
      external potlf
      real*8 deriv2
c
      vc = 0
c compute circular speed from radial derivative of the mid-plane potential
      if( r .gt. 0. )then
        iu = 1
        rad = 0
        vc = deriv2( r, potlf )
        if( vc .gt. 0.d0 )then
          vc = sqrt( r * vc )
        else
          vc = 0
        end if
      end if
      return
      end

      real*8 function vh( r )
      implicit none
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c externals
      external potlf
      real*8 deriv2, vdisc, vhalo
c
      integer jcmp
c
      vh = 0
      if( r .gt. 0. )then
        iu = 1
        rad = 0
        vh = r * deriv2( r, potlf )
        jcmp = icmp
        do icmp = 1, ncmp
          if( icmp .ne. jcmp )then
            if( disc( icmp ) )then
              vh = vh - vdisc( r )**2
            else
              vh = vh - vhalo( r )**2
            end if
          end if
        end do
        icmp = jcmp
        vh = sqrt( abs( vh ) )
      end if
      return
      end

      real*8 function mofr( r )
      implicit none
c calling argument
      real*8 r
c
c externals
      real*8 vh
c
      mofr = r * vh( r )**2
      return
      end
