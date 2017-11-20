      program genplt
c  Copyright (C) 2014, Jerry Sellwood
c
c program to plot physical properties of a selected composite or single
c   component mass model.  Additional plots are made if a DF is available
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
      implicit none
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
      include 'inc/orbanl.f'
c
      logical dblint, df, disk, sphd
      common / genloc / dblint, df, disk, sphd
c
c externals
      logical gtlogl
      real*8 gmassd, gmassh, gsigma, gsigmi, phitot, xtoomn
      external frtot, gmassd, gmassh, gsigma, lndbld, phitot
      external qtoom, lrhohal, lrhohli, sigmau, sigmav, vcirc, vmean
      external v2mean, gsigmi, qtoomi, sigmui, vmeani, v2meani, vdisc
      external vhalo, lgsigma, lgsigmi, akappa, xtoomn, kepl, vtot
      external frtotl, rhohal, rhohli
c
c local variables
      integer jcmp
      real a, b, rmin, xmax, ymax, ymin
c
c preliminaries
      lnag = .true.
      master = .true.
      call gtintg( 'Enter 0 for terminal input', ni )
      if( ni .ne. 0 )then
        call getset
        rewind no
      end if
c input selected model
      call msetup
      if( ncmp .gt. 1 )then
        print *, 'Multi-component model', ncmp
        call gtintg( 'Select which component to plot', icmp )
      else
        icmp = 1
      end if
c decide type of plot
      df = dist( icmp )
      disk = disc( icmp )
      sphd = sphrod( icmp )
      if( df )dblint = 
     +       gtlogl( 'Curves drawn from double integrals over distfn?' )
c set radial extent
      xmax = rmax
      if( dblint )then
        call tapset
      else
        xmax = 1.5 * xmax
        call gtreal( 'Enter new rmax', xmax )
        call cutoff( xmax )
      end if
c run some sanity checks
      call goofch
c
      rmin = 0
      if( ctype( 1 ) .eq. 'JAFF' .or. ctype( 1 ) .eq. 'NFW ' )rmin = .01
      xmax = rmax
c initialize plot
      call jsbgn
c surface density
      if( disk )then
        call newplot
        call yrange( rmin, xmax, ymin, ymax, lgsigma )
        call jscale( rmin, xmax, ymin, ymax )
        call jsaxis( 'x', 'radius', 1 )
        call jsaxis( 'y', 'ln\\Sigma', 1 )
        if( disk )then
          if( dblint )then
            call jsplt2( lgsigmi )
            call jsdash( 2, 2, 2, 2 )
            call jsplt2( lgsigma )
            call jsdash( 0, 0, 0, 0 )
          else
            call jsplt2( lgsigma )
          end if
        end if
      end if
c density profiles
      if( .not. sphd )then
        call newplot
        ymin = 0
        if( disc( icmp ) )then
          call yrange( rmin, xmax, ymin, ymax, gmassd )
        else
          call yrange( rmin, xmax, ymin, ymax, gmassh )
        end if
        call jscale( rmin, xmax, ymin, ymax )
        call jsaxis( 'x', 'radius', 1 )
        call jsaxis( 'y', 'm(r)', 1 )
        if( disc( icmp ) )then
          call jsplt2( gmassd )
        else
          call jsplt2( gmassh )
        end if
c
        if( .not. disc( icmp ) )then
          call newplot
c linear scaling of density
          call yrange( rmin, xmax, ymin, ymax, rhohal )
          call jscale( rmin, xmax, ymin, ymax )
          call jsaxis( 'x', 'radius', 1 )
          call jsaxis( 'y', '\\rho', 1 )
          call jsplt2( rhohal )
          if( df )then
            call pgsci( 2 )
            call jsdash( 2, 2, 2, 2 )
            call jsplt2( rhohli )
            call jsdash( 0, 0, 0, 0 )
            call pgsci( 1 )
          end if
c lograithmic scaling of density
          call newplot
          a = -3
          b = log10( xmax )
          ymin = -10
          call yrange( a, b, ymin, ymax, lrhohal )
c          ymin = max( ymin, ymax - 5. )
          call jscale( a, b, ymin, ymax )
          call jslaxs( 'x', 'radius', 1 )
          call jslaxs( 'y', '\\rho', 1 )
          call jsplt2( lrhohal )
          if( df )then
            call jsdash( 2, 2, 2, 2 )
            call pgsci( 2 )
            call jsplt2( lrhohli )
            call jsdash( 0, 0, 0, 0 )
            call pgsci( 1 )
          end if
        end if
        if( .not. disc( icmp ) )then
          call newplot
          a = -6
          b = log10( xmax )
          ymin = -10
          call yrange( a, b, ymin, ymax, frtotl )
          ymin = max( ymin, ymax - 5. )
          call jscale( a, b, ymin, ymax )
          call jslaxs( 'x', 'radius', 1 )
          call jslaxs( 'y', 'radial force', 1 )
          call jsplt2( frtotl )
        end if
      end if
c total potential
      call newplot
      call yrange( rmin, xmax, ymin, ymax, phitot )
      call jscale( rmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'radius', 1 )
      call jsaxis( 'y', '\\Phi', 1 )
      call jsplt2( phitot )
c radial attraction
      call newplot
      call yrange( rmin, xmax, ymin, ymax, frtot )
      call jscale( rmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'radius', 1 )
      call jsaxis( 'y', 'radial force', 1 )
      call jsplt2( frtot )
c circular speed
      call newplot
      call yrange( rmin, xmax, ymin, ymax, vcirc )
      call jscale( rmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'radius', 1 )
      call jsaxis( 'y', 'V_c', 1 )
      call jsplt2( vcirc )
      call jsdash( 2, 2, 2, 2 )
      call jsplot( kepl )
      call jsdash( 0, 0, 0, 0 )
c mean circular spped - if possible
      if( df )then
        if( gtlogl( 'Plot vmean?' ) )then
          call jsdash( 2, 2, 2, 2 )
          if( dblint )then
            call jsplt2( vmeani )
          else
            call jsplt2( vmean )
          end if
          call jsdash( 0, 0, 0, 0 )
        end if
      end if
c circular speeds of separate components
      if( ncmp .gt. 1 )then
        call jsdash( 2, 2, 2, 2 )
        jcmp = icmp
        do icmp = 1, ncmp
          if( disc( icmp ) )then
            call jsplt2( vdisc )
          else
            call jsplt2( vhalo )
          end if
        end do
        call jsdash( 0, 1, 0, 1 )
        call jsplt2( vtot )
        icmp = jcmp
      else
        call jsdash( 0, 1, 0, 1 )
        call jsplt2( vtot )
      end if
c radial velocity dispersion - if possible
      if( df )then
        if( gtlogl( 'Plot sigmau?' ) )then
          call newplot
          call yrange( rmin, xmax, ymin, ymax, sigmau )
          ymin = 0
          call jscale( rmin, xmax, ymin, ymax )
          call jsaxis( 'x', 'radius', 1 )
          call jsaxis( 'y', '\\sigma_u', 1 )
          call jsplt2( sigmau )
          call jsdash( 2, 2, 2, 2 )
          if( dblint )call jsplt2( sigmui )
          call jsdash( 0, 0, 0, 0 )
        end if
      end if
c principal frequencies
      call newplot
      morb = 1
      lorb = 0
      call yrange( rmin, xmax, ymin, ymax, lndbld )
      ymax = min( ymax, 10. )
      call jscale( rmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'radius', 1 )
      call jsaxis( 'y', 'frequency', 1 )
      call jsplt2( lndbld )
      call jsplt2( akappa )
      do morb = 2, 3
        call jscale( rmin, xmax, ymin, real( morb ) * ymax )
        lorb = -1
        if( morb .eq. 2 )call jsdash( 2, 2, 2, 2 )
        if( morb .eq. 3 )call jsdash( 0, 1, 0, 1 )
        call jsplt2( lndbld )
        lorb = 1
        call jsplt2( lndbld )
      end do
      call jsdash( 0, 0, 0, 0 )
c Toomre's X parameter
      if( disk )then
        call newplot
        call yrange( rmin, xmax, ymin, ymax, xtoomn )
        call jscale( rmin, xmax, ymin, ymax )
        call jsaxis( 'x', 'radius', 1 )
        call jsaxis( 'y', 'X * m', 1 )
        call jsplt2( xtoomn )
      end if
c Toomre's Q parameter
      if( df .and. disk )then
        call newplot
        if( dblint )then
          call yrange( rmin, xmax, ymin, ymax, qtoomi )
        else
          call yrange( rmin, xmax, ymin, ymax, qtoom )
        end if
        call jscale( rmin, xmax, ymin, ymax )
        call jsaxis( 'x', 'radius', 1 )
        call jsaxis( 'y', 'Toomre Q', 1 )
        if( dblint )then
          call jsplt2( qtoomi )
        else
          call jsplt2( qtoom )
        end if
      end if
c azimuthal velocity dispersion
      if( df )then
        if( gtlogl( 'Plot v2mean?' ) )then
          call newplot
          if( dblint )then
            call yrange( rmin, xmax, ymin, ymax, v2meani )
          else
            call yrange( rmin, xmax, ymin, ymax, v2mean )
          end if
          call jscale( rmin, xmax, ymin, ymax )
          call jsaxis( 'x', 'radius', 1 )
          call jsaxis( 'y', '<V^{2}>', 1 )
          if( dblint )then
            call jsplt2( v2meani )
          else
            call jsplt2( v2mean )
          end if
        end if
      end if
c
      call jsend
      stop
      end

      real*8 function xtoomn( r )
c evaluates Toomre's X parameter at the give radius
      implicit none
c
c calling argument
      real*8 r
c
c externals
      real*8 gsigma, xtoom
c
      xtoomn = 0
      if( gsigma( r ) .gt. 0.d0 )xtoomn = xtoom( r, 1 )
      return
      end

      subroutine goofch
c checks that the numerical derivatives of some functions match
c   the programmed gradients
      implicit none
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      logical dblint, df, disk, sphd
      common / genloc / dblint, df, disk, sphd
c
c externals
      external gmassd, gmassh, phitot, vcirc
      real*8 deriv2, frtot, gmassh, gsigma, gsigmi
      real*8 rhohal, rhohli, vcirc, vcgrad, vdisc, vhalo
c
c local variables
      integer jcmp
      real a, b
      real*8 r
      include 'inc/pi.f'
c
      r = .5
      do while ( r .gt. 0. )
        print *, 'Checking values for r =', sngl( r )
        if( disk )then
c surface density checks
          a = 2. * pi * r * gsigma( r )
          b = deriv2( r, gmassd )
          print *, '   2pi Sigma & dM/dr', a, b, a - b
          if( dist( icmp ) )then
            a = gsigma( r )
            b = gsigmi( r )
            print *, '       Sigma & Sigmi', a, b, a - b
          end if
        end if
c radial force and potential derivative
        a = frtot( r )
        b = -deriv2( r, phitot )
        print *, '    frtot & -dPhi/dr', a, b, a - b
c circ vel and radial force
        b = frtot( r )
        if( b .lt. 0. )then
          a = vcirc( r )
          b = sqrt( -r * b )
          print *, '  vc & sqrt( -r*fr )', a, b, a - b
          if( ncmp .gt. 1 )then
            b = 0
            jcmp = icmp
            do icmp = 1, ncmp
              if( disc( icmp ) )then
                b = b + vdisc( r )**2
              else
                b = b + vhalo( r )**2
              end if
            end do
            b = sqrt( b )
            print *, 'vc & sqrt( vd*vd + vh*vh )', a, b, a - b
            icmp = jcmp
          end if
c circ vel gradient
          a = vcgrad( r )
          b = deriv2( r, vcirc )
          print *, '    vcgrad & dv_c/dr', a, b, a - b
        end if
        if( ( .not. disk ) .and. ( .not. sphd ) )then
c radial force and m( r )
          if( ncmp .eq. 1 )then
            a = -r**2 * frtot( r )
            b = gmassh( r )
            print *, '  -r^2 frtot & GM(r)', a, b, a - b
            if( df )then
              a = rhohal( r )
              b = rhohli( r )
              if( b .gt. 0. )then
                print *, ' rho(r) & int f d^3v', a, b, a / b
              else
                print *, 'density from integrating DF =', b,
     +                   ' but should be', a
              end if
            end if
          end if
c halo density
          a = 4. * pi * r**2 * rhohal( r )
          b = deriv2( r, gmassh )
          print *, ' 4pi*r^2*rho & dM/dr', a, b, a - b
        end if
        call gtdble( 'Enter r', r )
      end do
      return
      end

      subroutine newplot
c sets up and labels a new panel
      implicit none
c
c common blocks
c
      include 'inc/jscmmn.f'
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      logical dblint, df, disk, sphd
      common / genloc / dblint, df, disk, sphd
c
c externals
      character*4 dftype, dtype, htype
c
c local variable
      integer jcmp
c
      jcmp = icmp
      if( npic .gt. 2 )call jschsz( .2 )
      call jspage
      call jssize( .18, .98, .18, .98 )
      if( ipic .gt. 1 )return
      call jscale( 0., 1., 0., 1. )
      do icmp = 1, ncmp
        if( disc( icmp ) )then
          call jsbldt( 'Disc type =' )
          call jsbldt( dtype( imtyp( icmp ) ) )
          if( ctype( icmp ) .eq. 'SOFT' )then
            call jsbldt( ' r_* =' )
            call jsbldf( sngl( rstar ), 5, 2 )
          end if
          if( cmpmas( icmp ) .lt. 1.d0 )then
            call jsbldt( ' active disc mass =' )
            call jsbldf( cmpmas( icmp ), 5, 2 )
          end if
        else
          call jsbldt( ' Halo type =' )
          call jsbldt( htype( imtyp( icmp ) ) )
          if( rscale( icmp ) .ne. 1.d0 )then
            call jsbldt( 'halo scale =' )
            call jsbldf( rscale( icmp ), 5, 2 )
          end if
          if( dfcns( 3, icmp ) .ne. 0.d0 )then
            call jsbldt( 'halo parameter =' )
            call jsbldf( dfcns( 3, icmp ), 5, 2 )
          end if
        end if
        if( dist( icmp ) )then
          call jsbldt( ' DF type =' )
          call jsbldt( dftype( idftyp( icmp ) ) )
          if( dfcns( 2, icmp ) .eq. 0.d0 )then
            call jsbldt( 'DF parameter =' )
            call jsbldf( dfcns( 1, icmp ), 5, 2 )
          else
            call jsbldt( 'DF parameters =' )
            call jsbldf( dfcns( 1, icmp ), 5, 2 )
            call jsbldf( dfcns( 2, icmp ), 5, 2 )
          end if
        end if
      end do
      call jswrit( .1, 1.5 )
      icmp = jcmp
      return
      end

      real function kepl( r )
c the circular speed of an equivalent point mass
      implicit none
c
      real r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      kepl = 100
      if( r .gt. 0. )kepl = mnorm( 1 ) * cmpmas( 1 ) / sqrt( r )
      return
      end

      real*8 function vtot( r )
c the circular speed in a composite model
      implicit none
c
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 frtot, vdisc, vhalo
c
c local variables
      real*8 a
c
      vtot = 0
      if( r .gt. 0.d0 )then
        if( ncmp .gt. 1 )then
          a = 0
          do icmp = 1, ncmp
            if( disc( icmp ) )then
              a = a + vdisc( r )**2
            else
              a = a + vhalo( r )**2
            end if
          end do
          vtot = sqrt( a )
        else
          a = -r * frtot( r )
          vtot = sqrt( max( a, 0.d0 ) )
        end if
      end if
      return
      end

      real*8 function frtotl( lr )
      implicit none
c
c calling argument
      real*8 lr
c
c external
      real*8 frtot
c
c local variables
      real*8 r
c
      r = 10.**lr
      frtotl = frtot( r )
      return
      end

c      include 'inc/pgiden.f'
