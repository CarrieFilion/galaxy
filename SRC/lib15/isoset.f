      subroutine isoset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c sets up parameters for spherical isotropic models described in ch 4 of BT
c   and then converts to my units which set G = M_tot = scale radius = 1
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
      external isopsf
      real*8 gmisph, isrhop, rhoisp, Gammaf
c
c local array
       real*8 y( 2 )
c
c local variables
      integer i, ifail
      real*8 cn, haloc, mc, mtot, psifac, rs, s, s0, tol, zero
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'ISOSET',
     +                                        'Nonsense value of icmp' )
      haloc = dfcns( 3, icmp )
c integral for radial density profile is in dimensionless units
      mc = cmpmas( icmp )
      rs = rscale( icmp )
      cmpmas( icmp ) = 1
      rscale( icmp ) = 1
      rhofac = 1
c
      iusedf = 0
      jusedf = 0
c
c King model
      if( ctype( icmp ) .eq. 'KING' )then
        if( haloc .le. 0.d0 )call crash( 'ISOSET',
     +                                     'Central potential not set' )
c set initial constants
        sigmak = 1
        y( 1 ) = haloc
c polytrope by integration of the Lane-Emden equation (4-108c) of BT
      else if( ctype( icmp ) .eq. 'POLY' )then
        dfcns( 1, icmp ) = haloc - 1.5
        y( 1 ) = 1
      else
        call crash( 'ISOSET', 'Unrecognized halo type' )
      end if
c find tidal radius
      rtidal = 0
      s = 100
      y( 2 ) = 0
      tol = 1.d-8
      zero = 0
      ifail = 0
      call ode_toy( rtidal, s, 2, y, tol, 0.d0, 1, zero, isopsf, ifail )
c find scaled total mass - initialize table of psi(r) first
      s0 = rhoisp( 0.d0 )
      mtot = gmisph( rtidal )
      if( ctype( icmp ) .eq. 'KING' )then
c set s0 to King radius in scaled units (see notes)
        s0 = 3. / sqrt( isrhop( haloc ) )
        conc = log10( rtidal / s0 )
        if( master )write( no,
     +       '( '' Tidal radius and concentration index are'', 2f10.4 )'
     +                                                )rtidal / s0, conc
c rescale to units in which G = 1, total mass = mc & King radius = rs
        psifac = 4.d0 * pi * s0 * mc / ( mtot * rs ) ! = sigmak**2
        sigmak = sqrt( psifac )
        rtidal = rs * rtidal / s0
        rhofac = ( mc / mtot ) * ( s0 / rs )**3
        dfnorm( icmp ) = rhofac / ( 2.d0 * pi * psifac )**1.5
      else if( ctype( icmp ) .eq. 'POLY' )then
c see notes (herschel:/home/sellwood/models/models.tex)
        if( haloc .gt. 0.5d0 )then
c note: \Gamma( x + 1 ) = x!
          ifail = 0
          cn = ( 2. * pi )**1.5 * Gammaf( haloc - 0.5d0, ifail )
     +                                   / Gammaf( haloc + 1.d0, ifail )
        else
          cn = 1
        end if
        s0 = rtidal
        psifac = 4 * pi * s0 * mc / ( rs * mtot ) ! = Psi(0)
        rhofac = mc / ( mtot * psifac**haloc ) * ( s0 / rs )**3
        dfnorm( icmp ) = rhofac / cn
        rtidal = rs
        if( master )write( no,
     +   '( ''ISOSET completed for polytrope of index'', f10.4 )' )haloc
      else
        call crash( 'ISOSET', 'Unrecognized halo type 2' )
      end if
      rtrunc( icmp ) = rtidal
c restore mass and scale length and rescale table
      cmpmas( icmp ) = mc
      rscale( icmp ) = rs
      do i = 1, nrdd
        rtab( i ) = rtab( i ) * rs / s0
        psitab( i ) = psitab( i ) * psifac
        mrtab( i ) = mrtab( i ) * mc / mtot
      end do
      return
      end
