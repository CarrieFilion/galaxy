      character*4 function dftype( idf )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the 4 character code of the DF having the input number
c
c calling argument
      integer idf
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local array
      integer mdfn
      parameter ( mdfn = 28 )
      character*4 distf( mdfn )
c
c distribution function types are:
c     1 - None
c     2 - Kalnajs
c     3 - Lia (generalised Kalnajs for KT only)
c     4 - Miyamoto (for KT only)
c     5 - Zang (for MTZ only)
c     6 - Omega - Kalnajs's bizarre function for Omega models
c     7 - Merritt - for spherical Jaffe model
c     8 - Dejonghe - for spherical Plummer model
c     9 - Lucy (disc) - generalised Newton-Binney method
c    10 - Neil - any other abitrary file of particles!
c    11 - Kalnajs's Composite B_O
c    12 - King function for a King model
c    13 - Isoptropic spherical polytrope
c    14 - Wyn Evans's functions for power law discs
c    15 - Dejonghe-deZeeuw functions for Kuz'min-Kutuzov spheroids
c    16 - Functions for Toomre model 0 given by Evans & Collett (MNRAS 1993)
c    17 - Sawamura's DFs for his polynomial discs (PASJ 1988 v40 p279)
c    18 - Halo DF from _DFITER_
c    19 - Polyachenko's lowered polytropic DFs (Astr Lett 1994 v20 p416)
c    20 - Hernquist's isotropic DF (ApJ 1990 v356 p359)
c    21 - adiabatically compressed halo - defined iteratively
c    22 - Henon's spherical isochrone (eq 4-144 of BT)
c    23 - Eddington's isotropic form for a sphere
c    24 - singular isothermal sphere
c    25 - adiabatically compressed halo
c    26 - Shu DF - epicyclic approx for a generic disk model
c    27 - Uniform sphere from Polyachenko & Shukhman - see BT08 problem 4.11
c    28 - Jeans equation - as described by Hernquist, ApJS, v86, p389
      data distf / 'NONE', 'KALN', 'LIA ', 'MIYA', 'ZANG', 'OMEG',
     +             'MERR', 'DEJO', 'LUCY', 'NEIL', 'CPB0', 'KING',
     +             'POLY', 'EVAN', 'DJDZ', 'EVCO', 'SAWA', 'DFIT',
     +             'PPLY', 'HERN', 'ADIA', 'HNON', 'EDDI', 'SISP',
     +             'COMP', 'SHUE', 'USPS', 'JEAN' /
c
      ndfns = mdfn
      if( ( idf .ge. 1 ) .and. ( idf .le. mdfn ) )then
        dftype = distf( idf )
      else
        print *, idf
        call crash( 'DFTYPE', 'Impossible calling argument' )
      end if
      return
      end
