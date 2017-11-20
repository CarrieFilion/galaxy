      character*4 function htype( ihal )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the 4 character code of the halo type having the input number
c
c calling argument
      integer ihal
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local array
      integer mhalos
      parameter ( mhalos = 34 )
      character*4 halo( mhalos )
c
c halo types are:
c    1 - NONE   none
c    2 - ASDI   frozen disc
c    3 - UNIS   uniform sphere
c    4 - PLUM   Plummer model
c    5 - FE     Fall/Efstathiou model
c    6 - ISOT   cored isothermal
c    7 - BSS    Bahcall/Schmidt/Soniera model for the Galaxy
c    8 - 3198   NGC 3198 fit
c    9 - JAFF   Jaffe model
c   10 - KEPL   Kepler potential
c   11 - POWE   power law
c   12 - RQUA   r**1/4
c   13 - SCHM   modified Schmidt model
c   14 - DFIT   iterative DF model (properties not predictable analytically)
c   15 - TEDS   Thomasson etal
c   16 - KING   King model
c   17 - POLY   Spherical polytrope
c   18 - POWC   power law with a harmonic core
c   19 - KKSP   Kuz'min-Kutuzov oblate isochrone
c   20 - HERN   Hernquist model
c   21 - UNKN   Unknown halo type (properties not predictable analytically)
c   22 - AGME   Aguilar-Merritt 1/r cutoff model
c   23 - ADIA   adiabatically compressed halo
c   24 - ISOC   Henon's isochrone
c   25 - NFW    Navarro, Frenk & White
c   26 - SISP   singular isothermal sphere
c   27 - VKA1   fit to Valenzuela & Klypin (2003) model A_1
c   28 - MYNG   Miyamoto-Nagai flattnd Plummer sphere or thicknd Kuzmin disk
c   29 - MTAB   tabulated M(r)
c   30 - EXPH   2 comp spherical exp used by Donner & Thomasson A&A v290 p785
c   31 - ISOG   cored isoth with Gssn cutoff: Hernquist 1993, ApJS, v86, p389
c   32 - EINA   Einasto profile
c   33 - BURK   Burkert profile: ApJL, v447, L25
c   34 - CUBI   a finite radius bulge with a cubic density profile
c
      data halo / 'NONE', 'ASDI', 'UNIS', 'PLUM', 'FE  ', 'ISOT',
     +            'BSS ', '3198', 'JAFF', 'KEPL', 'POWE', 'RQUA',
     +            'SCHM', 'DFIT', 'TEDS', 'KING', 'POLY', 'POWC',
     +            'KKSP', 'HERN', 'UNKN', 'AGME', 'ADIA', 'ISOC',
     +            'NFW ', 'SISP', 'VKA1', 'MYNG', 'MTAB', 'EXPH',
     +            'ISOG', 'EINA', 'BURK', 'CUBI' /
c
      nhalos = mhalos
      if( ( ihal .gt. 0 ) .and. ( ihal .le. mhalos ) )then
        htype = halo( ihal )
      else
        print *, 'ihal =', ihal
       call crash( 'HTYPE', 'Impossible calling argument' )
      end if
      return
      end
