      character*4 function dtype( idsc )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the 4 character code of the disk type having the input number
c
c calling argument
      integer idsc
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local array
      integer mdiscs
      parameter ( mdiscs = 22 )
      character*4 tdisc( mdiscs )
c
c disc types are:
c      1 - None
c      2 - KT          Kuz'min-Toomre
c      3 - SOFT        Softened Kuz'min-Toomre
c      4 - MFK         Maclaurin-Freeman-Kalnajs
c      5 - MTZ         Infinite Mestel-Toomre-Zang
c      6 - EXP         Exponential
c      7 - SC          Sellwood-Carlberg
c      8 - ISOC        Isochrone
c      9 - HUNT        Hunter's family
c     10 - SAND        Sanders
c     11 - POWE        Power law
c     12 - COSI        Cosine
c     13 - GAUS        Gaussian
c     14 - FMES        Finite Mestel
c     15 - ARBT        Arbitrary - given in tabular form
c     16 - RYBI        Rybicki disc - see Evans & Collett MNRAS v264 p353
c     17 - SPLY        Sawamura's polynomial discs - PASJ v40 p279
c     18 - PPLY        Polyachenko's polytropic discs - Astr Lett v20 p416
c     19 - UNKN        Unknown type - properties not predictable
c     20 - COMP        Composite exponential + Gaussian
c     21 - MTAB        M(r) table
c     22 - DOTH        Donner-Thomasson 2 differenced exps A&A v290 p785
c
      data tdisc / 'NONE', 'KT  ', 'SOFT', 'MFK ', 'MTZ ', 'EXP ',
     +             'SC  ', 'ISOC', 'HUNT', 'SAND', 'POWE', 'COSI',
     +             'GAUS', 'FMES', 'ARBT', 'RYBI', 'SPLY', 'PPLY',
     +             'UNKN', 'COMP', 'MTAB', 'DOTH' /
c
      ndiscs = mdiscs
      if( ( idsc .gt. 0 ) .and. ( idsc .le. mdiscs ) )then
        dtype = tdisc( idsc )
      else
        print *, idsc
        call crash( 'DTYPE', 'Impossible calling argument' )
      end if
      return
      end
