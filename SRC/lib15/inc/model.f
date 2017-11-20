c
c information about all mass components
c
      integer ndfns, ndiscs, nhalos, ncmp
      logical cmprssd, fixed, ldfit, lgfld, numcal, singlr
      real*8 epsiln, rstar
      common / model / ncmp, ndiscs, nhalos, ndfns, fixed, singlr,
     +                 numcal, cmprssd, ldfit, lgfld
      common / modelr / epsiln, rstar
c external perturbers
      integer isphp
      logical extbar, exgnrc, exsphr, exspir, quapp
      real abar, alpha2, bbar, beta2, bphase, cbar, epspi, eptbr, mptbr
      real Omegab, Omegsp, sphase
      common / extcmp / exsphr, extbar, mptbr, eptbr, abar, bbar, cbar,
     +                  alpha2, beta2, quapp, Omegab, bphase, exspir,
     +                  Omegsp, sphase, epspi, exgnrc, isphp
c
c information about separate mass components
c
      character*4 cdft( mcmp ), ctype( mcmp )
      integer icmp, idftyp( mcmp ), impar( mcmp ), imtyp( mcmp )
      integer jdfcs( 4, mcmp ), iztyp( mcmp )
      logical disc( mcmp ), dist( mcmp ), Lztapr( mcmp )
      logical sphrod( mcmp ), trnctd( mcmp )
      real cmpmas( mcmp ), comi( 6, mcmp ), dfcns( 3, mcmp )
      real fmass( mcmp ), indexi( mcmp ), indexo( mcmp )
      real rscale( mcmp ), rtrunc( mcmp ), z0init( mcmp )
      common / cmps / icmp, imtyp, impar, cmpmas, rscale, rtrunc,
     +                fmass, disc, z0init, iztyp, idftyp, dist, dfcns,
     +                jdfcs, comi, sphrod, trnctd, Lztapr, indexi,
     +                indexo
      common / cmpsc /ctype, cdft
      real*8 DFnorm( mcmp )
      real*8 Lzcrt( mcmp ), Lzmno( mcmp ), Lzrmn( mcmp )
      real*8 Lztmx( mcmp ), Lztmn( mcmp ), mnorm( mcmp ), znorm( mcmp )
      common / cmps2 / Lztmn, Lztmx, Lzcrt, Lzrmn, Lzmno, znorm,
     +       DFnorm, mnorm
      real cmppar( 2, mcmp )
      common / cmps3 / cmppar
c radial and energy bounds for the entire system
      real*8 Emaxe, Emine, rhole, rmax, Phimax
      common / bounds / rmax, Phimax, Emaxe, rhole, Emine
c constants required for various DFs
      integer mmiy, nmiy, mmiy1, nmiy1
      logical kcold, retract
      real*8 cmiy, conc, evbeta, evsigma, omegak, omeg21, Phi0
      real*8 rtidal, rhofac, sigmai2, sigmak, amiy( 2, 20 )
      common / dfcons / Phi0, kcold, retract, omegak, omeg21,
     +                  sigmai2, rtidal, conc, sigmak, rhofac,
     +                  mmiy, nmiy, evbeta, evsigma, cmiy, mmiy1, nmiy1,
     +                  amiy
c tables for isotropic DFs
      integer iusedf, jusedf
      real*8 mrtab( nrdd ), rtab( nrdd ), psitab( nrdd )
      common / dftabs / iusedf, jusedf, mrtab, rtab, psitab
c adiabatic halo compression
      integer iccmp( 2 ), icomp, idfn0( 2 ), ihalo0( 2 ), ircmp( 2 )
      integer irigid, iterad, iusead, ncomp, nradad, nrigid
      logical r0core, r1cusp, r2cusp, snglph
      real*8 cc0( 2 ), dfb0( 2 ), dfm0( 2 ), Emaxo( 2 )
      common / adiabi / iusead, iterad, nradad, ncomp, nrigid, ihalo0,
     +                  idfn0, iccmp, ircmp, icomp, irigid
      common / adiabl / r0core, r1cusp, r2cusp, snglph
      common / adiabr / cc0, dfm0, dfb0, Emaxo
