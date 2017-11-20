c
c variables for parallel processing
c
      integer numprocs, myid
      logical parallel, master
      common / archit / parallel, numprocs, myid, master
c
c sizes of some allocatable arrays in module aarrays
c
      integer bhmt, bhnlev, lptcls, lsfpl, lsfpr, ndrct
      common / bsizes / lptcls, lsfpl, lsfpr, ndrct, bhmt, bhnlev
c
c method for determining accelerations
c
      integer igrid( mgrids ), ncode, ncoor, ndimen, ngrid
      logical bht, c2d, c3d, gtype( 0:ncodes ), hybrid, dr3d, ldrct
      logical lgrd, lhyb2, lnag, lsfp, noslfg, p2d, p3a, p3d
      logical s3d, scf, sf2d, sf3d, threed, twod, twogrd
      common / code / twod, threed, ncoor, lgrd, hybrid, ngrid, igrid,
     +                ncode, noslfg, p2d, c3d, sf2d, p3a, p3d, s3d,
     +                sf3d, c2d, scf, dr3d, bht, twogrd, ndimen, lsfp,
     +                ldrct, lnag, lhyb2
      equivalence ( gtype( 0 ), noslfg )
c
c run number, step count analysis intervals, etc
c
      integer iphys, ipict, irun, isfld, istep, isubst !, jactiv( ncmp )
      logical lprint, phys, potl
      real cvers
      common / print / cvers, irun, istep, lprint, phys, iphys, ipict,
     +                 potl, isubst, isfld !, jactiv
c
c mesh pointers
c
c      integer lpt( mzones, mgrids )
      integer lstep( 2, 2:mzones ), mstep( mzones )
      common / point / lstep, mstep !, lpt
c
c particles and linked lists
c
      integer ilist, inext, iplane, izone, jgrid, jlist, jplane, lpf
      integer mzone, nin, nlists, noff, nout, nshort, nwpp
      integer islist( 2, lists, mnodes ), ncen( mcmp )
      logical lzones, tcent, uqmass
      common / pstore / lpf, nwpp, uqmass, noff, ncen, nshort, tcent,
     +                  lzones, nin, nout, nlists, izone, inext,
     +                  mzone, jgrid, ilist, jlist, iplane, jplane,
     +                  islist
      integer ncenm( mcmp ), noffm, nshortm
      logical lurand
      common / pstore2 / noffm, nshortm, ncenm, lurand
c
c locations and nature of particle populations
c
      integer igrd( mcmp ), npr( mcmp ), nsp( mcmp )
      logical quiet( mcmp ), rigidp( mcmp ), smr( mcmp ), testp( mcmp )
      logical netcom, netlmo
      common / grdpop / nsp, igrd, npr, quiet, smr, rigidp, testp,
     +                  netcom, netlmo
c
c time step scheme and zoning
c
      logical lfrv
      common / intgtr / lfrv
c
c zones and guard radii
c
      integer jrad( mzones ), krad( mzones ), nzones, nstep( mzones )
      logical sphbnd
      real tzone( mzones ), zbr( mzones )
      common / zones / nzones, nstep, zbr, tzone, jrad, krad, sphbnd
c
      integer nguard, nsub( mguard )
      real gbr( mguard )
      common / gbound / nguard, nsub, gbr
c
c scale factors between internal and physical units and other constants
c
      integer nbod, tsoft
      logical fixaxi, heavies, pertbn, rigidh, suppl
      real gvfac, lscale, pmass, rguard, softl, ts
      common / scales / nbod, fixaxi, rigidh, suppl, pertbn, heavies,
     +                  tsoft
      common / scalesr / lscale, ts, pmass, softl, gvfac, rguard
c
c external perturber
c
      integer nptrb, npred
      logical gshift, sumptb
      real accptb( 3 ), peptrb
      real*8 accpz( 3, 2, mzones ), pzacc( 3, mzones )
c      real epsptb, mptrb
      real vgal( 3 ), vptrb( 3 ), xptrb( 3 ), xshift( 3 ), pertbr( 16 )
      common / prtbn / sumptb, nptrb, xptrb, vptrb, accptb, xshift,
     +                 vgal, peptrb, npred, gshift, accpz, pzacc
c, mptrb, epsptb
      equivalence ( pertbr( 1 ), xptrb( 1 ) )
c
c recentering grid options
c
      integer icstp, jebind( mgrids ), kci, loch( mgrids ), nebind
      logical centrd, lbind, ltroid, lheavy
      real xcen( 3, mgrids )
      common / centroid / centrd, icstp, kci, xcen, lbind, ltroid,
     +                    nebind, jebind, lheavy, loch
c
      integer ipast( npast )
      real cenfit( 3, 3, mgrids ), pastp( 3, npast, mgrids )
      real xcpred( 3, mzones, mgrids )
      common / cnpath / pastp, cenfit, xcpred, ipast
c
c analysis
c
      logical angm, danl, dnst, frqs, grey, intg, lgsp, lsave, lval
      logical moit, moni, orbs, plot, pntl, rhor, rngs, s3dc, sate
      logical scfc, sfpc, sigr, sphb, vfld, vflh, vlgs, zanl, zprf
c, mode
      logical dattyp( ndattyp )
      common / reslts / lsave, plot, dnst, pntl, frqs, lgsp, angm, vfld,
     +                  intg, orbs, grey, danl, moni, sfpc, zanl, vflh,
     +                  sphb, zprf, rngs, scfc, rhor, sate, vlgs, moit,
     +                  sigr, s3dc, lval
c, mode
      equivalence ( plot, dattyp( 1 ) )
c
      integer nprop, mprho
      real drfac, drfrqs, hzfac, Lzfac, rbess, thfac
      common / anlcns / drfac, thfac, hzfac, Lzfac, rbess, nprop, mprho,
     +                  drfrqs
c
      real unit_L, unit_M, unit_T, unit_V
      logical scale_set
      common / physcl / scale_set, unit_L, unit_M, unit_T, unit_V
