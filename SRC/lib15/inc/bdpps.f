c Global constants and variables for the Hernquist SCF code
      integer lbasmx, nbasmx
      parameter ( nbasmx = 100, lbasmx = 20 )
      real*8 anltilde( 0:nbasmx, 0:lbasmx ), c1( 1:nbasmx, 0:lbasmx )
      real*8 c2( 1:nbasmx, 0:lbasmx ), c3( 1:nbasmx )
      real*8 coeflm( 0:lbasmx, 0:lbasmx ), twoalpha( 0:lbasmx )
      real*8 ultrasp( 0:nbasmx, 0:lbasmx )
      real*8 ultraspt( 0:nbasmx, 0:lbasmx )
      real*8 ultrasp1( 0:nbasmx, 0:lbasmx )
      common / scfarrays / anltilde, c1, c2, c3, coeflm,
     +                     twoalpha, ultrasp, ultraspt, ultrasp1
      integer nbas, lmin, lskip
      common / scfint / nbas, lmin, lskip
c
c Global constants and variables for the Smooth-Field-Particle (SFP) code.
c NB `bdpps' stands for `biorthonormal density-potential pair set'.
c Originally created by David Earn
c
c parameters
c
c maximum allowed azimuthal- (m) and radial-order (n) in the basis expansion
      integer mostm, mostn
      parameter ( mostm = 4, mostn = 100 )
c Abel-Jacobi functions described by Kalnajs (1976)
c maximum permitted value of k
      integer mostkaj
      parameter ( mostkaj = 30 )
c maximum permitted number of planes for Hankel function solver
      integer mplanes
      parameter ( mplanes = 21 )
c maximum number of basis functions
      integer mostf
      parameter ( mostf = ( mostn + 1 ) * ( mostm + 1 ) )
c
c Abel-Jacobi functions described by Kalnajs (1976)
c The greatest power of r that may be required is therefore 
c   mostpower = 2*mostn + mostm + 2*mostkaj  (cf Kalnajs 1976, equation 22)
c
c
c common variables
c
      integer basis, lastf, maxm, maxn, minm, minn, ncontrib, nrtab
      integer msel( mostf ), nsel( mostf )
      logical linear
      real deltak, maxmaxr, maxr, minmaxr, newmaxr, threshold
      real sfpwt( 0:mostn )
      common / bdpps / basis, minm, maxm, minn, maxn, lastf, minmaxr,
     +                 maxmaxr, linear, threshold, ncontrib, maxr,
     +                 newmaxr, nsel, msel, sfpwt, nrtab
      equivalence ( deltak, basis )
c
      character*4 basset
      common / btype / basset
c
      integer maj, kaj, naj
      common / abljac / maj, kaj, naj
c
c GLOBAL VARIABLES IN COMMON BLOCK / bdpps /

c basset indicates the name of the basis set to be used

c Maxr is the maximum radius (in grid units) out to which the bdpps
c functions are defined (the `function edge').  Minmaxr is the minimum
c allowed value of maxr, and newmaxr is the grid radius of the outermost
c particle at the current step.  Somewhat confusingly, there is a
c different variable called rmax in common block /elims/.  Maxmaxr is
c the maximum allowed value of maxr: beyond maxmaxr particles are always
c dropped off the edge of the world.

c The number that uniquely specifies which basis set is being used.
c 0 <= basis <= mostkaj yields Abel-Jacobi functions with k=basis.
c Other sets will have higher values of basis because this number
c is used as a quick branch in basfun.f and basder.f .

c Max and min n and m selected for the current run.  Do not confuse
c maxm with mmax in common block /grid/.

c List of m>=0 coefficients being used.  (m<0 coefficients are always
c used as well.)  m>0 refers to real parts (cos terms) and m<0 to 
c imaginary parts (sin terms).

c Coefficients of bdpps basis functions.

c For specifying basis functions in a list.  The coefficient of
c the current function is coeff( nsel(currentf), msel(currentf) ).

c A flag indicating if no axisymmetric (m=0) basis functions have been
c selected.  If  fixed_axi=.true.  then the axisymmetric component of the
c force is calculated from a fixed potential using halfrc.

c A flag that is .true. during the linear epoch, as defined by
c subroutine cmonitor.

c The threshold value of the diagnostic calculated in cmonitor.
c When the diagnostic first exceeds this value the linear epoch
c is deemed over (value is set in sfpinput).

c The number of particles that contribute to the coefficients computed
c by findc at the current step.  This excludes particles in the
c central hole and particles outside maxmaxr.
c
c roots and coefficients for accurate Abel-Jacobi functions
c   the coefficients are of the polynomial in  x=(1-r^2)
c NB these arrays are really wasteful of space since only values below the
c   body diagonal are non-zero
      real*8 root( 0:mostm, 0:mostn, 1:mostn )
      real*8 xcoeff( 0:mostm, 0:mostn, 0:mostkaj )
      common / ajnumbers / root, xcoeff
c
c basis function values and derivatives at some radius
      real*8 dervals( mostf ), funvals( mostf )
      common / valfun / funvals, dervals
