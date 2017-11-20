      subroutine resrad( omega, m, res )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Routine to find co-rotation and 3 Lindblad resonance radii in the
c   theoretical mass model for a given Omega_p  = omega / m
c The resonances are:
c    (1) co-rotation
c    (2) outer Lindblad
c    (3) OILR
c    (4) IILR
c The routine flags the non-existence of a resonace by setting the
c    relevant radius to large or negative values
c    uses Burkardt's routine local_min
c
c calling arguments
      integer m
      real*8 omega, res( 4 )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      common / rdfres / am, rmn, ltheory
      logical ltheory
      real rmn
      real*8 am
c
c externals
      real*8 akappa, local_min, omegac, resrdf
      external resrdf
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir, ires, k
      logical lr
      real dr, r
      real*8 epsa, epsr, err, fmin, om, omegap, rm, rmx, sgn, x1, x2
c
      parameter ( rmx = 100. )
c
      ltheory = .false.
      if( .not. ltheory )then
        rmn = rinh( jgrid )
        call rewnd
        time = -1
        do while ( time .lt. tme( kt ) )
          call nextrec( 'FRQS', ifail )
          if( ifail .ne. 0 )then
            print *, 'Using theoretical frequencies in RESRAD'
            time = tme( kt )
            ltheory = .true.
          end if
        end do
      end if
c set default values
      res( 1 ) = rmx
      res( 2 ) = rmx
      res( 3 ) = -1.
      res( 4 ) = -1.
c skip if retrograde
      if( omega .ge. 0. )then
        am = m
        omegap = omega / am
        lr = .true.
c attempt to find resonances
        do ires = 1, 4
          x1 = .1 * rhole
          x1 = max( rm, 1.d-5 )
          x2 = max( 10. * rmax, rmx )
          if( ires .eq. 3 )then
c no ILRs possible within a uniformly rotating disc
c   and problematic for polynomial disc
            lr = ( ncmp .gt. 1 ) .or. ( ctype( 1 ) .ne. 'MFK ' ) .or.
     +                                ( ctype( 1 ) .ne. 'SPLY' )
            if( lr )then
c find maximum of omega - kappa / m
              if( ltheory )then
                epsr = 1.e-6
                epsa = epsr
                fmin = local_min( x1, x2, epsr, epsa, resrdf, rm )
              else
c just inspect tabulated values
                fmin = 0
                do while ( x1 .lt. x2 )
                  epsa = resrdf( x1 )
                  fmin = min( fmin, epsa )
                  x1 = x1 + drfrqs / lscale
                end do
              end if
c determine whether LRs are present - n.b. fmin is -ve!
              lr = omegap + fmin .lt. 0
            end if
          end if
c find radius of resonance
          if( lr )then
c revise range for ILRs only
            if( ires .gt. 2 )then
              x1 = rm
              x2 = rmx
              if( ires. eq. 4 )x2 = 0.
            end if
            epsa = 1.e-8
            sgn = 1.
            if( ires .gt. 2 )sgn = -1.
            ifail = 1
            ind = 1
            ir = 0
            do while ( ind .ne. 0 )
              call fndzro( x1, x2, err, epsa, ir, w, ind, ifail )
              if( ltheory )then
                om = omegac( x1 )
              else
c get Omega from data array
                r = x1 * lscale
                k = nint( ( r - rmn ) / drfrqs )
                if( k .lt. mr - 1 )then
                  dr = r - k * drfrqs
                  om = ( 1. - dr ) * wres( k + 1 ) + dr * wres( k + 2 )
                else
                  dr = real( mr ) * drfrqs / lscale
                  om = wres( mr ) * dr / x1
                end if
              end if
              err = omegap - om
              if( ires. gt. 1 )then
                if( ltheory )then
                  err = err - sgn * akappa( x1 ) / am
                else
c get kappa from data array
                  if( k .lt. mr - 2 )then
                    k = k + mr
                    om = ( 1. - dr ) * wres( k + 1 ) +
     +                          dr   * wres( k + 2 )
                  else
                    dr = real( mr - 1 ) * drfrqs / lscale
                    om = wres( 2 * mr - 1 ) * dr / x1
                  end if
                  err = err - sgn * om / am
                end if
              end if
            end do
c return only successful results
            if( ( ifail .eq. 0 ) .or. ( ifail .gt. 3 ) )res( ires ) = x1
            if( ifail .ne. 0 )print *,
     +  'ifail = ', ifail, ' from FNDZRO for ires =', ires, ' in RESRAD'
          end if
        end do
      end if
      return
      end
