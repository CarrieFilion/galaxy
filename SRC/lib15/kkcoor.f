      subroutine kkcoor( r, z, lambda, nu )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes spheroidal coordinates from cylindrical polars for the
c   Kuz'min-Kutuzov oblate spheroid.  The value of a is hrad and the
c   axis ratio c/a is haloc; both are stored in / model / .  The
c   notation is that used by Dejonghe & de Zeeuw in Ap J v333 p90 (1988).
c
c calling arguments
      real*8 lambda, nu, r, z
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real*8 a2, b, c, c2, haloc, hrad
c
      hrad = rscale( icmp )
      haloc = dfcns( 3, icmp )
      a2 = hrad * hrad
      c2 = a2 * haloc * haloc
c do special cases first
      if( r * z .eq. 0.d0 )then
        lambda = c2 + z * z
        nu = a2 + r * r
      else
c general case
        b = a2 + c2 + r * r + z * z
        c = c2 * ( a2 + r * r ) + a2 * z * z
        if( 4. * c .gt. b * b )then
          print *, 'Error in KKCOOR'
          print *, 'Calling arguments are:', r, z
          print *, 'Derived values are:', b * b, 4.* c
          stop
        end if
        c = sqrt( b * b - 4. * c )
        lambda = .5 * ( b + c )
        nu = .5 * ( b - c )
      end if
      return
      end
