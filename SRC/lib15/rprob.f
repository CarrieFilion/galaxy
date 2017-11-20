      real*8 function rprob( rin )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns a radius from the inverse mass distribution of the desired model
c
c calling argument
      real*8 rin
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gmassd, gmassh, gmassi, gsigmt, ranuni
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir
      real*8 err, fmax, reqf, r1, r2, s, tol, x, y
c
      rprob = -1
      if( disc( icmp ) )then
c no extra tapers
        if( .not. Lztapr( icmp ) )then
c Kuz'min/Toomre disc
          if( ( ctype( icmp ) .eq. 'KT  ' ) .or.
     +        ( ctype( icmp ) .eq. 'SOFT' ) )then
            do while ( ( rprob .gt. rin ) .or. ( rprob .lt. rhole ) )
              rprob = 1. / ( 1. - ranuni( 0. ) )**2
              rprob = sqrt( rprob - 1. )
            end do
            rprob = rprob * rscale( icmp )
c Maclaurin/Freeman/Kalnajs disc
          else if( ctype( icmp ) .eq. 'MFK ' )then
            do while ( rprob .lt. rhole )
              rprob = ranuni( 0. )
              rprob = ( 1. - rprob )**.666666666667
              rprob = sqrt( 1. - rprob )
            end do
            rprob = rprob * rscale( icmp )
c Mestel/Toomre/Zang disc
          else if( ctype( icmp ) .eq. 'MTZ ' )then
            do while ( rprob .lt. rhole )
              rprob = rin * ranuni( 0. )
            end do
            rprob = rprob * rscale( icmp )
c Gaussian disc
          else if( ctype( icmp ) .eq. 'GAUS' )then
            x = gmassd( rin ) * ranuni( 0. )
            rprob = sqrt( -2. * log( 1. - x ) )
            rprob = rprob * rscale( icmp )
c Rybicki disc
          else if( ctype( icmp ) .eq. 'RYBI' )then
            do while ( rprob .lt. rhole )
              fmax = sqrt( 1. + rin**2 ) - 1.
              rprob = sqrt( ( fmax * ranuni( 0. ) + 1. )**2 - 1. )
            end do
            rprob = rprob * rscale( icmp )
c discs for which m(r) is known but cannot be inverted for r
          else if( ( ctype( icmp ) .eq. 'EXP ' ) .or.
     +             ( ctype( icmp ) .eq. 'SAND' ) .or.
     +             ( ctype( icmp ) .eq. 'POWE' ) .or.
     +             ( ctype( icmp ) .eq. 'COSI' ) .or.
     +             ( ctype( icmp ) .eq. 'DOTH' ) )then
            reqf = ranuni( 0. )
            r1 = rhole
            r2 = rin
            fmax = gmassi( r2 )
            tol = 1.e-6
            ind = 1
            ir = 0
            ifail = 0
            do while ( ind .ne. 0 )
              call fndzro( r1, r2, err, tol, ir, w, ind, ifail )
              err = gmassd( r1 ) / fmax - reqf
            end do
            rprob = r1
          end if
c discs for which only the surface density is known
        else
          fmax = gsigmt( rhole )
          s = -1
          do while ( s .lt. 0. )
            rprob = -1
            do while ( ( rprob .lt. rhole ) .or. ( rprob .gt. rin ) )
              x = rin * ranuni( 0. )
              y = rin * ranuni( 0. )
              rprob = sqrt( x * x + y * y )
            end do
            s = gsigmt( rprob ) / fmax - ranuni( 0. )
          end do
        end if
      else
c uniform sphere
        if( ctype( icmp ) .eq. 'UNIS' )then
          rprob = rin * ranuni( 0. )**.3333333333333333333
          rprob = rprob * rscale( icmp )
c Hernquist halo
        else if( ctype( icmp ) .eq. 'HERN' )then
          rprob = 2 * rin
          do while ( rprob .gt. rin )
            x = ranuni( 0. )
            rprob = ( sqrt( x ) + x ) / ( 1. - x )
          end do
          rprob = rprob * rscale( icmp )
c generalized Aguilar-Merritt sphere
        else if( ctype( icmp ) .eq. 'AGME' )then
          x = 1. / ( 3. + dfcns( 3, icmp ) )
          rprob = rin * ranuni( 0. )**x
          rprob = rprob * rscale( icmp )
c halos for which m(r) is known but cannot be inverted for r
        else if( ( ctype( icmp ) .eq. 'NFW ' ) .or.
     +           ( ctype( icmp ) .eq. 'ADIA' ) .or.
     +           ( ctype( icmp ) .eq. 'ISOG' ) )then
          reqf = ranuni( 0. )
          r1 = rhole
          r2 = rin
          fmax = gmassh( r2 )
          tol = 1.e-6
          ind = 1
          ir = 0
          ifail = 0
          do while ( ind .ne. 0 )
            call fndzro( r1, r2, err, tol, ir, w, ind, ifail )
            err = gmassh( r1 ) / fmax - reqf
          end do
          rprob = r1
        end if
      end if
      if( rprob .lt. 0.d0 )call crash( 'RPROB',
     +                                        'Unknown mass component' )
      return
      end
