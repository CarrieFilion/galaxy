      real*8 function qtoom( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns Toomre's Q parameter at the requested radius in the currently
c   defined model
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 akappa, gsigma, sigmau
c
c local variables
      real*8 a, s
c
      if( ( cdft( icmp ) .eq. 'SHUE' ) )then
        qtoom = dfcns( 1, icmp )
      else
        s = gsigma( r )
        if( s .eq. 0. )then
          qtoom = 0
c Newton's Bayesian prescription
        else
          if( cdft( icmp ) .eq. 'LUCY' )then
            a = dfcns( 2, icmp ) * ( dfcns( 1, icmp ) - r )
            a = exp( a )
            qtoom = 1.5 + .5 * ( ( a - 1. / a ) / ( a + 1. / a ) )
c the defining expression
          else
            qtoom = sigmau( r ) * akappa( r ) / ( 3.36 * s )
          end if
        end if
      end if
      return
      end
