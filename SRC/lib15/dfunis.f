      real*8 function dfunis( E, L )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c DF for an uniform sphere (from Polyachenko & Shukhman - see BT08 prob. 4.11)
c
c calling arguments
      real*8 E, L
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real*8 arg, tiny
      parameter ( tiny = 1.d-10 )
c
      dfunis = 0
      arg = ( L / rscale( icmp ) )**2
      if( arg .lt. -Phimax )then
        arg = arg + 2. * ( E - Phi0 )
        if( arg .ge. tiny )then
          arg = sqrt( arg )
          dfunis = cmpmas( icmp )**2 * dfnorm( icmp ) / arg
        end if
      end if
      return
      end
