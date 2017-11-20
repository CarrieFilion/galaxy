      subroutine lsmon2( nm2, xc, sumsq, gc, istate, gpjrnm, cond,
     +                   posdef, niter, nf, iw, liw, w, lw )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c part of mode fitting software
c
c Monitoring routine called by E04KDF
c It reports the current position and value of the residuals and inquires
c    at infrequent intervals whether the user wishes to abort the search
c
c calling arguments
      integer liw, lw, nf, niter, nm2, istate( nm2 ), iw( liw )
      logical posdef
      real*8 cond, gpjrnm, sumsq, gc( nm2 ), w( lw ), xc( nm2 )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      include 'inc/jscmmn.f'
c
c external
      logical gtlogl
c
      print 200, niter, nf, sumsq, xc
 200  format( i4, 'th iteration', i6, ' function calls, current sumsq',
     +        e12.4 / ' freqs', 40f8.4 )
      if( null )then
c this many iterations implies fit is not converging
        if( niter .ge. 50 * nm2 )iabort = -1
      else
c check whether to continue
        if( ( mod( niter, 5 * nm2 ) .eq. 0 ) .and.
     +      ( niter .ne. 0 ) )then
          if( .not. gtlogl( 'Do you want to continue?' ) )iabort = -1
        end if
      end if
      return
      end
