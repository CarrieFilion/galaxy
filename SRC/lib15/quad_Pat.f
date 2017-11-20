      real*8 function quad_Pat( fun, a, b, epsr, ifail )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c adaptive quadrature of a given function - should work for most
c   well-behaved functions
c
c the result is returned through the function value
c
c arguments: fun   the function to be integrated, which must be a
c                  real*8 function of a single variable
c            a, b  the range of integration
c            epsr  the relative precision required for the result
c            ifail an error flag that is returned, a value of 1 indicates
c                  the requested accuracy was not achieved
c                      see Burkardt's source code
c
c   uses Burkardt's routine QUAD that implements Patterson's method
c
c calling arguments
      integer ifail
      external fun
      real*8 a, b, epsr
c
c local variables
      integer k, npts
      real*8 res( 8 )
c
      call quad( a, b, res, k, epsr, npts, ifail, fun )
c      print *, 'additional output from quad'
c      print *, 'order and number of function evaluations', k, npts
      quad_Pat = res( k )
      return
      end
