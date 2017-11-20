      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c    D1MACH ( 3) = B**(-T), the smallest relative spacing.
c    D1MACH ( 4) = B**(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer i

      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      else if ( i .eq. 1 ) then
        d1mach = 4.450147717014403D-308
      else if ( i .eq. 2 ) then
        d1mach = 8.988465674311579D+307
      else if ( i .eq. 3 ) then
        d1mach = 1.110223024625157D-016
      else if ( i .eq. 4 ) then
        d1mach = 2.220446049250313D-016
      else if ( i .eq. 5 ) then
        d1mach = 0.301029995663981D+000
      else if ( 5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      end if

      return
      end
