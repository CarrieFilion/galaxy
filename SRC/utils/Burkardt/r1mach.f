      function r1mach ( i )

c*********************************************************************72
c
cc R1MACH returns single precision real machine constants.
c
c  Discussion:
c
c    Assume that single precision real numbers are stored with a mantissa 
c    of T digits in base B, with an exponent whose value must lie 
c    between EMIN and EMAX.  Then for values of I between 1 and 5, 
c    R1MACH will return the following values:
c
c      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
c      R1MACH(3) = B**(-T), the smallest relative spacing.
c      R1MACH(4) = B**(1-T), the largest relative spacing.
c      R1MACH(5) = log10(B)
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
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 5.
c
c    Output, real R1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      real r1mach

      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      else if ( i .eq. 1 ) then
        r1mach = 1.1754944E-38
      else if ( i .eq. 2 ) then
        r1mach = 3.4028235E+38
      else if ( i .eq. 3 ) then
        r1mach = 5.9604645E-08
      else if ( i .eq. 4 ) then
        r1mach = 1.1920929E-07
      else if ( i .eq. 5 ) then
        r1mach = 0.3010300E+00
      else if ( 5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      end if

      return
      end
