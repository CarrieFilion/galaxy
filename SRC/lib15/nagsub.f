      subroutine nagsub
c utility routine to satisfy a linker by providing dummy externals
c   of the NAG library routines that could be used by the code
c The program will compile and execute, but if any external named in
c   this routine is called during execution, then the program will
c   abort with a detailed error message
      implicit none
c  Copyright (C) 2014, Jerry Sellwood
c
      character*6 routine
      integer n
      parameter ( n = 49 )
      integer ir, icall( n )
      save icall, ir
      data icall / n * 0 /, ir / 0 /
c
      entry c06fcf
c fft of complex data
      routine = 'C06FCF'
      ir = 1
      go to 1
c
      entry d01ahf
c 1D quadrature, adaptive, finite interval, strategy due to Patterson,
c   suitable for well-behaved integrands
      routine = 'D01AHF'
      ir = 2
      go to 1
c
      entry d01ajf
c 1D quadrature, adaptive, finite interval, strategy due to Piessens and
c   de Doncker, allowing for badly behaved integrands
      routine = 'D01AJF'
      ir = 3
      go to 1
c
      entry d01akf
c 1D quadrature, adaptive, finite interval, method suitable for oscillating
c   functions
      routine = 'D01AKF'
      ir = 4
      go to 1
c
      entry d01amf
c 1D quadrature, adaptive, infinite or semi-infinite interval
      routine = 'D01AMF'
      ir = 5
      go to 1
c
      entry d01apf
c 1D quadrature, adaptive, finite interval, weight function with end-point
c   singularities of algebraico-logarithmic type
      routine = 'D01APF'
      ir = 6
      go to 1
c
      entry d01aqf
c 1D quadrature, adaptive, finite interval, weight function $1/(x-c)$,
c   Cauchy principal value (Hilbert transform)
      routine = 'D01AQF'
      ir = 7
      go to 1
c
      entry d01baf
c 1D Gaussian quadrature
      routine = 'D01BAF'
      ir = 8
      go to 1
c
      entry d01baz
c 1D Gaussian quadrature utility
      routine = 'D01BAZ'
      ir = 9
      go to 1
c
      entry d01bbf
c Pre-computed weights and abscissae for Gaussian quadrature rules,
c   restricted choice of rule
      routine = 'D01BBF'
      ir = 10
      go to 1
c
      entry d01bcf
c Calculation of weights and abscissae for Gaussian quadrature rules,
c   general choice of rule 
      routine = 'D01BCF'
      ir = 11
      go to 1
c
      entry d01bdf
c 1D quadrature, non-adaptive, finite interval
      routine = 'D01BDF'
      ir = 12
      go to 1
c
      entry d01fbf
c Multi-D Gaussian quadrature over hyper-rectangle
      routine = 'D01FBF'
      ir = 13
      go to 1
c
      entry d01gaf
c 1D quadrature, integration of function defined by data values,
c    Gill–Miller method 
      routine = 'D01GAF'
      ir = 14
      go to 1
c
      entry d02baf
c solves the initial value problem for a first-order system of ordinary
c   differential equations using Runge–Kutta methods
      routine = 'D02BAF'
      ir = 15
      go to 1
c
      entry d02bbf
c solves the initial value problem for a first-order system of ordinary
c   differential equations using Runge–Kutta methods
      routine = 'D02BBF'
      ir = 16
      go to 1
c
      entry d02bgf
c ODEs, IVP, Runge–Kutta–Merson method, until a component attains given
c   value (simple driver) 
      routine = 'D02BGF'
      ir = 17
      go to 1
c
      entry d02bjf
c ODEs, IVP, Runge–Kutta method, until function of solution is zero,
c   integration over range with intermediate output (simple driver)
      routine = 'D02BJF'
      ir = 18
      go to 1
c
      entry d02bjw
c 1D Gaussian quadrature utility
      routine = 'D02BJW'
      ir = 19
      go to 1
c
      entry d02bjx
c 1D Gaussian quadrature utility
      routine = 'D02BJX'
      ir = 20
      go to 1
c
      entry d02cjf
c ODEs, IVP, Adams method, until function of solution is zero,
c   intermediate output (simple driver)
      routine = 'D02CJF'
      ir = 21
      go to 1
c
      entry d02pcf
c ODEs, IVP, Runge–Kutta method, integration over range with output
      routine = 'D02PCF'
      ir = 22
      go to 1
c
      entry d02pvf
c 1D Gaussian quadrature
      routine = 'D02PVF'
      ir = 23
      go to 1
c
      entry e01aaf
c Interpolated values, Aitken's technique, unequally spaced data, one variable 
      routine = 'E01AAF'
      ir = 24
      go to 1
c
      entry e01baf
c Interpolating functions, cubic spline interpolant, one variable
      routine = 'E01BAF'
      ir = 25
      go to 1
c
      entry e02bbf
c Evaluation of fitted cubic spline, function only
      routine = 'E02BBF'
      ir = 26
      go to 1
c
      entry e02bef
c Least-squares cubic spline curve fit, automatic knot placement
      routine = 'E02BEF'
      ir = 27
      go to 1
c
      entry e02caf
c Least-squares surface fit by polynomials, data on lines
      routine = 'E02CAF'
      ir = 28
      go to 1
c
      entry e02cbf
c Evaluation of fitted polynomial in two variables
      routine = 'E02CBF'
      ir = 29
      go to 1
c
      entry e02dcf
c Least-squares surface fit by bicubic splines with automatic knot placement,
c   data on rectangular grid
      routine = 'E02DCF'
      ir = 30
      go to 1
c
      entry e02def
c Evaluation of fitted bicubic spline at a vector of points
      routine = 'E02DEF'
      ir = 31
      go to 1
c
      entry e04abf
c Minimum, function of one variable using function values only
      routine = 'E04ABF'
      ir = 32
      go to 1
c
      entry e04jaf
c Minimum, function of several variables using function values only
      routine = 'E04JAF'
      ir = 49
      go to 1
c
      entry f02faf
c eigenvalues of a matrix
      routine = 'F02FAF'
      ir = 33
      go to 1
c
      entry f04aef
c Solution of real simultaneous linear equations with multiple right-hand
c   sides using iterative refinement (Black Box)
      routine = 'F04AEF'
      ir = 34
      go to 1
c
      entry f04jgf
c Least squares solution of real simultaneous linear equations with single
c   right-hand side
      routine = 'F04JGF'
      ir = 48
      go to 1
c
      entry f04atf
c Solution of real simultaneous linear equations, one right-hand side using
c   iterative refinement (Black Box) 
      routine = 'F04ATF'
      ir = 35
      go to 1
c
      entry g05caf
c Pseudo-random real numbers, uniform distribution over (0,1)
      routine = 'G05CAF'
      ir = 36
      go to 1
c
      entry g05cbf
c Initialize random number generating functions to give repeatable sequence
      routine = 'G05CBF'
      ir = 37
      go to 1
c
      entry g05cff
c Save state of random number generating routines
      routine = 'G05CFF'
      ir = 38
      go to 1
c
      entry g05cgf
c Restore state of random number generating routines
      routine = 'G05CGF'
      ir = 39
      go to 1
c
      entry m01daf
c Rank a vector, real numbers
      routine = 'M01DAF'
      ir = 40
      go to 1
c
      entry s14aaf
c $\Gamma$ function
      routine = 'S14AAF'
      ir = 41
      go to 1
c
      entry s14baf
c $\Gamma$ function
      routine = 'S14BAF'
      ir = 47
      go to 1
c
      entry s15aef
c Error function ${\rm erf}(x)$
      routine = 'S15AEF'
      ir = 42
      go to 1
c
      entry s15aff
c Dawson's integral
      routine = 'S15AFF'
      ir = 43
      go to 1
c
      entry s17aef
c Bessel function $J_0(x)$
      routine = 'S17AEF'
      ir = 44
      go to 1
c
      entry s17aff
c Bessel function $J_1(x)$
      routine = 'S17AFF'
      ir = 45
      go to 1
c
      entry x02amf
      routine = 'X02AMF'
      ir = 46
c
    1 if( icall( ir ) .eq. 0 )then
      print *, 'The program wants to use the external ' // routine
      print *, 'which is from the NAG library'
      print *
      print *, 'If the user has access to the NAG library, then the'
      print *, 'flag to link to the NAG library needs to be reset and'
      print *, 'the software reinstalled as described in the'
      print *, 'documentation'
      print *
      print *, 'If the user does not have NAG, then an equvalent'
      print *, 'routine will be needed.  The calling arguments and'
      print *, 'functionality of the routine to be substituted are'
      print *, 'described in the NAG on-line documentation:'
      print *,
     + 'http://www.nag.co.uk/numeric/fl/manual/html/FLlibrarymanual.asp'
      call
     + crash( 'NAGSUB', 'NAG library routine '// routine //' required' )
      end if
      icall( ir ) = icall( ir ) + 1
      return
      end
