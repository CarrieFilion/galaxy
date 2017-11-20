      subroutine jstrns( str, symb, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c look up table to "translate" TEX style Greek character strings to
c   the format used by PGPLOT
c
c calling arguments
      character*(*) str
      character*1 symb
      integer n
c
c local arrays
      character ch( 48 )*1, greek( 48 )*7
      integer l( 24 )
c
c local variables
      integer i, j, m
      logical ok
c
c translation tables - SGS
c      data ch / 'a', 'b', 'c', 'd', 'e', 'f', 'g', '{', 'i', 'j', 'k',
c     +          'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', '|', 'v',
c     +          'w', 'x',
c     +          'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
c     +          'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
c     +          'W', 'X' /
c translation tables - pgplot
      data ch / 'a', 'b', 'g', 'd', 'e', 'z', 'y', 'h', 'i', 'k', 'l',
     +          'm', 'n', 'c', 'o', 'p', 'r', 's', 't', 'u', 'f', 'x',
     +          'q', 'w',
     +          'A', 'B', 'G', 'D', 'E', 'Z', 'Y', 'H', 'I', 'K', 'L',
     +          'M', 'N', 'C', 'O', 'P', 'R', 'S', 'T', 'U', 'F', 'X',
     +          'Q', 'W' /
      data greek / 'alpha  ', 'beta   ', 'gamma  ', 'delta  ',
     +  'epsilon', 'zeta   ', 'eta    ', 'theta  ', 'iota   ',
     +  'kappa  ', 'lambda ', 'mu     ', 'nu     ', 'xi     ',
     +  'omicron', 'pi     ', 'rho    ', 'sigma  ', 'tau    ',
     +  'umicron', 'phi    ', 'chi    ', 'psi    ', 'omega  ',
     +  'Alpha  ', 'Beta   ', 'Gamma  ', 'Delta  ',
     +  'Epsilon', 'Zeta   ', 'Eta    ', 'Theta  ', 'Iota   ',
     +  'Kappa  ', 'Lambda ', 'Mu     ', 'Nu     ', 'Xi     ',
     +  'Omicron', 'Pi     ', 'Rho    ', 'Sigma  ', 'Tau    ',
     +  'Umicron', 'Phi    ', 'Chi    ', 'Psi    ', 'Omega  ' /
      data l / 5, 4, 5, 5, 7, 4, 3, 5, 4, 5, 6, 2, 2, 2, 7, 2, 3, 5, 3,
     +         7, 3, 3, 3, 5 /
c
      m = len( str )
      ok = .false.
c
      do i = 1, 24
        j = l( i )
        if( j .le. m )then
          if( str( 1:j ) .eq. greek( i )( 1:j ) )then
            write( symb, '( a1 )' )ch( i )
            n = j
            ok = .true.
          else if( str( 1:j ) .eq. greek( i + 24 )( 1:j ) )then
            write( symb, '( a1 )' )ch( i + 24 )
            n = j
            ok = .true.
          end if
        end if
      end do
c check that symbol was found
      if( .not. ok )then
        print *, 'string is ' // str( 1:m )
        call crash( 'JSTRNS', 'Symbol not found' )
      end if
      return
      end
