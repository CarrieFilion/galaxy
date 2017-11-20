      real*8 function dfwght( E, Lz )
c function to weight DF when unequal particle masses are requested
c  value of this function is the inverse of the desired particle weights
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      real*8 E, Lz
c
c dummy version for the linker
      call crash( 'DFWGHT', 'Need to define desired weight fn' )
c weight as inverse of angular momentum with a small offset
      dfwght = 1. / ( .01 + Lz )
      return
      end
