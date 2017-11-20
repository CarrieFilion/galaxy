      subroutine new_unit( n )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to supply an unused logical unit number that can be associated
c   with a new file
c
c calling argument
      integer n
c
c local variable
      integer k
      save k
c
      data k / 6 /
c
      k = k + 1
      n = k
      return
      end
