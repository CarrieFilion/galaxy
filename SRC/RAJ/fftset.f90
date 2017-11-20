      subroutine fftset(n, m)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Routine identifier                                            tk

!     Purpose

!       To set up in /factor/ a table of multipliers for the fft package
!       and (optionally) to set up a re-ordering table for each distinct
!       length of transform required

!     Parameters

!     n   -  the number of lengths of transform required

!     m   -  an integer array of length  m,  holding the indices for the
!            transfer lengths on entry.  its values do not change.

!     Note:  only one multiplier table is set up, its length being that
!            for the longest transform.  shorter transforms use a subset
!            of this table, starting at its first entry.

      use constants
      
      use interfac, local => fftset

      use twiddle_factors
      
      integer              :: n, m(:), i, j, k, l, p, q, r, st, pt
      integer, save        :: number = 0
      double precision     :: x

!     Find maximum m and tag repeated values

      pt = 1
      if(n.gt.1) then
         do i = 2, n
            if(m(pt).lt.m(i)) pt = i
            do j = i, n
               if((m(j).eq.m(i - 1)).and.(m(j).gt.0)) m(j) = 1 - i
            end do
         end do
      end if
      st = n + 1

!     Cycle over required index arrays

      do i = 1, n
         if(m(i).lt.0) then

!           Plant link if already available

            j = - m(i)
            m(i) = m(j)

	 end if
      end do

!     Check space for multipliers

      l = m(pt)
      r = 2**(l - 1)
      if(allocated(fac)) deallocate(fac)
      allocate(fac(r))

!     Set up multiplier array

      p = 1
      q = 2
      st = 3
      fac(1) = real(2**(m(pt) - 2), kind=kin)
      fac(2) = fac(1)
      x = 2.0d0*fac(1)
      do j = 3, l
         do k = 1, p
            fac(st) = 0.5d0*fac(q)
            fac(st + 1) = x - fac(st)
            st = st + 2
            q = q + 1
         end do
         p = 2*p
      end do
      x = pi_over_2/x
      st = st - 1
      do j = 1, st
         fac(j) = cos(x*fac(j))
      end do

!     Finish.

      return
      end
