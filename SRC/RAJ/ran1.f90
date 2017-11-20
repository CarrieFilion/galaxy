      double precision function ran1(idum)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Purpose:

!        Returns a uniform random deviate between 0.0 and 1.0.

!     Note:  This routine is adapted from the routine of the same name
!            given on pages 196-197 of Numerical Recipes (Fortran), Press, 
!            Flannery, Teukolsky and Vetterling, 1989, Cambridge University
!            Press.

!     Parameter:

!     iseed - seed value on entry - reset by routine.  set to any
!             negative value to initialise or re-initialise the
!             sequence.

      use interfac, local => ran1
      
!     Local variables:

      real               :: r(97)
      integer            :: idum, iff, ix1, ix2, ix3, j
      integer, save      :: number = 0
      integer, parameter :: m1 = 259200, ia1 = 7141, ic1 = 54773
      integer, parameter :: rm1 = 1.0/m1, m2 = 134456, ia2 = 8121
      integer, parameter :: ic2 = 28411, rm2 = 1.0d0/m2
      integer, parameter :: m3 = 243000, ia3 = 4561, ic3 = 51349

!     Initialise on first call even if idum is not negative.

      iff = 0

!     Test for initialisation.

      if((idum.lt.0).or.(iff.eq.0)) then
        iff = 1

!       Seed the first routine.

        ix1 = mod(ic1 - idum, m1)
        ix1 = mod(ia1*ix1 + ic1, m1)

!       Use it to seed the second and third.

        ix2 = mod(ix1, m2)
        ix1 = mod(ia1*ix1 + ic1, m1)
        ix3 = mod(ix1, m3)

!       Fill the table with uniform sequential deviates generated
!       by the first two routines.

        do j = 1, 97
           ix1 = mod(ia1*ix1 + ic1, m1)
           ix2 = mod(ia2*ix2 + ic2, m2)

!          Combine low and high order pieces.

           r(j) = (float(ix1) + float(ix2)*rm2)*rm1
        end do
        idum = 1
      end if

!     Except when initialising this is where we start.  Generate the
!     next number for each sequence.

      ix1 = mod(ia1*ix1 + ic1, m1)
      ix2 = mod(ia2*ix2 + ic2, m2)
      ix3 = mod(ia3*ix3 + ic3, m3)

!     Use the third sequence to get a number between 1 and 97.

      j = 1 + (ix3 + 97)/m3
      if((j.gt.97).or.(j.lt.1)) then

!       Error condition.

        write(*, '('' Error in random number generator  ran1'')')
        stop 'Error in random number generator found in ran1'
      end if

!     Return the table entry.

      ran1 = r(j)

!     Refill table.

      r(j) = (float(ix1) + float(ix2)*rm2)*rm1

!     Finish.

      return
      end



