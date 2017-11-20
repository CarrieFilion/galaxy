      subroutine fftcsa(n, no, group, w, unity, halve)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Routine identifier                                            tl

!     Purpose

!        To perform a fourier cosine analysis along one axis of a
!        2 or 3 dimensional mesh.

!     Parameters:

!     n      -  (type integer) gives the transform length.

!     no     -  (type integer) gives the total number of transforms required.

!     group  -  (type integer) gives the number of transforms to be
!               evaluated in parallel at each stage.

!     w      -  a double precision array containing the data.

!     unity  -  (type logical) specifies if the data for each individual
!               transform is laid out sequentially in memory.

!     halve  -  (type logical) specifies if the end values should be halved
!               after transforming.

      use constants
      
      use interfac, local => fftcsa

      use twiddle_factors
      
!     Parameter definitions:
      
      integer, intent(in)             :: n, no, group
      integer, save                   :: number = 0
      logical, intent(in)             :: unity, halve
      double precision, intent(inout) :: w(n*no)
      
      integer st, st1, step, iw, i, j, p, q, pii, qi
      integer incr, ii, it, jump, lim, na, nb
      double precision c, s, u, v

!     Local variables:
      
      logical rev

!     Set initial values.

16    if(unity) then
        incr = n
        step = 1
      else
        incr = 1
        step = group
      end if
      st = 1

!     Cycle over transforms.

      do 15 ii = 1, no, group
      it = min0(group, no - ii + 1)
      do 1 i = 1, it
      na = n - 1
      lim = st + na*step
      if(halve) go to 17
      w(st) = 2.0d0*w(st)
      w(lim) = 2.0d0*w(lim)
17    lim = step + lim

!     Fold real section.

 2    p = st
      jump = na*step
      q = p + jump
      nb = na/2
      do 3 j = 1, nb
      u = w(p)
      w(p) = u + w(q)
      w(q) = u - w(q)
      p = p + step
3     q = q - step
      w(p) = 2.0d0*w(p)

!     Test scan complete.

      st1 = st + step + jump
      if(st1.eq.lim) go to 4

!     Re-order real group to complex form.

      p = st1
      q = p + jump - step
      do 5 j = 1, nb
      u = w(p)
      w(p) = w(q)
      w(q) = - u
      p = p + step
5     q = q - step

!     Adjust last real part.

      w(p) = - w(p)

!     Test scan complete.

      st1 = st1 + jump
      if(st1.eq.lim) go to 4

!     Fold first complex group.

!     Initialise and set first terms.

      p = st1
      q = p + jump
      qi = q
      pii = q + jump
      u = root2*w(q)
      w(q) = w(p) - u
      w(p) = w(p) + u
      if(nb.eq.1) go to 10

!     Fold intermediate terms.

      do 7 j = 2, nb
      p = p + step
      q = q - step
      pii = pii - step
      qi = qi + step
      u = recip_root2*(w(q) - w(qi))
      v = recip_root2*(w(q) + w(qi))
      w(q) = w(pii) - v
      w(qi) = w(p) - u
      w(p) = w(p) + u
 7    w(pii) = - w(pii) - v

!     Calculate end terms.

10    p = p + step
      pii = pii - step
      u = recip_root2*(w(p) - w(pii))
      w(pii) = w(p) - u
      w(p) = w(p) + u

!     Test scan complete.

      st1 = st1 + 2*jump
      if(st1.eq.lim) go to 4

!     Fold remaining groups.

      iw = 1
      rev = .true.
8     rev = .not.rev
      if(rev) go to 13
      iw = iw + 2
      c = fac(iw)
      s = fac(iw + 1)
      go to 14
13    u = c
      c = s
      s = u

!     Pointers and first term.

14    p = st1
      q = p + jump
      pii = q + jump
      qi = q
      u = w(q)/c
      w(q) = w(p) - u
      w(p) = w(p) + u
      if(nb.eq.1) go to 11

!     Fold intermediate terms.

      do 9 j = 2, nb
      p = p + step
      q = q - step
      pii = pii - step
      qi = qi + step
      u = c*w(q) - s*w(qi)
      v = s*w(q) + c*w(qi)
      w(q) = w(pii) - v
      w(pii) = - w(pii) - v
      w(qi) = w(p) - u
 9    w(p) = w(p) + u

!     Calculate end terms.

11    p = p + step
      pii = pii - step
      u = c*w(p) - s*w(pii)
      w(pii) = w(p) - u
      w(p) = w(p) + u

!     Next group if not end of scan.

      st1 = st1 + 2*jump
      if(st1.lt.lim) go to 8

!     Jump back for next fold.

 4    na = nb
      if(na.ne.1) go to 2

!     Final reduction.

6     st1 = st + 3*step
      iw = 1
      p = st + step
      u = w(p)
      w(p) = w(st) - u
      w(st) = w(st) + u
      na = (n - 1)/2
      nb = 2*step
      do 12 j = 2, na
      iw = iw + 1
      p = st1 + step
      u = w(p)/fac(iw)
      w(p) = w(st1) - u
      w(st1) = w(st1) + u
12    st1 = st1 + nb
      if(.not.halve) go to 1
      w(st) = 0.5d0*w(st)
      st1 = st + step
      w(st1) = 0.5d0*w(st1)
1     st = st + incr
15    st = st + group*n - it*incr

      return
      end
