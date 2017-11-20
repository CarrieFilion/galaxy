      subroutine fftsna(n, no, group, w, unity)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Routine identifier                                            tm

!     Purpose:

!        To perform a fourier sine analysis along one axis of a
!        2 or 3 dimensional mesh.

!     Parameters:

!     n      -  (type integer) gives the transform length.

!     no     -  (type integer) gives the total number of transforms required.

!     group  -  (type integer) gives the number of transforms to be
!               evaluated in parallel at each stage.

!     w      -  a double precision array containing the data.

!     unity  -  (type logical) specifies if the data for each individual
!               transform is laid out sequentially in memory.

      use constants

      use twiddle_factors
      
!     Parameter definitions:

      integer, intent(in)             :: n, no, group
      logical, intent(in)             :: unity
      double precision, intent(inout) :: w(n*no)

!     Local variables

      integer                       :: st, st1, step, p, q, pii, qi
      integer                       :: i, incr, ii, ip, iq, ipi, iqi, it, iw
      integer                       :: j, jump, lim, na, nb, ntot
      integer, save                 :: number = 0
      double precision              :: c, c1, s, u
      double precision, allocatable :: work1(:), work2(:), work3(:)
      logical                       :: rev

!     Acquire temporary memory.

      ntot = n*no
      allocate(work1(ntot), work2(ntot), work3(ntot))

!     Set initial values.

      if(unity) then
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
      na = n - 1
      lim = st + na*step

!     Fold real section.

 2    p = st + step
      jump = na*step
      q = st + jump - step
      nb = na/2
      if(nb.eq.1) go to 6
      do 3 j = 2, nb
      ip = p
      iq = q
      do 1 i = 1, it
      work1(i) = w(ip)
      w(ip) = work1(i) - w(iq)
      w(iq) = work1(i) + w(iq)
      ip = ip + incr
1     iq = iq + incr
      p = p + step
3     q = q - step
      iq = q
      do 5 i = 1, it
      w(iq) = 2.0d0*w(iq)
5     iq = iq + incr

!     No action required to convert real group to complex form.

!     Test scan complete.

 6    st1 = st + 2*jump
      if(st1.ge.lim) go to 4

!     Fold first complex group.

!     Initialise and set first terms.

      pii = st1
      q = pii + jump
      qi = q
      p = q + jump
      ipi = pii
      iqi = qi
      do 16 i = 1, it
      work1(i) = root2*w(iqi)
      w(iqi) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ipi = ipi + incr
16    iqi = iqi + incr
      if(nb.eq.1) go to 10

!     Fold intermediate terms.

      do 7 j = 2, nb
      pii = pii + step
      qi = qi - step
      p = p - step
      q = q + step
      ip = p
      iq = q
      ipi = pii
      iqi = qi
      do 7 i = 1, it
      work1(i) = recip_root2*(w(iq) - w(iqi))
      work2(i) = recip_root2*(w(iq) + w(iqi))
      w(iq) = work2(i) - w(ipi)
      w(ipi) = work2(i) + w(ipi)
      w(iqi) = w(ip) - work1(i)
      w(ip) = w(ip) + work1(i)
      ip = ip + incr
      iq = iq + incr
      ipi = ipi + incr
7     iqi = iqi + incr

!     Calculate end terms.

10    p = p - step
      pii = pii + step
      ip = p
      ipi = pii
      do 17 i = 1, it
      work1(i) = recip_root2*(w(ipi) + w(ip))
      w(ip) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ip = ip + incr
17    ipi = ipi + incr

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

14    pii = st1
      qi = pii + jump
      q = qi
      p = q + jump
      ipi = pii
      iqi = qi
      c1 = 1.0d0/c
      do 18 i = 1, it
      work1(i) = c1*w(iqi)
      w(iqi) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ipi = ipi + incr
18    iqi = iqi + incr
      if(nb.eq.1) go to 11

!     Fold intermediate terms.

      do 9 j = 2, nb
      pii = pii + step
      qi = qi - step
      q = q + step
      p = p - step
      ip = p
      iq = q
      ipi = pi
      iqi = qi
      do 9 i = 1, it
      work1(i) = c*w(iq) - s*w(iqi)
      work2(i) = s*w(iq) + c*w(iqi)
      w(iq) = work2(i) - w(ipi)
      w(ipi) = work2(i) + w(ipi)
      w(iqi) = w(ip) - work1(i)
      w(ip) = w(ip) + work1(i)
      ip = ip + incr
      iq = iq + incr
      ipi = ipi + incr
9     iqi = iqi + incr

!     Calculate end terms.

11    pii = pii + step
      p = p - step
      ip = p
      ipi = pii
      do 19 i = 1, it
      work1(i) = s*w(ip) + c*w(ipi)
      w(ip) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ip = ip + incr
19    ipi = ipi + incr

!     next group if not end of scan.

      st1 = st1 + 2*jump
      if(st1.lt.lim) go to 8

!     Jump back for next fold.

 4    na = nb
      if(na.ne.1) go to 2

!     Final reduction.

      st1 = st + 2*step
      iw = 1
      na = (n - 1)/2
      nb = 2*step
      do 12 j = 2, na
      iw = iw + 1
      p = st1 + step
      ip = p
      iq = st1
      c1 = 1.0d0/fac(iw)
      do 20 i = 1, it
      work1(i) = c1*w(ip)
      w(ip) = work1(i) - w(iq)
      w(iq) = work1(i) + w(iq)
      ip = ip + incr
20    iq = iq + incr
12    st1 = st1 + nb
      p = st + step
      ip = p
      do 21 i = 1, it
      w(ip) = 2.0d0*w(ip)
21    ip = ip + incr
15    st = st + group*n

!     Release temporary memory.

      deallocate(work1, work2, work3)

      return
      end
