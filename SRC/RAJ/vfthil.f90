      subroutine vfthil(n, no, x, stc)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Routine identifier				   xd

!     Purpose:

!	 To perform a hilbert transformation on a 2-dimensional
!	 mesh held in a 1-dimensional array.  The transform direction
!	 is orthogonal to the store direction.

!     Parameters

!     n    =  transform length.

!     no   =  number of transforms to be evaluated in parallel.

!     x    -  a real array holding the data.

!     stc  -  (type logical) - requests that the hilbert transform
!	      be used to change a sine transform into a cosine
!	      transform.

!     Note:  Vectors are assumed to be contiguous in memory, and so
!     orthogonal to the transform direction.

      use interfac, local => vfthil
      
      use mesh_specification
      
      use twiddle_factors
      
!     Local variables:

      logical                               :: stc
      integer                               :: n, no
      integer, save                         :: number = 0
      double precision                      :: recip, u, v, xx1, xx2
      double precision, target              :: x(n*no)
      double precision, allocatable         :: ws1(:), ws2(:), ws3(:), ws4(:)
      double precision, allocatable         :: ws5(:), ws6(:), ws7(:)
      double precision, allocatable, target :: work(:)
      double precision, allocatable         :: sum(:), diff(:)
      double precision, pointer             :: pt0(:), pt1(:), pt2(:), pt3(:)
      double precision, pointer             :: pt4(:), pt5(:)
      integer                               :: st0, st1, st2, st3, st4, st5, set
      integer, allocatable                  :: w(:)
      integer                               :: i, i0, i1, i2, i3, i4, iw, iw2
      integer                               :: j, jj, k, nb, nm1, m

!     Set initial values and check transform length.

      nm1 = n - 1
      recip = 1.0/float(nm1)
      nb = nm1/2

!     Acquire working memory.

      allocate(ws1(no), ws2(no), ws3(no), ws4(no), ws5(no), ws6(no))
      allocate(sum(no), diff(no), work(3*no*nb/2), w(3*nb/2))

!     Set up pointers for work areas in stack in blank common.

      do i = 0, 3*nb/2 - 1
         w(i + 1) = i*no
      end do

!     Select prologue according to transform required.

      if(.not.stc) then

!        Prepare cos-to-sine transformation.

         st2 = nb*no
         ws1 = x(st2+1:st2+no)
         ws2 = 0.0d0

!        Compute end results.

         st1 = - no
         do j = 1, nb
            st1 = st1 + no
            st2 = st2 + no
            ws1 = ws1 + x(st1+1:st1+no)
            ws2 = ws2 + x(st2+1:st2+no)
         end do
         sum   = recip*(ws1 + ws2)
         diff = recip*(ws1 - ws2)
      
      else

!        Prepare sine-to-cos transformation.

!        Preserve sum and difference of end values.

         st0  =  nm1*no
	 pt0  => x(st0+1:st0+no)
         sum  =  2.0d0*(x(1:no) + pt0)
	 diff =  2.0d0*(x(1:no) - pt0)

!        Initial pointer for twiddle factors.

         iw = nb
      
!        Cycle over data for adjustment.

         do j = 9, n, 4
            st1 = st0 - no
            st2 = st1 - no
            st3 = st2 - no
            u = fac(iw - 1)
            v = fac(iw)
            iw = iw - 2

!           Adjust first five sets of values in current set.

            i0 = st0
            st0 = st3 - no
            x(i0+1:i0+no)   = u*x(st1+1:st1+no)
            x(st1+1:st1+no) = u*x(st2+1:st2+no)
            x(st2+1:st2+no) = v*x(st3+1:st3+no)
            x(st3+1:st3+no) = v*x(st0+1:st0+no)

         end do

!        Adjust first five values in first set.

         st1 = st0 - no
         st2 = st1 - no
         st3 = st2 - no
         u = fac(2)
         x(st0+1:st0+no) = u*x(st1+1:st1+no)
         x(st1+1:st1+no) = u*x(st2+1:st2+no)
         x(st2+1:st2+no) =   x(st3+1:st3+no)
         x(st3+1:st3+no) =   0.0d0
         x(1:no) = 0.0d0

      end if

!     Evaluation of s-alphas.

!     Calculate t-alphas from odd data

      st1 = no*(nb + 1)
      do jj = 1, nb
         work(w(jj)+1:w(jj)+no) = x(st1+1:st1+no)
      st1 = st1 + no
      end do

!     Initial number of sets for folding

      set = 1

!     Cycle over folds.

fld1: do

!        Advance to next fold

         nb = nb/4
         st1 = 1
         if((16*set + 1).gt.n) exit                                    fld1
         iw = nb + 1
         iw2 = iw + nb

!        Cycle over sets for folding

         do j = 1, nb

!           Step pointers

            st2 = st1 + set
            st3 = st2 + set
            st4 = st3 + set

!           Find twiddle factors

            u = fac(iw)
            xx1 = fac(iw2)
            xx2 = fac(iw2 + 1)
            iw = iw + 1
            iw2 = iw2 + 2

!           Fold current set

            do k = 1, set
               m = w(st1); pt1 => work(m+1:m+no)
               m = w(st2); pt2 => work(m+1:m+no)
               m = w(st3); pt3 => work(m+1:m+no)
               m = w(st4); pt4 => work(m+1:m+no)
               ws2 = xx1*(pt1 - pt2)
               ws4 =	   pt1 + pt2
               ws3 =  xx2*(pt3 - pt4)
               ws5 =	   pt3 + pt4
               pt1 =	u*(ws2 - ws3)
               pt3 =	   ws2 + ws3
               pt2 =	u*(ws4 - ws5)
               pt4 =	   ws4 + ws5
               st1 = st1 + 1
               st2 = st2 + 1
               st3 = st3 + 1
               st4 = st4 + 1
            end do
            st1 = st4

         end do

!        Increase set length and advance to next fold

         set = 4*set
      
      end do                                                           fld1

!     Check for single fold.

         if((8*set + 1).le.n) then

!        Perform single fold.

         st2 = set + 1
         u = fac(3)
         do j = 1, 2
            do k = 1, set			  
               m   = w(st1); pt1 => work(m+1:m+no)
               m   = w(st2); pt2 => work(m+1:m+no)
               ws2 = u*(pt1 - pt2)
               pt2 =	pt1 + pt2
               pt1 =	ws2
               st1 = st1 + 1			  	    
               st2 = st2 + 1
            end do
            st1 = st2				  	  
            st2 = st1 + set			  	  
            u = fac(4)  			  	  
         end do
      
      end if
  
!     Initialise counts for unfolding.

      set = nm1/4
      nb = 1

!     Generate lowest level t-alphas for unfolding.

      st1 = 1
      st2 = set + 1
      st3 = set + st2
      u = fac(2)
      do j = 1, set
         k   = w(st1); pt1 => work(k+1:k+no)
         k   = w(st2); pt2 => work(k+1:k+no)
         k   = w(st3); pt3 => work(k+1:k+no)
         pt3 = u*(pt2 - pt1)
         pt2 = pt2 + pt1
         pt1 = - pt3
         st1 = st1 + 1
         st2 = st2 + 1
         st3 = st3 + 1
      end do

!     Initial pointers for unfolding process.

      st0 = nm1/2 + 1
      set = set/2

!     Cycle over folding stages.

fld2: do

!     Initial pointers for current unfolding stage.

         st1 = 1
         st2 = set + 1
         st3 = set + st2
         st4 = set + st3
         st5 = st0 - set

!        Unfold first quartet.

         do j = 1, set
            m	= w(st1); pt1 => work(m+1:m+no)
            m	= w(st2); pt2 => work(m+1:m+no)
            m	= w(st3); pt3 => work(m+1:m+no)
            m	= w(st4); pt4 => work(m+1:m+no)
            m	= w(st5); pt5 => work(m+1:m+no)
            ws2 = pt1 + pt3
            ws3 = pt2 + pt4
            pt5 = pt1 - pt3
            pt1 = ws2 + ws3
            pt2 = ws2 - ws3
            st1 = st1 + 1
            st2 = st2 + 1
            st3 = st3 + 1
            st4 = st4 + 1
            st5 = st5 + 1
         end do

!        Adjust zero and pack odd section at end.

         st0 = st0 - set - (set + 1)/2
         st2 = st0
         do j = 1, set, 2
            m = w(st1); pt1 => work(m+1:m+no)
            m = w(st2); pt2 => work(m+1:m+no)
            pt2 = pt1
            st1 = st1 + 1
            st2 = st2 + 1
         end do

!        Next set requires single fold only.

         st1 = st2 + set
         st2 = st1 + set
         u = fac(2)
         do j = 1, set
            m = w(st1); pt1 => work(m+1:m+no)
            m = w(st2); pt2 => work(m+1:m+no)
            ws2 = u*pt2
            pt2 = pt1 - ws2
            pt1 = pt1 + ws2
            st1 = st1 + 1
            st2 = st2 + 1
         end do

!        Jump if initial unfolding.

         if(nb.ne.1) then

!           Initial pointers

            iw = 2
   	    iw2 = 3
   	    st1 = st2

   !	    Cycle over sets for unfolding.

   	    do jj = 2, nb

   !	       Find twiddle factors.

   	       u = fac(iw)
   	       xx1 = fac(iw2)
   	       xx2 = fac(iw2 + 1)
   	       iw = iw + 1
   	       iw2 = iw2 + 2

   !	       Initial pointers.

   	       st2 = st1 + set
   	       st3 = st2 + set
   	       st4 = st3 + set

   !	       Fold quartet.

   	       do j = 1, set
   		  m   = w(st1); pt1 => work(m+1:m+no)
   		  m   = w(st2); pt2 => work(m+1:m+no)
   		  m   = w(st3); pt3 => work(m+1:m+no)
   		  m   = w(st4); pt4 => work(m+1:m+no)
   		  ws2 = pt1 + u*pt3
   		  ws4 = pt1 - u*pt3
   		  ws3 = pt2 + u*pt4
   		  ws5 = pt2 - u*pt4
   		  pt1 = ws2 + xx1*ws3
   		  pt2 = ws2 - xx1*ws3
   		  pt3 = ws4 + xx2*ws5
   		  pt4 = ws4 - xx2*ws5
   		  st1 = st1 + 1
   		  st2 = st2 + 1
   		  st3 = st3 + 1
   		  st4 = st4 + 1
   	       end do

   !	       Pointer for next set.

   	       st1 = st4

   	    end do

   	 end if

!   	 Test for completion.

   	 if(set.le.2) exit                                             fld2
   	 set = set/4
   	 nb = 4*nb
      
      end do                                                           fld2

!     Single fold.

      if(set.eq.2) then
      
         st0 = st0 - 2
         m = w(st0);   pt0 => work(m+1:m+no)
         m = w(st0+1); pt5 => work(m+1:m+no)
         m = w(1);     pt1 => work(m+1:m+no)
         m = w(2);     pt2 => work(m+1:m+no)
         pt0 = pt1 + pt2
         pt5 = pt1 - pt2
         iw = 2
         st1 = st0 + 3
         do j = 9, n, 4
            u = fac(iw)
            m = w(st1);   pt1 => work(m+1:m+no)
            m = w(st1+1); pt2 => work(m+1:m+no)
            ws2 = u*pt2
            pt2 = pt1 - ws2
            pt1 = pt1 + ws2
            iw = iw + 1
            st1 = st1 + 2
         end do
      
      else

         st0 = st0 - 1
         m = w(st0); pt0 => work(m+1:m+no)
         m = w(1);   pt1 => work(m+1:m+no)
         pt0 = pt1
      
      end if

!     Return t-alphas to result area and pick up data for u-alphas.
!     This is stacked in the working area with blanks after the first
!     three elements, to allow for shifting down.

      nb = nm1/2 + 1
      st1 = no*nb - no
      st3 = 3*nm1/4
      do j = 1, nb
         m               = w(st3); pt3 => work(m+1:m+no)
         ws2             = x(st1+1:st1+no)
         x(st1+1:st1+no) = pt3
         pt3             = ws2
         st1             = st1 - no
         st3             = st3 - 1
      end do
      m = w(4);     pt0 => work(m+1:m+no)
      m = w(1);     pt1 => work(m+1:m+no)
      m = w(2);     pt2 => work(m+1:m+no)
      m = w(st3+1); pt3 => work(m+1:m+no)
      m = w(st3+2); pt4 => work(m+1:m+no)
      m = w(st3+3); pt5 => work(m+1:m+no)
      pt1 = pt3
      pt2 = pt4
      pt0 = pt5

!     Initialise folding of even data.

      set = 1
      nb = nm1/8
      st0 = 6*nb + 1
      
!     Cycle over folding stages.

fld3: do
      
         if(nb.lt.2) exit                                              fld3

!        Fold to next level.

         st1 = st0
         st2 = st0 + set
         st3 = st2 + set
         st4 = st3 + set
         iw = nb + 1
         iw2 = iw + iw

!        Cycle over sets for folding.

         do jj = 2, nb

!           no pointers.

            i3 = 3*set
            st1 = st1 - i3
            st2 = st2 - i3
            st3 = st3 - i3
            st4 = st4 - i3
            iw = iw - 1
            iw2 = iw2 - 2

!           Pick up twiddle factors.

            u = fac(iw)
            xx1 = fac(iw2 - 1)
            xx2 = fac(iw2)

!           Fold set.

            do j = 1, set

               st1 = st1 - 1
               st2 = st2 - 1
               st3 = st3 - 1
               st4 = st4 - 1
               m   = w(st1); pt1 => work(m+1:m+no)
               m   = w(st2); pt2 => work(m+1:m+no)
               m   = w(st3); pt3 => work(m+1:m+no)
               m   = w(st4); pt4 => work(m+1:m+no)
               ws2 = xx1*(pt1 - pt2)
               ws4 =	  pt1 + pt2
               ws3 = xx2*(pt3 - pt4)
               ws5 =	  pt3 + pt4
               pt1 = u*(ws2 - ws3)
               pt3 =	ws2 + ws3
               pt2 = u*(ws4 - ws5)
               pt4 =	ws4 + ws5

            end do
      
         end do

!        Fold top set and move back.

         i2 = set + set
         i4 = i2 + i2
         st1 = st1 - i2
         st2 = st1 + set
         st5 = i4 + 1
         st3 = st5 + i2
         st4 = st3 + set
         u = fac(2)
         do j = 1, set
            m = w(st1); pt1 => work(m+1:m+no)
            m = w(st2); pt2 => work(m+1:m+no)
            m = w(st3); pt3 => work(m+1:m+no)
            m = w(st4); pt4 => work(m+1:m+no)
            m = w(st5); pt5 => work(m+1:m+no)
            m = w(st5+1); pt0 => work(m+1:m+no)
            ws2 = pt1 - pt2
            pt4 = pt1 + pt2
            pt3 = u*ws2
            pt5 = 0.0d0
            pt0 = 0.0d0
            st1 = st1 + 1
            st2 = st2 + 1
            st3 = st3 + 1
            st4 = st4 + 1
            st5 = st5 + 2
         end do

!        Move next set.

         if(nb.ne.2) then
            st1 = st2
            st2 = st4 + i4
            do j = 1, i4
               m = w(st1); pt1 => work(m+1:m+no)
               m = w(st2); pt2 => work(m+1:m+no)
               m = w(st4); pt4 => work(m+1:m+no)
               pt2 = pt1
               pt4 = 0.0d0
               st1 = st1 + 1
               st2 = st2 + 1
               st4 = st4 + 1
            end do
         end if

!        Fold first quartet.

         st1 = 1
         st2 = st1 + set
         st3 = st2 + set
         st4 = st3 + set
         do j = 1, set
            m	= w(st1); pt1 => work(m+1:m+no)
            m	= w(st2); pt2 => work(m+1:m+no)
            m	= w(st3); pt3 => work(m+1:m+no)
            m	= w(st4); pt4 => work(m+1:m+no)
            ws4 = pt1 + pt2
            pt1 = pt1 - pt2
            pt3 = pt1
            pt2 = ws4 - pt4
            pt4 = ws4 + pt4
            st1 = st1 + 1
            st2 = st2 + 1
            st3 = st3 + 1
            st4 = st4 + 1
         end do
         set = i4

!        Test for folding complete.

         if(nb.le.4) exit                                              fld3

!        Advance to next fold.

         nb = nb/4
      
      end do                                                           fld3

!     Single fold if necessary.

      if((nb.eq.4).or.(nb.eq.1)) then	!  nb = 1 if transform length < 17.
      
         st1 = st0 - set
         st2 = st0
         u = fac(2)

!        Fold set.

         do j = 1, set
            st1 = st1 - 1
            st2 = st2 - 1
            m = w(st1); pt1 => work(m+1:m+no)
            m = w(st2); pt2 => work(m+1:m+no)
            ws2 = u*(pt1 - pt2)
            pt2 =    pt1 + pt2
            pt1 =    ws2
         end do

!        Clear zero multiplier locations.

         i2 = set + set
         st2 = i2 + 1
         do j = 1, set
            m = w(st2); pt2 => work(m+1:m+no)
            pt2 = 0.0d0
            st2 = st2 + 1
         end do

!        Fold first set.

         st1 = 1
         st2 = st1 + set
         do j = 1, set
            m	= w(st1); pt1 => work(m+1:m+no)
            m	= w(st2); pt2 => work(m+1:m+no)
            ws2 = pt1 - pt2
            pt2 = pt1 + pt2
            pt1 = ws2
            st1 = st1 + 1
            st2 = st2 + 1
         end do
         set = i2
      
      end if

!     Calculate lowest level u-alphas.

      st1 = 1
      st2 = set + 1
      st3 = set + st2
      do j = 1, set
         m   = w(st1); pt1 => work(m+1:m+no)
         m   = w(st2); pt2 => work(m+1:m+no)
         m   = w(st3); pt3 => work(m+1:m+no)
         ws1 = pt1 + pt2
         pt1 = pt2 - pt1
         pt2 = pt3 - ws1

         st1 = st1 + 1
         st2 = st2 + 1
         st3 = st3 + 1
      end do

!     Initialise unfolding.

      nb = 1
      iw = 1

!     Jump if one fold only

      if(set.gt.1) then

!        Cycle over folding stages.

fld4:    do

!           Initial pointers.

            st5 = set/2
            st4 = 1
            st3 = 1 - st5
            st2 = st3 - st5
            st1 = st2 - st5
            iw = nb
            iw2 = nb + nb - 1

!           Unfold two levels.

            i3 = set + st5
            do jj = 1, nb

!              Advance pointers.

               st1 = st1 + i3
               st2 = st2 + i3
               st3 = st3 + i3
               st4 = st4 + i3

!              Pick up twiddle factors.

               iw = iw + 1
               u = fac(iw)
               iw2 = iw2 + 2
               xx1 = fac(iw2)
               xx2 = fac(iw2 + 1)

!              Unfold quartet of sets.

               do j = 1, st5

            	  m = w(st1); pt1 => work(m+1:m+no)
            	  m = w(st2); pt2 => work(m+1:m+no)
            	  m = w(st3); pt3 => work(m+1:m+no)
            	  m = w(st4); pt4 => work(m+1:m+no)
            	  ws6 = u*pt3
            	  ws2 = pt1 + ws6
            	  ws3 = pt1 - ws6
            	  ws6 = u*pt4
            	  ws4 = pt2 + ws6
            	  ws5 = pt2 - ws6
            	  ws6 = xx1*ws4
            	  pt1 = ws2 + ws6
            	  pt2 = ws2 - ws6
            	  ws6 = xx2*ws5
            	  pt3 = ws3 + ws6
            	  pt4 = ws3 - ws6
            	  st1 = st1 + 1
            	  st2 = st2 + 1
            	  st3 = st3 + 1
            	  st4 = st4 + 1

               end do
      
            end do

!           Test for completion.

            if(set.le.4) exit                                          fld4

!           Advance to next fold.

            nb = 4*nb
            set = set/4
	    
	 end do                                                        fld4
	 
      end if
      
!     Single fold.

      if(set.eq.4) then

         st1 = - 1
         iw = nm1/4
         nb = iw
         do j = 1, nb
            st1 = st1 + 2
            iw  = iw + 1
            u	= fac(iw)
            m	= w(st1);   pt1 => work(m+1:m+no)
            m	= w(st1+1); pt2 => work(m+1:m+no)
            ws2 = u*pt2
            pt2 = pt1 - ws2
            pt1 = pt1 + ws2
         end do
      
      end if
      nb = nm1/2

!     Transfer back to data area.

      st1 = nb*no + no
      st2 = 1
      do j = 1, nb
         m = w(st2); pt2 => work(m+1:m+no)
         x(st1+1:st1+no) = pt2
         st1 = st1 + no
         st2 = st2 + 1
      end do

!     Epilogue

      if(.not.stc) then

!        Finish calculation of sine transform.

         st1 = no
         st2 = st1 + no
         u = - fac(2)
         st3 = st2 + no
         st0 = st3 + no
         x(1:no) = sum(1:no)
         x(st1+1:st1+no) = - x(st2+1:st2+no)
         x(st2+1:st2+no) = u*x(st3+1:st3+no)
         x(st3+1:st3+no) = u*x(st0+1:st0+no)


!        Cycle over rest of data.

         iw = 3
         do j = 9, n, 4
            st1 = st0 + no
            st2 = st1 + no
            st3 = st2 + no
            i0 = st0
            st0 = st3 + no
            u = - fac(iw + 1)
            v = - fac(iw)
            iw = iw + 2
            x(i0+1:i0+no) = u*x(st1+1:st1+no)
            x(st1+1:st1+no) = u*x(st2+1:st2+no)
            x(st2+1:st2+no) = v*x(st3+1:st3+no)
            x(st3+1:st3+no) = v*x(st0+1:st0+no)
         end do
         x(st0+1:st0+no) = diff
      
      else

!        Complete calculation of cosine transform.
 
         st1 = 0
         st2 = nm1*no/2 + no
         do j = 3, n, 2
            x(st1+1:st1+no) = x(st1+1:st1+no) + sum
            x(st2+1:st2+no) = x(st2+1:st2+no) + diff
            st1 = st1 + no
            st2 = st2 + no
         end do
         x(st1+1:st1+no) = x(st1+1:st1+no) + sum
 
      end if
      
!     Finish.
 
      return
      end
