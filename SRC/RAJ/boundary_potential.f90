      subroutine boundary_potential
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added do loops and indices to lines 551-553 in order to satisfy pgf90
      
!     Purpose:

!	 To calculate the potential on the mesh boundary due to the
!	 boundary corrective charges (ie the original boundary charges
!        plus the corrective charges).  The boundary corrective charges
!        are provided by poiss in bdy_hor, bdy_sn and bdy_we as bordered
!        sine-sine transforms.  The edge charges are themselves bordered
!        sine transforms along the edges.

!     NOTE 1:  The results returned in bdy_horiz etc and to the mesh are
!              the actual values on the boundary planes, not the sums and
!              differences used in calculations in most of the subroutine.

!     NOTE 2:  The various Fourier transform routines produce transforms
!              in scrambled order with values for even indices preceding
!              those for odd indices.  In one dimension, sine transform
!              of length (N+1), N even, pick up values from the position
!              range (1, N-1) and the results are returned to the same
!              values, leaving the index order of the results as:
!              (0, N/2, N/4, 3N/4, N/8, 7N/8, 3N/8, 5N/8, ..., N).  Cosine
!              transforms are generated with index order:
!              (0, N, N/2, N/4, 3N/4, .....).  For consistency in
!              processing the Greens function, these must be coerced into
!              the same order as the sine transforms.  This is done by
!              subroutine `greens_function' for the Greens function values
!              and in this routine for the CS, SC and CC transforms
!              generated from the SS transforms provided by subroutine
!             `solve_poisson'.

      use constants
      
      use greens_fn
      
      use interfac, local => boundary_potential
      
      use mesh_data
      
      use mesh_specification

      use workspace
      
      use verify
      
!     Type definitions.

      type face
         double precision, allocatable :: sc(:,:), cs(:,:), cc(:,:)
	 double precision, pointer     :: ss(:,:)
      end type face
      
      type bdy
         type(face)                    :: horizontal, southnorth, westeast
      end type bdy

!     Local variables for general working.
      
      integer                          :: i, j, k, p1, p2, p3
      integer                          :: ch1, ch2, ch3
      integer, save                    :: number = 0
      double precision                 :: contribution, normalise
      double precision, allocatable    :: ptr_we(:), contrib(:)
      double precision, allocatable    :: cc(:,:), cs(:,:), sc(:,:)
      double precision, pointer        :: green_value(:)
      double precision, pointer        :: pt(:), pt1(:)
      type(bdy), target                :: source, values
      character*5                      :: part

!     Local integers used as pointers to the Greens function and the boundary
!     areas.  The pointers ptr_horiz and ptr_sn are used to pick out
!     the even or odd terms as required and change from 1 to 2 on encountering
!     the first transform value with odd index.  The pointer array ptr_we
!     is set to control use of axis 3 transforms in the same way.

      integer                       :: ptr_horiz, ptr_sn
      integer                       :: ip, ip1, ip2, ip3, or2
      integer                       :: index(4)
      
!     *********************************************************************

!     Part 1:  Initial setting up and acquisition of the boundary sine-sine
!              transforms.  

!     *********************************************************************
      
!     Acquire working memory.  The transforms for the source and the
!     results are held in the structures `source' and `values' with the
!     obvious layout.
      
      allocate(source%horizontal%sc(n23,2), source%horizontal%cs(n23,2),  &
&              source%horizontal%cc(n23,2), source%southnorth%sc(n13,2),  &
&              source%southnorth%cs(n13,2), source%southnorth%cc(n13,2),  &
&              source%westeast%sc(n12,2),   source%westeast%cs(n12,2),    &
&              source%westeast%cc(n12,2))
      allocate(values%horizontal%sc(n23,2), values%horizontal%cs(n23,2),  &
&              values%horizontal%cc(n23,2), values%southnorth%sc(n13,2),  &
&              values%southnorth%cs(n13,2), values%southnorth%cc(n13,2),  &
&              values%westeast%sc(n12,2),   values%westeast%cs(n12,2),    &
&              values%westeast%cc(n12,2))

!     Set up the SS pointers to the CS areas in `source'.  These are used
!     as temporary work space in part 2.

      source%horizontal%ss => source%horizontal%cs
      source%southnorth%ss => source%southnorth%cs
      source%westeast%ss   => source%westeast%cs
      
!     Calculate sums and differences of the sine-sine transforms of the
!     boundary data in the SS areas in the source arrays.  These transforms
!     are bordered by edge charges which are required in calculating CS, SC
!     and CC transforms in part 3 but are not required elsewhere.

      if(symmetric) then
      
!$OMP  PARALLEL DO
      
         do i = 1, n23
            source%horizontal%ss(i,1) = bdy_hor(i,1) + bdy_hor(i,2)
            source%horizontal%ss(i,2) = bdy_hor(i,1) - bdy_hor(i,2)
            source%southnorth%ss(i,1) = bdy_sn(i,1)  + bdy_sn(i,2)
            source%southnorth%ss(i,2) = bdy_sn(i,1)  - bdy_sn(i,2)
            source%westeast%ss(i,1)   = bdy_we(i,1)  + bdy_we(i,2)
            source%westeast%ss(i,2)   = bdy_we(i,1)  - bdy_we(i,2)
         end do

!$OMP  END PARALLEL DO
      
      else

!$OMP  PARALLEL DO
	 
	 do i = 1, n23
	    source%horizontal%ss(i,1) = bdy_hor(i,1) + bdy_hor(i,2)
            source%horizontal%ss(i,2) = bdy_hor(i,1) - bdy_hor(i,2)
	 end do

!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO
	 
	 do i = 1, n13
            source%southnorth%ss(i,1) = bdy_sn(i,1)  + bdy_sn(i,2)
            source%southnorth%ss(i,2) = bdy_sn(i,1)  - bdy_sn(i,2)
	 end do
	 
!$OMP  PARALLEL DO

	 do i = 1, n12
            source%westeast%ss(i,1)   = bdy_we(i,1)  + bdy_we(i,2)
            source%westeast%ss(i,2)   = bdy_we(i,1)  - bdy_we(i,2)
	 end do

!$OMP  END PARALLEL DO

      end if


!     Clear the horizontal edges in the south and north planes.  These
!     contributions are duplicated in and picked up from the horizontal
!     planes.

      do i = 1, 2
         source%southnorth%ss(1:n3,i)	      = 0.0d0
         source%southnorth%ss(n13-n3+1:n13,i) = 0.0d0
      end do

!     Clear the horizontal edges in the west and east boundaries.
!     Contributions for horizontal edges are duplicated in and picked up
!     from the horizontal planes and those for vertical edges in and from
!     the southnorth planes.
      
      do i = 1, 2
         source%westeast%ss(1:n2,i)	    = 0.0d0
         source%westeast%ss(n12-n2+1:n12,i) = 0.0d0
         source%westeast%ss(1:n12:n2,i)     = 0.0d0
         source%westeast%ss(n2:n12:n2,i)    = 0.0d0
      end do
      
!     Copy sums and differences to the SC areas for use in part 3 and later.

      source%horizontal%sc      = source%horizontal%ss
      source%southnorth%sc      = source%southnorth%ss
      source%westeast%sc        = source%westeast%ss
      
!     *********************************************************************

!     Part 2:  Calculate the CSS, SCS and SSC terms which are parallel
!     terms for each opposed pair of boundaries.

!     *********************************************************************

!     Clear the edges as necessary to avoid over-weighting the edge charges.
!     Note that the edges for the west and east planes are already clear.

      do i = 1, 2

!        Clear the vertical edges of the southnorth planes.  The horizontal
!        edges have already been cleared.

         source%southnorth%ss(1:n13:n3,i)  = 0.0d0
         source%southnorth%ss(n3:n13:n3,i) = 0.0d0
      
!        Clear the edges of the horizontal planes.
      
	 source%horizontal%ss(1:n3,i)	      = 0.0d0
	 source%horizontal%ss(n23-n3+1:n23,i) = 0.0d0
	 source%horizontal%ss(1:n23:n3,i)     = 0.0d0
	 source%horizontal%ss(n3:n23:n3,i)    = 0.0d0

      end do
      
!     Calculate the parallel contributions to the sums and differences
!     of the boundary potential.  We calculate the sums in planes
!     corresponding to second index 1 in bdy_hor etc and differences 
!     in planes corresponding to second index 2.  Sums and differences
!     of the sine transforms of the Greens function planes on the
!     boundary are provided by subroutine greens_function.

      bdy_hor(:,1) = source%horizontal%ss(:,1)*green_horiz(:,1)
      bdy_hor(:,2) = source%horizontal%ss(:,2)*green_horiz(:,2)
      bdy_sn(:,1)  = source%southnorth%ss(:,1)*green_sn(:,1)
      bdy_sn(:,2)  = source%southnorth%ss(:,2)*green_sn(:,2)
      bdy_we(:,1)  = source%westeast%ss(:,1)*green_we(:,1)
      bdy_we(:,2)  = source%westeast%ss(:,2)*green_we(:,2)
	 
!     The parallel pair calculation is now complete and the ss variables
!     are no longer required.

      nullify(source%horizontal%ss, source%southnorth%ss, source%westeast%ss)
     
!     ***********************************************************************

!     Part 3:  Set up the sc, cs and cc transforms of the boundary charges,
!              made up of boundary charges present on the original mesh and
!              corrective charges calculated in subroutine solve_poisson.
!              This routine uses the convention that the last indicator of
!              the pair sc or cs refers to the store direction.  After
!              generating each cosine transform the order of indices is
!              coerced to match the sine transforms.

!     ***********************************************************************

!     Cycle over the sum and difference terms.

p_m3: do i = 1, 2

!        Calculate the CS transforms from the SS planes held in the
!        SC variables as bordered SS transforms.
      
         pt => source%horizontal%cs(:,i)
	 pt = source%horizontal%sc(:,i)
	 call vfthil(n2, n3, pt, .true.)
	 
         pt => source%southnorth%cs(:,i)
	 pt = source%southnorth%sc(:,i)
	 call vfthil(n1, n3, pt, .true.)
	 
         pt => source%westeast%cs(:,i)
	 pt = source%westeast%sc(:,i)
         call vfthil(n1, n2, pt, .true.)
	 
!        Coerce CS planes to the required order.
	 
	 pt => source%horizontal%cs(n3+1:n23,i); pt = cshift(pt, n3)
	 pt => source%southnorth%cs(n3+1:n13,i); pt = cshift(pt, n3)
	 pt => source%westeast%cs(n2+1:n12,i);   pt = cshift(pt, n2)
	 
!        Calculate the CC transforms for the horizontal planes, coercing
!        to the required order after transformation.
      
         pt => source%horizontal%cc(:,i)
	 pt =  source%horizontal%cs(:,i)
         call transp(pt, n2, n3)
         call vfthil(n3, n2, pt, .true.)
	 pt1 => source%horizontal%cc(n2+1:n23,i); pt1 = cshift(pt1, n2)	 
         call transp(pt, n3, n2)
	 
!        Calculate the CC transforms for the south and north planes
      
         pt => source%southnorth%cc(:,i)
	 pt =  source%southnorth%cs(:,i)
         call transp(pt, n1, n3) 
         call vfthil(n3, n1, pt, .true.)
	 pt1 => source%southnorth%cc(n1+1:n13,i); pt1 = cshift(pt1, n1)
         call transp(pt, n3, n1)
      
!        Calculate the CC transforms for the west and east planes.

         pt => source%westeast%cc(:,i)
	 pt =  source%westeast%cs(:,i)
         call transp(pt, n1, n2)
         call vfthil(n2, n1, pt, .true.)
	 pt1 => source%westeast%cc(n1+1:n12,i); pt1 = cshift(pt1, n1)
         call transp(pt, n2, n1)
      
!        Calculate the SC transforms for the horizontal planes.
      
         pt => source%horizontal%sc(:,i)
	 call transp(pt, n2, n3)
         call vfthil(n3, n2, pt, .true.)
	 pt1 => source%horizontal%sc(n2+1:n23,i); pt1 = cshift(pt1, n2)
         call transp(pt, n3, n2)
      
!        Calculate the SC transforms for the south and north planes.
      
         pt => source%southnorth%sc(:,i)
	 call transp(pt, n1, n3)
         call vfthil(n3, n1, pt, .true.)
	 pt1 => source%southnorth%sc(n1+1:n13,i); pt1 = cshift(pt1, n1)
         call transp(pt, n3, n1)
      
!        Calculate the SC transforms for the west and east planes.
      
         pt => source%westeast%sc(:,i)
	 call transp(pt, n1, n2)
         call vfthil(n2, n1, pt, .true.)
	 pt1 => source%westeast%sc(n1+1:n12,i); pt1 = cshift(pt1, n1)
         call transp(pt, n2, n1)
      
!        Clear the bordering charges for the sine transformed directions.
!        Subsequent calculations get these contributions come from the
!        cosine transforms in those directions.

!        SC transforms.  The sine transforms are perpendicular to the
!        store direction so the edges for clearing are parallel to it.
      
         source%horizontal%sc(1:n3,i)         = 0.0d0
         source%horizontal%sc(n23-n3+1:n23,i) = 0.0d0
	 source%southnorth%sc(1:n3,i)         = 0.0d0
	 source%southnorth%sc(n13-n3+1:n13,i) = 0.0d0
	 source%westeast%sc(1:n2,i)           = 0.0d0
	 source%westeast%sc(n12-n2+1:n12,i)   = 0.0d0
	 
!        CS transforms.  The sine transforms are parallel to the store
!        direction so the edges for clearing are perpendicular to it.
	 
	 source%horizontal%cs(1:n23:n3,i)     = 0.0d0
	 source%horizontal%cs(n3:n23:n3,i)    = 0.0d0
	 source%southnorth%cs(1:n13:n3,i)     = 0.0d0
	 source%southnorth%cs(n3:n13:n3,i)    = 0.0d0
	 source%westeast%cs(1:n12:n2,i)       = 0.0d0
	 source%westeast%cs(n2:n12:n2,i)      = 0.0d0
      
      end do                                                           p_m3
      
!     *********************************************************************

!     Part 4: Calculate contributions to the various potential transforms.
!             This requires a scan over the cosine transformed Greens
!             function to evaluate the cross terms.  We calculate the
!             CCC, CCS, CSC and SCC terms at each point of the mesh and
!             find it convenient to calculate parallel boundary
!             contributions from these transforms during the scan.
!             Values of the SSS transform of interior charges are not used
!             or modified at this stage.  The corrective contributions to
!             the boundary charges were calculated by solve_poisson and
!             combined with existing boundary charges at that stage.

!     *********************************************************************
      
!     Set up the initial bottom/top selector to select sums or differences
!     of boundary terms for the scan over the mesh.  Note that these are
!     used for selecting in a cosine transformed direction in which sums
!     contribute to even index values and differences to odd index values.
!     This contrasts with the case of sine transformed directions in 
!     subroutines solve_poisson and solve_laplace where sums contribute to
!     odd index values and differences to even index values.
      
!     Clear the result areas.
      
      values%horizontal%sc = 0.0d0
      values%horizontal%cs = 0.0d0
      values%horizontal%cc = 0.0d0
      values%southnorth%sc = 0.0d0
      values%southnorth%cs = 0.0d0
      values%southnorth%cc = 0.0d0
      values%westeast%sc   = 0.0d0
      values%westeast%cs   = 0.0d0
      values%westeast%cc   = 0.0d0
      
!     Set up the pointer change values.  Recall that values are laid out
!     in the order 0, N/2, N/4, 3N/4, ... , N where N stands for one of
!     n1, n2, n3.
      
      ch1 = (n1 + 1)/2
      ch2 = (n2 + 1)/2
      ch3 = (n3 + 1)/2
      p3  = (n3 - 1)/2
      
!     Acquire memory for the west-east pointers and set it up.

      allocate(ptr_we(n3))
      ptr_we = 1
      ptr_we(ch3:n3-1) = 2
      
!$OMP  PARALLEL                                                          &
!$OMP  DEFAULT(SHARED)                                                   &
!$OMP  PRIVATE(cc, cs, sc, contrib, green_value, ip, ip1, ip2, ip3, j,   &
!$OMP          or2, p1, p2, ptr_horiz, ptr_sn)

!     Acquire memory for the contributions vector.

      allocate(contrib(n3))

!     Acquire memory for the working areas for contributions to values on
!     the horizontal planes and clear them.  Use of these working areas
!     imposes a slight penalty in serial mode but is essential in parallel
!     mode.

      allocate(cc(n23,2), cs(n23,2), sc(n23,2))
      sc	  = 0.0d0
      cs	  = 0.0d0
      cc	  = 0.0d0

!$OMP  DO

!     Cycle over the planes of the mesh.

pln:  do i = 1, n1

!        Set the initial value of the Greens function pointer.

         ip = (i-1)*n23 + 1
	 
!	 Set the initial local origin pointer for the south and north planes.
	 
	 or2 = (i-1)*n3 + 1
	 
!	 Set the pointer for the westeast planes.
	 
	 ip3 = (i-1)*n2 + 1
	 
!	 Set the top/bottom pointer for the horizontal planes.
	 
	 ptr_horiz = 1
	 if((i.ge.ch1).and.(i.ne.n1)) ptr_horiz = 2

!        Cycle over the lines in a plane.

!        Set the initial value for the south-north sum or difference selector.
	 
	 ptr_sn   = 1
	 
!        Set the initial pointer for the horizontal planes.
	 
	 ip1 = 1
	 
!        Cycle over lines.
	 
lin: 	 do j = 1, n2
	 
!	    Set the pointer for the south and north planes to the local origin.
	    
	    ip2   = or2
	    
!	    Change the even/odd south/north pointer if necessary.

            if(j.eq.ch2) ptr_sn = 2
	    if(j.eq.n2)  ptr_sn = 1
	    
!           Pick up local Greens function values.

            green_value => green(ip:ip+n3-1)
	    
!           Contributions from the CCC transform on the boundary.  These
!           collect contributions from and accumulate results to the CC
!           transforms on all the boundaries.  Values for these and all
!           other transforms are picked up from the even or odd boundary
!           terms in the area `source' using the selectors ptr_horiz,
!           ptr_sn, ptr_we.  The results are returned to the
!           corresponding regions in the area `values'.  Values of the
!           sums(differences) are in the planes with second index 1(2)
!           in the sources and the results.

            contrib = green_value*                                       &
&		     (source%horizontal%cc(ip1:ip1+n3-1,ptr_horiz) +	 &
&		      source%southnorth%cc(ip2:ip2+n3-1,ptr_sn)    +	 &
&		      source%westeast%cc(ip3,ptr_we))
 
            p1 = ip1 + n3 - 1
	    cc(ip1:p1,ptr_horiz) = cc(ip1:p1,ptr_horiz) + contrib

            p2 = ip2 + n3 - 1
	    values%southnorth%cc(ip2:p2,ptr_sn) =                        &
&	          values%southnorth%cc(ip2:p2,ptr_sn)    + contrib
	    
	    values%westeast%cc(ip3,1)          =                         &
&	          values%westeast%cc(ip3,1) + contrib(n3) +              &
&                                             sum(contrib(1:p3))
            values%westeast%cc(ip3,2)          =                         &
&	          values%westeast%cc(ip3,2) + sum(contrib(p3+1:n3-1))

!	    Contributions from the CCS transform.  These draw values
!           from and return results to CS transforms on the horizontal,
!           south and north boundaries.

            contrib = green_value*                                       &
&	             (source%horizontal%cs(ip1:ip1+n3-1,ptr_horiz) +     &
&                     source%southnorth%cs(ip2:ip2+n3-1,ptr_sn))

            cs(ip1:p1,ptr_horiz) = cs(ip1:p1,ptr_horiz) + contrib

            values%southnorth%cs(ip2:p2,ptr_sn)     =                    &
&	          values%southnorth%cs(ip2:p2, ptr_sn)   + contrib
	       
!           Contributions from the CSC transform.  These interact with
!           the CS source and result planes on the west and east
!           boundaries and the SC planes on the horizontal boundaries.

            contrib = green_value*                                       &
&	             (source%horizontal%sc(ip1:p1,ptr_horiz) +           &
&                     source%westeast%cs(ip3,ptr_we))

            sc(ip1:p1,ptr_horiz) = sc(ip1:p1,ptr_horiz) + contrib

            values%westeast%cs(ip3,1)              =                     &
&	          values%westeast%cs(ip3,1) + contrib(n3) +              &
&                                             sum(contrib(1:p3))
            values%westeast%cs(ip3,2)              =                     &
&	          values%westeast%cs(ip3,2) + sum(contrib(p3+1:n3-1))

!           Contributions from the SCC transform.  These interact with
!           the SC source and result planes on the southnorth and
!           westeast boundaries.

            contrib = green_value*                                       &
&	             (source%southnorth%sc(ip2:p2,ptr_sn) +              &
&                     source%westeast%sc(ip3,ptr_we))

            values%southnorth%sc(ip2:p2,ptr_sn) =                        &
&	          values%southnorth%sc(ip2:p2,ptr_sn) + contrib

            values%westeast%sc(ip3,1)           =                        &
&	          values%westeast%sc(ip3,1) + contrib(n3) +              &
&                                             sum(contrib(1:p3))
            values%westeast%sc(ip3,2)           =                        &
&	          values%westeast%sc(ip3,2)  + sum(contrib(p3+1:n3-1))

!           Advance the pointers for the horizontal and south and north
!           planes.

            ip1 = ip1 + n3
	    ip2 = ip2 + n3

!	    Advance the pointer for the west and east planes.
	    
	    ip3 = ip3 + 1
	    
!	    Advance the pointer for the Greens function.
	    
	    ip  = ip  + n3
	    
	 end do			                                        lin
	 
!        Omit accumulation of contributions if not at end of loop.

      end do                                                            pln
      
	 
!$OMP  END DO

!     Accumulate the horizontal plane contributions to the result area.   

!  do loops added by JAS, without which these lines would not compile in pgf90!
      do i = 1, 2
      do j = 1, n23
      values%horizontal%cc(j,i) = values%horizontal%cc(j,i) + cc(j,i)
      values%horizontal%cs(j,i) = values%horizontal%cs(j,i) + cs(j,i)
      values%horizontal%sc(j,i) = values%horizontal%sc(j,i) + sc(j,i)
      end do
      end do

!!$OMP CRITICAL(coscos)
!      values%horizontal%cc = values%horizontal%cc + cc
!!$OMP END CRITICAL(coscos)
!!$OMP CRITICAL(cossine)
!      values%horizontal%cs = values%horizontal%cs + cs
!!$OMP CRITICAL(sinecos)
!      values%horizontal%sc = values%horizontal%sc + sc
!!$OMP END CRITICAL(sinecos)

!     Release the horizontal contribution working areas.	          

      deallocate(cc, cs, sc, contrib) 	          

!$OMP END PARALLEL

!     Release the pointer array for the west and east boundaries.

      deallocate(ptr_we)

!     ***************************************************************

!     Part 5:  Consolidate the results into bordered SS transforms.
!              Border terms have contributions from all transforms
!              and these must not be cleared.

!     ***************************************************************

!     Set up pointers to the areas in `values' to receive the SS
!     transforms.  These contain CS transforms at this stage and the
!     first step is to convert them to bordered SS transforms.

      values%horizontal%ss => values%horizontal%cs
      values%southnorth%ss => values%southnorth%cs
      values%westeast%ss   => values%westeast%cs
      
!     Cycle over the even and odd parts of the transforms.

      do i = 1, 2

!        Coerce CS transforms to the standard form for the cosine transform
!        as required for vfthil.  This order for a transform of length
!        (N+1) is: 0, N, N/2, N/4, 3N/4, ...

	 pt => values%horizontal%ss(n3+1:n23,i); pt = cshift(pt, -n3)
	 pt => values%southnorth%ss(n3+1:n13,i); pt = cshift(pt, -n3)
	 pt => values%westeast%ss(n2+1:n12,i);   pt = cshift(pt, -n2)

!        CS transforms to bordered SS.
      
         call vfthil(n2, n3, values%horizontal%ss(:,i), .false.)
         call vfthil(n1, n3, values%southnorth%ss(:,i), .false.)
         call vfthil(n1, n2, values%westeast%ss(:,i), .false.)
	 
      end do
      
      do i = 1, 2
	 
!        Coerce CC transforms to standard form perpendicular to store
!        direction.

         pt => values%horizontal%cc(n3+1:n23,i); pt = cshift(pt, -n3)
	 pt => values%southnorth%cc(n3+1:n13,i); pt = cshift(pt, -n3)
	 pt => values%westeast%cc(n2+1:n12,i);   pt = cshift(pt, -n2)

!        CC transforms to bordered SC.
      
         call vfthil(n2, n3, values%horizontal%cc(:,i), .false.)
         call vfthil(n1, n3, values%southnorth%cc(:,i), .false.)
         call vfthil(n1, n2, values%westeast%cc(:,i), .false.)
	 
!	 Coerce SC transforms to standard form and convert to bordered SS.
	 
         pt => values%horizontal%cc(:,i)
	 call transp(pt, n2, n3)
	 pt1 => values%horizontal%cc(n2+1:n23,i); pt1 = cshift(pt1, -n2)
         call vfthil(n3, n2, pt, .false.)
         call transp(pt, n3, n2)

         pt => values%southnorth%cc(:,i)
	 call transp(pt, n1, n3)
	 pt1 => values%southnorth%cc(n1+1:n13,i); pt1 = cshift(pt1, -n1)
         call vfthil(n3, n1, pt, .false.)
         call transp(pt, n3, n1)

         pt => values%westeast%cc(:,i)
	 call transp(pt, n1, n2)
	 pt1 => values%westeast%cc(n1+1:n12,i); pt1 = cshift(pt1, -n1)
         call vfthil(n2, n1, pt, .false.)
         call transp(pt, n2, n1)
	 
      end do
      
!     Coerce SC transforms to standard form and convert to bordered SS.
	 
      do i = 1, 2
	 
	 pt => values%horizontal%sc(:,i)
	 call transp(pt, n2, n3)
	 pt1 => values%horizontal%sc(n2+1:n23,i); pt1 = cshift(pt1, -n2)
         call vfthil(n3, n2, pt, .false.)
         call transp(pt, n3, n2)

         pt => values%southnorth%sc(:,i)
	 call transp(pt, n1, n3)
	 pt1 => values%southnorth%sc(n1+1:n13,i); pt1 = cshift(pt1, -n1)
         call vfthil(n3, n1, pt, .false.)
         call transp(pt, n3, n1)

         pt => values%westeast%sc(:,i)
	 call transp(pt, n1, n2)
	 pt1 => values%westeast%sc(n1+1:n12,i); pt1 = cshift(pt1, -n1)
         call vfthil(n2, n1, pt, .false.)
         call transp(pt, n2, n1)
	 
      end do
      
!     Combine bordered SS transforms.
      
      do i = 1, 2
     
         values%horizontal%ss(:,i) = values%horizontal%ss(:,i) +          &
&        			     values%horizontal%sc(:,i) +          &
&                                    values%horizontal%cc(:,i)
      
         values%southnorth%ss(:,i) = values%southnorth%ss(:,i) +          &
&        			     values%southnorth%sc(:,i) +          &
&                                    values%southnorth%cc(:,i)
      
         values%westeast%ss(:,i)   = values%westeast%ss(:,i) +	   	  &
&        			     values%westeast%sc(:,i) +            &
&                                    values%westeast%cc(:,i)

      end do
      
!     Add these values to the result areas.  Also normalise them to remove
!     the factor 2 generated by adding/subtracting sums and differences.

!     Horizontal boundaries.

      values%horizontal%ss(:,1) = bdy_hor(:,1) + values%horizontal%ss(:,1)
      values%horizontal%ss(:,2) = bdy_hor(:,2) + values%horizontal%ss(:,2)
      bdy_hor(:,1)              = values%horizontal%ss(:,1) +             &
&                                 values%horizontal%ss(:,2)	   
      bdy_hor(:,2)              = values%horizontal%ss(:,1) -             &
&                                        values%horizontal%ss(:,2)

!     South and north boundaries.

      values%southnorth%ss(:,1) = bdy_sn(:,1)  + values%southnorth%ss(:,1) 
      values%southnorth%ss(:,2) = bdy_sn(:,2)  + values%southnorth%ss(:,2) 
      bdy_sn(:,1)               = values%southnorth%ss(:,1) +             &
&				  values%southnorth%ss(:,2)	           
      bdy_sn(:,2)               = values%southnorth%ss(:,1) -             &
&				  values%southnorth%ss(:,2)

!     West and east boundaries.

      values%westeast%ss(:,1)  = bdy_we(:,1)  + values%westeast%ss(:,1)
      values%westeast%ss(:,2)  = bdy_we(:,2)  + values%westeast%ss(:,2)
      bdy_we(:,1)              = values%westeast%ss(:,1) +	          &
&				 values%westeast%ss(:,2)
      bdy_we(:,2)              = values%westeast%ss(:,1) -	          &
&				 values%westeast%ss(:,2)

!     Normalise the results.

      normalise = recip32*real((n2-1)*(n3-1), kind = kin)
      bdy_hor   = normalise*bdy_hor 
      normalise = recip32*real((n1-1)*(n3-1), kind = kin)
      bdy_sn    = normalise*bdy_sn	    
      normalise = recip32*real((n1-1)*(n2-1), kind = kin)
      bdy_we    = normalise*bdy_we	    

!     ****************************************

!     Part 6:  return the results to the mesh.

!     ****************************************

!     Copy values on the west and east boundaries.

      w(1:n123-n3+1:n3)  = bdy_we(:,1)
      w(n3:n123:n3)      = bdy_we(:,2)

!     Copy values on the south and north boundaries.

      j = 0
      do i = 0, n13-n3, n3
         w(j+1:j+n3) = bdy_sn(i+1:i+n3,1)
	 j           = j + n23
	 w(j-n3+1:j) = bdy_sn(i+1:i+n3,2)
      end do

!     Copy values on the horizontal boundaries.

      w(1:n23)           = bdy_hor(:,1)
      w(n123-n23+1:n123) = bdy_hor(:,2)

!     Release working memory.

      nullify(values%horizontal%ss, values%southnorth%ss,                 &
&             values%westeast%ss, pt)
      
      deallocate(source%horizontal%sc, source%horizontal%cs,              &
&                source%horizontal%cc, source%southnorth%sc,	   	  &
&                source%southnorth%cs, source%southnorth%cc,	   	  &
&                source%westeast%sc,   source%westeast%cs,  	   	  &
&                source%westeast%cc)

      deallocate(values%horizontal%sc, values%horizontal%cs,              &
&                values%horizontal%cc, values%southnorth%sc,	   	  &
&                values%southnorth%cs, values%southnorth%cc,	   	  &
&                values%westeast%sc,   values%westeast%cs,  	   	  &
&                values%westeast%cc)

!     Finish.
      
      return
      end
