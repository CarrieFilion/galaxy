      module mesh_specification
      
!        This module holds general paramters which are needed in the
!        potential determination.
      
         integer	   :: formul = 3
         integer	   :: n1, n2, n3, m1, m2, m3
	 integer           :: cells(3)
	 integer           :: n23, n13, n12, n123
         double precision  :: w1, w2, w3, smooth = 0.0d0
	 double precision  :: h1 = 1.0d0, h2 = 1.0d0, h3 = 1.0d0
         logical	   :: symmetric, smoothing = .false.
	 
      end module
