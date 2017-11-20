c radius arrays for Qsetup etc
      integer mradq
      parameter ( mradq = 100 )
      integer nradq
      real ake, bke, drq, qset
      common / qset / bke, ake, nradq, drq, qset( 6, mradq )
