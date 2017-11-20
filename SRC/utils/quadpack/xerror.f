      SUBROUTINE XERROR (XMESS, NMESS, NERR, LEVEL)

c*********************************************************************72
c
cc XERROR replaces the SLATEC XERROR routine.
c
      CHARACTER*120 XMESS

      IF (LEVEL.GE.1) THEN
      IERR=I1MACH(3)
      WRITE(IERR,'(1X,A)') XMESS(1:NMESS)
      WRITE(IERR,'('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)')
     .      NERR,LEVEL
      END IF
      RETURN
      END

