 SUBROUTINE fixsld(x,euler)
 !***********************************************************************
 !
 !     updates configuration and factors of dependant displacements
 !
 !***********************************************************************
 IMPLICIT NONE
 REAL (kind=8),INTENT(IN OUT) :: x(:,:)
 REAL (kind=8), POINTER :: euler(:,:)
 END SUBROUTINE fixsld
