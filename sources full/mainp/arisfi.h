 SUBROUTINE arisfi(naris,x,nndpd,euler)
 !***********************************************************************
 !
 !     calculates dependant displacements
 !
 !***********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: naris,nndpd(3,naris)
 REAL (kind=8),INTENT(IN OUT) :: x(:,:)
 REAL (kind=8),POINTER :: euler(:,:)
 END SUBROUTINE arisfi
