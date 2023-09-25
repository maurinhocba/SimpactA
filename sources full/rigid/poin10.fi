      SUBROUTINE poin10 (nelem, nnode, lnods, matno, flag)

      IMPLICIT NONE
      INTEGER (kind=4), INTENT(IN) :: nelem, nnode
      INTEGER (kind=4), POINTER :: lnods(:,:), matno(:)
      LOGICAL (kind=4) :: flag

      ALLOCATE ( lnods(nnode,nelem) )
      lnods = 0
      IF( flag )THEN
        ALLOCATE ( matno(nelem) )
        matno = 0
      ELSE 
        NULLIFY (matno)
      END IF

      RETURN
      END SUBROUTINE poin10
