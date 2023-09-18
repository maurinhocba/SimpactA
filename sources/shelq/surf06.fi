SUBROUTINE surf06(heade, nelem)
!******************************************************************
!
!*** get the surface definition from the element set
!
!******************************************************************
USE surf_db
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nelem    !number of elements in the set
  TYPE(ele06),POINTER:: heade
  !--- Local variables
  INTEGER(kind=4):: ie
  TYPE(ele06),POINTER:: e
  TYPE(srf_seg),POINTER:: seg

  !Initialize empty list
  IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
  NULLIFY(surfa%head,surfa%tail)

  e => heade
  DO ie=1,nelem      !for each element in the set
    CALL new_seg(seg)   !get memory for a surface segment
    seg%nodes = e%lnods
    CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list
    e => e%next
  END DO

  surfa%nelem = nelem      !number of segments in surface

RETURN
END SUBROUTINE surf06
