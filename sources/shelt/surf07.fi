SUBROUTINE surf07(heade, nelem,ngaus)
!******************************************************************
!
!*** get the surface definition from the element set
!
!******************************************************************
USE surf_db
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: nelem      !number of elements in the set
  INTEGER(kind=4):: ngaus      !number of integration points (flag)
  TYPE(ele07),POINTER:: heade
  !--- Local variables
  INTEGER(kind=4):: ie
  TYPE(srf_seg),POINTER:: seg
  TYPE(ele07),POINTER:: e

  !Initialize empty list
  IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
  NULLIFY(surfa%head,surfa%tail)

  e => heade
  DO ie=1,nelem      !for each element in the set

    CALL new_seg(seg)   !get memory for a surface segment
    IF( ngaus /= 1 )THEN
      seg%nodes = (/ e%lnods((/ 1,4,6 /)), 0 /)
      CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list

      CALL new_seg(seg)   !get memory for a surface segment
      seg%nodes = (/ e%lnods((/ 4,2,5 /)), 0 /)
      CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list

      CALL new_seg(seg)   !get memory for a surface segment
      seg%nodes = (/ e%lnods((/ 5,6,4 /)), 0 /)
      CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list

      CALL new_seg(seg)   !get memory for a surface segment
      seg%nodes = (/ e%lnods((/ 6,5,3 /)), 0 /)
    ELSE
      seg%nodes = (/ e%lnods(1:3), 0 /)
    END IF
    CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list

    e => e%next
  END DO

  IF( ngaus /= 1 )THEN
    surfa%nelem = 4*nelem    !number of segments in surface
  ELSE
    surfa%nelem = nelem      !number of segments in surface TLLL
  END IF

RETURN
END SUBROUTINE surf07
