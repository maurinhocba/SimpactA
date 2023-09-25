SUBROUTINE surf08( heade, nelem)
!******************************************************************
!
!*** get the element connectivities from the element set
!    used for contact  ONLY for 2-NODE elements presently
!
!******************************************************************
USE surf_db !, ONLY : srf_seg,surfa,delete_seg,add_seg
IMPLICIT NONE

  !--- dummy arguments
  INTEGER(kind=4),INTENT(IN):: nelem      !number of elements
  TYPE(ele08),POINTER:: heade   !pointer to first element in the set
  !TYPE(cont_srf),POINTER:: surfa   !INTENT(OUT) surface data
  !--- Local variables
  INTEGER(kind=4):: ie
  TYPE(srf_seg),POINTER:: seg
  TYPE(ele08),POINTER:: e

  !GENERATE SURFACE data_base
  !Initialize empty list
  IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
  NULLIFY(surfa%head,surfa%tail)
  e => heade        !point to first element
  DO ie=1,nelem     !for each element
    ! store segment in list
    CALL new_seg(seg)                         !allocate surface segment
    seg%nodes(1) = e%lnods(1)                 !first node
    seg%nodes(2) = e%lnods(2)                 !second node
    CALL add_seg(seg, surfa%head, surfa%tail)  !add segment to the d_b
    e => e%next
  END DO
  surfa%nelem = nelem                         !keep number of segments

RETURN
END SUBROUTINE surf08
