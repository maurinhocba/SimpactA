 SUBROUTINE surf25 ( heade, nelem)
 !******************************************************************
 !
 !*** get the surface definition from the element set
 !
 !******************************************************************
 USE surf_db
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem         !number of elements in the set
 TYPE (ele25), POINTER :: heade    !pointer to first element

 INTEGER (kind=4) :: ie
 TYPE (ele25), POINTER :: elem
 TYPE (srf_seg), POINTER  :: seg

 !Initialize empty list
 IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
 NULLIFY(surfa%head,surfa%tail)

 elem => heade      !point to first element in the set

 DO ie=1,nelem      !for each element in the set
   CALL new_seg(seg)   !get memory for a surface segment
   seg%nodes(1:4) = elem%lnods(1:4) !transfer nodal connectivities
   CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list

   elem => elem%next      !point to next element
 END DO
 surfa%nelem = nelem      !number of segments in surface
 RETURN
 END SUBROUTINE surf25
