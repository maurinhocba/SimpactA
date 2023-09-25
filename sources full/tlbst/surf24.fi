 SUBROUTINE surf24 ( heade, nelem)
 !******************************************************************
 !
 !*** get the surface definition from the element set
 !
 !******************************************************************
 USE surf_db
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem         !number of elements in the set
 TYPE (ele24), POINTER :: heade    !pointer to first element

 INTEGER (kind=4) :: in,ie
 TYPE (ele24), POINTER :: elem
 TYPE (srf_seg), POINTER  :: seg

 !Initialize empty list
 IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
 NULLIFY(surfa%head,surfa%tail)

 elem => heade      !point to first element in the set

 DO ie=1,nelem      !for each element in the set
   CALL new_seg(seg)   !get memory for a surface segment
   DO in = 1,3      !transfer nodal connectivities
     seg%nodes(in) = elem%lnods(in)
   END DO
   seg%nodes(4) = seg%nodes(3)    !make fourth node equal to third
   !        seg%nodes(4) = 0
   !        WRITE(lures, '(i8,2x,4i8)') elem%numel,seg%nodes(1:3)
   CALL add_seg(seg, surfa%head, surfa%tail) !add to surface list

   elem => elem%next      !point to next element
 END DO
 surfa%nelem = nelem      !number of segments in surface
 RETURN
 END SUBROUTINE surf24
