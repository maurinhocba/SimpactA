 SUBROUTINE elmd04 (nelem, heade, taile, iwrit)
 !******************************************************************
 !
 !*** READ element DATA for 3-D tetrahedral element (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem, iwrit
 TYPE (ele04), POINTER :: heade,taile !pointers to first and last elements

 TYPE(section), POINTER :: sec        !pointer to a section
 LOGICAL found
 INTEGER (kind=4) :: i,ip,chnode,isec,osec,secty,mtype
 INTEGER (kind=4) :: mats(5) = (/ 1,5,6, 7,92 /)  !acceptable material types
 TYPE (ele04), POINTER :: elem      !auxiliar pointer


 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE(lures,"(//2X,'ELEMENT MAT           NODE NUMBERS ',  &
                           &    10x,'angles')",ERR=9999)

 osec = -1  !initializes

 DO
   ! loop over elements

   CALL listen('ELMD04')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   ALLOCATE (elem)           !reserve memory for data
   NULLIFY(elem%gausv)

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   IF( isec /= osec )THEN        !checks
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELMD04: Invalid Section Number used')
     secty = sec%secty        !check section type for element
     IF( secty /= 0 ) CALL runend('ELMD16: Invalid sect. type for SOLID')
     IF( sec%mabas >= 0 )THEN
       mtype = sec%mtbas%mtype  !check material type for element
       IF( .NOT. ANY(mtype == mats) )  &
         CALL runend('ELMD16: Invalid Mat. type for SOLAG')
     END IF
     osec = isec
   END IF

   elem%matno = ip        !store section position in array

   elem%lnods(1:4) = INT (param(3:6))    !connectivities
   elem%angle(1:3) = param(7:9)  ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i13,i8,4i7,3F10.4)",ERR=9999)elem%numel,isec, &
                                                   elem%lnods(1:4),elem%angle
   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   CALL add_ele04e (elem, heade, taile) !add element to list

 END DO
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd04
