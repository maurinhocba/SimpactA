 SUBROUTINE elmd04 (nelem, heade, taile, iwrit, ngaus)
 !******************************************************************
 !
 !*** READ element DATA for 3-D SOLID element  (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4) :: nelem, iwrit, ngaus
 TYPE (ele04), POINTER :: heade,taile !pointers to first and last elements
 ! local variables
 TYPE(section), POINTER :: sec        !pointer to a section
 LOGICAL found
 INTEGER (kind=4) :: i,ip,j,chnode,isec,osec,secty,mtype

 INTEGER (kind=4) :: mats(3) = (/ 1,5,6 /)  !acceptable material types
 TYPE (ele04), POINTER :: elem      !auxiliar pointer

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE(lures,"(//2X,'ELEMENT MAT           NODE NUMBERS ',&
                              &    4x,'angle')",ERR=9999)

 osec = -1  !initializes
 j = 2+nnode !last integer value to read
 DO
   ! loop over elements

   CALL listen('ELMD04')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele04e (elem)    !reserve memory for data

   elem%numel = INT (param(1))   !element label

   isec = INT (param(2))         !associated material

   IF( isec /= osec )THEN        !checks
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELMD04: Invalid Section Number used')
     secty = sec%secty        !check section type for element
     IF( secty /= 0 ) CALL runend('ELMD04: Invalid sect. type for SOLIDS')
     IF( sec%mabas >= 0 )THEN
       mtype = sec%mtbas%mtype  !check material type for element
       IF( .NOT. ANY(mtype == mats) )  &
         CALL runend('ELMD04: Invalid Mat. type for SOLIDS')
     END IF
     osec = isec
   END IF

   elem%matno = ip        !store section position in array

   elem%lnods = INT (param(3:j))  !connectivities
   elem%angle = param(j+1)        ! Euler angle between standart and Ortho
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i6,i3,8i6,F7.4)",ERR=9999) &
                  elem%numel,isec,elem%lnods,elem%angle
   !     original node labels are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO
   ! reserve memory for arrays and initializes
   ALLOCATE( elem%stint(6,ngaus))
   elem%stint = 0d0
   CALL add_ele04e (elem, heade, taile) !add element to list

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd04
