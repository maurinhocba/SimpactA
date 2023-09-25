 SUBROUTINE elmd18 (nelem, heade, taile, iwrit, ngaus, shell)
 !******************************************************************
 !
 !*** READ element DATA for 3-D SOLID element  (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem, iwrit, ngaus
 TYPE (ele18), POINTER :: heade,taile !pointers to first and last elements
 TYPE(section), POINTER :: sec        !pointer to a section
 LOGICAL, INTENT(IN) :: shell         !assumed shear strain is active

 LOGICAL found
 INTEGER (kind=4) :: i,ip,chnode,isec,osec,secty,mtype

 INTEGER (kind=4) :: mats(6) = (/ 1,5,6,7,8,92 /)  !acceptable material types
 TYPE (ele18), POINTER :: elem


 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE(lures,"(//2X,'ELEMENT MAT           NODE NUMBERS ',  &
                           &    10x,'angles')",ERR=9999)

 osec = -1  !initializes
 DO
   ! loop over elements

   CALL listen('ELMD18')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele18e (elem)    !reserve memory for data
   ALLOCATE (elem%cartd(nnode,3,ngaus), elem%dvol(ngaus), elem%stint(6,ngaus))    !reserve memory for data

   elem%numel = INT (param(1))    !element label

   isec = INT (param(2))         !associated material

   IF( isec /= osec )THEN        !checks
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELMD18: Invalid Section Number used')
     secty = sec%secty        !check section type for element
     IF( secty /= 0 ) CALL runend('ELMD18: Invalid sect. type for SOLID')
     IF( sec%mabas >= 0 )THEN
       mtype = sec%mtbas%mtype  !check material type for element
       IF( .NOT. ANY(mtype == mats) )  &
         CALL runend('ELMD18: Invalid Mat. type for SOLAG')
     END IF
     osec = isec
   END IF

   elem%matno = ip        !store section position in array

   elem%lnods(1:nnode) = INT (param(3:10))  !connectivities
   elem%angle = param(11:13)   ! Euler angles between standart and Ortho
   IF( shell )THEN
     ALLOCATE(elem%nfdas(8,4),elem%jacin(2,2,ngaus))
     elem%angle(2:3) = 0d0
   ELSE
     NULLIFY( elem%nfdas,elem%jacin )
   END IF
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i10,i3,8i7,3F10.4)",ERR=9999) elem%numel,isec,  &
                               elem%lnods(1:nnode),elem%angle
   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   CALL add_ele18e (elem, heade, taile) !add element to list

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd18
