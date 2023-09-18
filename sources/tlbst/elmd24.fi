 SUBROUTINE elmd24 (nelem, heade, taile, iwrit, shear)
 !******************************************************************
 !
 !*** READ element DATA for shell element BST-CST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER ::  nnode=3
 INTEGER (kind=4) :: nelem, iwrit, shear
 TYPE (ele24), POINTER :: heade,taile !pointers to first and last elements

 INTEGER (kind=4) :: i,ip,chnode,isec,secty,nlayr,osec,mtype
 INTEGER (kind=4) :: secs(1) = (/ 12 /)      !only acceptable section type
 TYPE (ele24), POINTER :: elem      !auxiliar pointer
 TYPE(section), POINTER :: sec
 LOGICAL :: found


 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//6X,'ELEMENT     MAT       NODE NUMBERS ', &
                             & 10x,'angle')",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('ELMD24')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele24e(elem)     !reserve memory for data

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))          !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('TLBST: Invalid Section Number used')
     secty = sec%secty       !check section type for element
     IF( ALL(secty /= secs) ) CALL runend('TLBST: Invalid sect. type for TLBST ')
     nlayr = sec%iprop(1)     !number of layers
     IF( sec%mabas >= 0 )THEN
       mtype = sec%mtbas%mtype  !material type
       IF( mtype /= 4 )CALL runend('ELTBSTL: Only Material 4 can be used')
     ELSE
       CALL runend('TBSTL: Only Material 4 can be used')
     END IF
     osec = isec
     IF( sec%iprop(4) == 1 .AND. .NOT.shear) THEN
       WRITE(*,"(' WARNING: transverse shear not computed for section',i3)")isec
       WRITE(55,"(' WARNING: transverse shear not computed for section',i3)")isec
     END IF
   END IF
   elem%matno = ip        !store section position in array

   elem%lnods(1:nnode) = INT (param(3:5))  !connectivities
   elem%lnods(4:6) = 0                     !initializes LSIDE
   elem%angle = param(6)   ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i13,i8,3i7,F14.4)",ERR=9999) elem%numel,isec, &
                                                    elem%lnods(1:3),elem%angle

   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   CALL add_ele24e (elem, heade, taile) !add element to list

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd24
