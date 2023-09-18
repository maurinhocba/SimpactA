 SUBROUTINE elmd15 ( nelem, heade, taile, iwrit )
 !******************************************************************
 !
 !*** READ element DATA for shell element RBST
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER ::  nnode=3
 INTEGER (kind=4) :: nelem, iwrit
 TYPE (ele15), POINTER :: heade,taile !pointers to first and last elements

 INTEGER (kind=4) :: i,chnode,isec,secty,nlayr,osec,mtype,nvare,ip
 INTEGER (kind=4) :: secs(2) = (/ 12,15 /)      !only acceptable section type
 TYPE (ele15), POINTER :: elem      !auxiliar pointer
 TYPE(section), POINTER :: sec
 LOGICAL :: found


 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//6X,'ELEMENT     MAT       NODE NUMBERS ', &
                             & 10x,'angle')",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('ELMD15')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele15e(elem)     !reserve memory for data

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))          !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELRBST: Invalid Section Number used')
     secty = sec%secty       !check section type for element
     IF( ALL(secty /= secs) ) CALL runend('ELRBST: Invalid sect. type for RBST ')
     nlayr = sec%iprop(1)     !number of layers
     mtype = sec%mtbas%mtype  !material type
     IF( mtype == 4) CALL runend('ELRBST: Invalid mater type for RBST ')
     nvare = 0   !number of internal variables per integration point (plastic)
     IF( mtype == 30 ) nvare = 3 + sec%mtbas%matdef(7) !user defined
     osec = isec
   END IF
   elem%matno = ip        !store section Label

   elem%lnods(1:nnode) = INT (param(3:5))  !connectivities
   elem%lnods(4:6) = 0                     !initializes LSIDE
   elem%angle = param(6)   ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i15,i8,3i7,F14.4)",ERR=9999) elem%numel,isec, &
                                                    elem%lnods(1:3),elem%angle

   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   IF( nvare > 0 )THEN   !reserve space for Gauss variables
     ALLOCATE (elem%gausv(nvare,nlayr))
     elem%gausv = 0d0                     !initializes
   ELSE
     NULLIFY( elem%gausv )
   END IF

   elem%mems = .FALSE.                  !initializes

   CALL add_ele15e (elem, heade, taile) !add element to list

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd15
