 SUBROUTINE elmd25 (nelem, heade, taile, iwrit)
 !******************************************************************
 !
 !*** READ element DATA for shell element BSQ  (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER ::  nnode=4
 INTEGER (kind=4) :: nelem, iwrit
 TYPE (ele25), POINTER :: heade,taile !pointers to first and last elements

 INTEGER (kind=4) :: i,j,ip,chnode,isec,secty,nlayr,osec,mtype,nvare
 INTEGER (kind=4) :: secs(2) = (/ 12,13 /)      !only acceptable section type
 TYPE (ele25), POINTER :: elem      !auxiliar pointer
 TYPE(section), POINTER :: sec
 LOGICAL :: found


 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//6X,'ELEMENT     MAT       NODE NUMBERS ', &
                             & 10x,'angle')",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('BSQ   ')         !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele25e(elem)     !reserve memory for data

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))          !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('EL_BSQ : Invalid Section Number used')
     secty = sec%secty       !check section type for element
     IF( ALL(secty /= secs) ) CALL runend('EL_BSQ : Invalid sect. type for  BSQ ')
     nlayr = sec%iprop(1)     !number of layers
     nvare = 0   !number of internal variables per integration point (plastic)
     IF( sec%mabas >= 0 )THEN
       mtype = sec%mtbas%mtype  !material type
       IF( mtype == 4) CALL runend('EL_BSQ: Invalid mater type for BSQ ')
       IF( mtype == 30 ) nvare = 3 + sec%mtbas%matdef(7) !user defined
     END IF
     osec = isec
   END IF
   elem%matno = ip        !store section position in array

   elem%lnods(1:nnode) = INT (param(3:6))  !connectivities
   elem%lnods(5:12) = 0                     !initializes LSIDE
   elem%angle = param(7)   ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i13,i8,4i7,F30.4)",ERR=9999)elem%numel,isec, &
                                                   elem%lnods(1:nnode),elem%angle

   !     original nodes are changed
   DO i=1,nnode
     j = elem%lnods(i)
     IF( j > 0 ) elem%lnods(i) = chnode(j)
   END DO

   IF( nvare > 0 )THEN   !reserve space for Gauss variables
     ALLOCATE (elem%gausv(nvare,nlayr))
     elem%gausv = 0d0                     !initializes
   END IF

   CALL add_ele25e (elem, heade, taile) !add element to list

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd25
