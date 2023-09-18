 SUBROUTINE elmd17 (nelem, heade, taile, iwrit, ngaus)
 !******************************************************************
 !
 !*** READ element DATA for 2-D QUADL element  (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem, iwrit, ngaus
 TYPE (ele17), POINTER :: heade,taile !pointers to first and last elements

 LOGICAL found,plast
 INTEGER (kind=4) :: i,ip,chnode,isec,nvare,mtype,osec
 INTEGER (kind=4) :: secs(3) = (/ 1,5,6 /)  !acceptable material types
 TYPE (ele17), POINTER :: elem
 TYPE (section), POINTER :: sec

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE(lures,"(//6X,'ELEMENT     MAT       NODE NUMBERS ',&
                           &    10x,'angle')",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('ELMD17')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   ALLOCATE (elem)           !reserve memory for data

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   IF( isec /= osec )THEN        !checks
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELMD17: Invalid Section Number used')
     mtype = sec%secty        !check section type for element
     IF( mtype /= 0 ) CALL runend('ELMD17: Invalid sect. type for SOLID')
     mtype = sec%mtbas%mtype       !check material type for element
     IF( .NOT. ANY(mtype == secs) )  &
       CALL runend('ELMD17: Invalid Mat. type for QUADL')
     plast = sec%mtbas%matdef(3) > 1  !plasticity
     IF( plast )THEN
       nvare = 5   !number of internal variables per integration point (plastic)
       IF(sec%mtbas%matdef(5) > 1) nvare = nvare + 4  !kinematic hardening
     ELSE
       nvare = 1
     END IF
     IF( mtype == 5 .AND. plast ) &
       CALL runend('QUADL: Orthotropic and plastic Not possible')
     osec = isec
   END IF

   elem%matno = ip       !store position in section array

   elem%lnods(1:nnode) = INT (param(3:6))  !connectivities
   elem%angle = param(7)   ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i13,i8,4i7,F17.4)",ERR=9999) elem%numel,isec, &
                               elem%lnods(1:4),elem%angle
   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   ALLOCATE (elem%dvol(ngaus),elem%stint(4,ngaus))  !reserve space for variables
   ALLOCATE (elem%cartd(nnode,2,MAX(ngaus,2)))  !reserve space for variables
   IF( plast )THEN
     ALLOCATE (elem%gausv(nvare,ngaus))   !reserve space for Gauss variables
   ELSE
     ALLOCATE (elem%gausv(1,1))   !reserve space for Gauss variables
   END IF
   elem%gausv = 0d0                     !initializes

   CALL add_ele17e (elem, heade, taile) !add element to list

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd17
