 SUBROUTINE elmda6(nelem,heade,taile,iwrit,nstre)
 !******************************************************************
 !
 !*** READ element DATA for shell elements
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4) :: nelem,iwrit,nstre
 TYPE (ele06), POINTER :: heade,taile !pointers to first and last elements

 !Local variables
 INTEGER (kind=4) isec,i,nlayr,mtype,osec,stype,nucom,ip,secty
 TYPE (ele06), POINTER :: elem      !auxiliar pointer
 TYPE (section), POINTER :: sec      !auxiliar pointer

 LOGICAL found, plast
 INTEGER(kind=4) :: chnode
 INTEGER (kind=4) :: secs(2) = (/ 14,15 /)

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//'      Element     Secc     ', &
                        &        '    Node Numbers             angle  ')",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('ELMDA6')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele06e (elem)           !reserve memory for data
   ALLOCATE(elem%strsg(nstre,ngaus))
   elem%strsg = 0d0                 !initializes integrated stresses
   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('SHELQ:  Invalid Section Number used')
     mtype = sec%secty       !check section type for element
     IF( .NOT. ANY(mtype == secs) ) CALL runend('SHELQ : Invalid sect. type for SHELQ')

     nlayr = sec%iprop(1)     !number of layers
     nucom = sec%iprop(2)     !number of material components
     secty = sec%secty        !section constitutive type
     stype = sec%secdef(3)    !section constitutive type
     plast = .FALSE.
     IF( stype <= 2 .AND. secty == 14) plast = sec%mtbas%matdef(3) > 1   !plasticiy
     osec = isec
   END IF
   elem%matno = ip      !store section position in array

   elem%lnods(1:nnode) = INT (param(3:6))  !connectivities
   elem%angle = param(7)   ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i13,i8,4i7,F14.4)",ERR=9999) elem%numel,isec, &
                                                   elem%lnods,elem%angle
   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   IF( plast )THEN   !reserve space for Gauss variables
     IF( stype == 0 )THEN  !layered large strain plasticity
       ALLOCATE( elem%ehist(5,ngaus), elem%strap(6*nlayr,ngaus) )
       ALLOCATE( elem%stres(5,nlayr,ngaus) )
       elem%ehist = 0d0
       elem%strap = 0d0
       elem%stres = 0d0
     ELSE   !solid isotropic
       ALLOCATE( elem%ehist(5,ngaus), elem%strap(nstre,ngaus) )
       elem%ehist = 0d0
       elem%strap = 0d0
     END IF
   END IF
   IF( nstre > 8 .AND. secty == 15) ALLOCATE(elem%nangl(2,nnode),elem%jas(2,nnode))
   CALL add_ele06e (elem, heade, taile) !add element to list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmda6
