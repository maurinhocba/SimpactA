 SUBROUTINE elmda7(nelem,nnode,heade,taile,rtype,iwrit,ngaus,ngamm,nstre)
 !******************************************************************
 !
 !*** READ element DATA for shell elements
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4) :: nelem,nnode,rtype,iwrit,ngaus,ngamm,nstre
 TYPE (ele07), POINTER :: heade,taile !pointers to first and last elements

 !Local variables
 INTEGER (kind=4) isec,i,nlayr,mtype,osec,stype,nucom,ip,secty
 TYPE (ele07), POINTER :: elem      !auxiliar pointer
 TYPE (section), POINTER :: sec      !auxiliar pointer

 LOGICAL found, plast
 INTEGER(kind=4) :: chnode
 INTEGER (kind=4) :: secs(2) = (/ 14,15 /)


 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//'   Element  Secc            Node ', &
                        &        'Numbers                  angle  ')",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('ELMDA7')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele07e (elem)    !reserve memory for data
   ALLOCATE( elem%dvolu(ngaus),elem%lnods(nnode),elem%jacin(2,2,ngaus), &
             elem%stra0(6,ngaus),elem%gamm0(ngamm),elem%strsg(nstre,ngaus),     &
             elem%qbar(ngamm),elem%ambda(2,ngaus))
   IF(ngaus == 3)THEN
     ALLOCATE(elem%cartd(nnode,2,ngaus))
   ELSE
     ALLOCATE(elem%cartd(3,2,ngaus))
     IF(nnode == 9)THEN
       ALLOCATE(elem%cd(4,2,3))
       elem%lnods(7:9) = 0
     END IF
   END IF
   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('SHELQ:  Invalid Section Number used')
     mtype = sec%secty       !check section type for element
     IF( .NOT. ANY(mtype == secs) ) CALL runend('SHELT : Invalid sect. type for SHELT')

     nlayr = sec%iprop(1)     !number of layers
     nucom = sec%iprop(2)     !number of material components
     secty = sec%secty        !section constitutive type
     stype = sec%secdef(3)    !section constitutive type
     plast = .FALSE.
     IF( stype <= 2 .AND. secty == 14) plast = sec%mtbas%matdef(3) > 1   !plasticiy
     osec = isec
   END IF
   elem%matno = ip      !store section position in array

   elem%lnods(1:6) = INT (param(3:8))  !connectivities
   elem%angle = param(9)   ! angle between Default Ort_1 and local Ort_1
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i10,i5,6i6,F14.4)",ERR=9999) elem%numel,isec, &
                                                   elem%lnods,elem%angle
   !     original nodes are changed
   DO i=1,6
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   IF( ngaus == 1 )THEN   ! switch nodes of TLLL version
     i = elem%lnods(4)
     elem%lnods(4) = elem%lnods(5)
     elem%lnods(5) = elem%lnods(6)
     elem%lnods(6) = i
   END IF

   IF( rtype > 1 )ALLOCATE( elem%beta(nbeta) )
   IF( stype == 0 .AND. rtype /= 0) & !check
     CALL runend('ELMDA7: Invalid mat.-Shell type com')

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
       NULLIFY (elem%stres)   !
     END IF
   ELSE
     NULLIFY (elem%ehist,elem%strap,elem%stres)   !
   END IF
   IF( nstre > 8 .AND. sec%secty == 15)THEN
     IF( ngaus == 1 )THEN
       ALLOCATE(elem%nangl(2,3),elem%jas(2,ngamm))
     ELSE
       ALLOCATE(elem%nangl(2,6),elem%jas(2,ngamm))
     END IF
   END IF
   elem%strsg = 0d0
   elem%qbar = 0d0
   CALL add_ele07e (elem, heade, taile) !add element to list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmda7
