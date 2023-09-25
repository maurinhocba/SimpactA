 SUBROUTINE elmda9(nelem,nnode,nstre,ngaus,axesc,heade,taile,iwrit)
 !******************************************************************
 !
 !*** READ element DATA for 2-D beam/shell elements
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nnode,nstre,ngaus,axesc,iwrit
 INTEGER (kind=4), INTENT(IN OUT) :: nelem
 TYPE (ele09), POINTER :: heade,taile !pointers to first and last elements
 TYPE (section), POINTER :: sec      !auxiliar pointer

 ! local variables
 INTEGER (kind=4) i,isec,mtype,osec,ip
 TYPE(ele09), POINTER :: elem
 LOGICAL found,plast
 INTEGER(kind=4) chnode
 INTEGER (kind=4) :: secs = 43

! *** READ the element nodal connections, and the property numbers.

 osec = -1
 IF(iwrit == 1) WRITE (lures,"(/' ELEMENT  SECC      NODE NUMBERS ')",ERR=9999)
 DO
   ! loop over elements

   CALL listen('ELMDA6')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   ALLOCATE (elem)           !reserve memory for data
   ALLOCATE (elem%lnods(nnode),       &
             elem%phil(MAX(ABS(axesc),1)),  &
             elem%jac(ngaus),         &
             elem%r0(ngaus),          &
             elem%stra0(nstre,ngaus), &
             elem%stres(nstre,ngaus), &
             elem%sede0(ngaus),       &
             elem%sedef(ngaus)        )
   ! initializes
   elem%phil = 0d0

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('SHREV:  Invalid Section Number used')
     mtype = sec%secty       !check section type for element
     IF( mtype /= secs ) CALL runend('SHREV : Invalid sect. type for SHREV')

     !nlayr = sec%iprop(1)     !number of layers
     !nucom = sec%iprop(2)     !number of material components
     !stype = sec%secdef(3)    !section constitutive type
     plast = sec%mtbas%matdef(3) > 1   !plasticiy
     IF( plast .AND. sec%mtbas%matdef(4) > 4 )THEN  !isotropic hardening type
       WRITE(lures,"('For SHREV element type',/,                    &
       & '   use only the following isotropic hardening types:',/,  &
       & '                      IHNONE,IHLINE,IHLUDW or IHSATU')",ERR=9999)
       WRITE(55   ,"('For SHREV element type',/,                    &
       & '   use only the following isotropic hardening types:',/,  &
       & '                      IHNONE,IHLINE,IHLUDW or IHSATU')",ERR=9999)
       CALL runend('SHREV:  Invalid Isotropic hard type')
     END IF
     osec = isec
   END IF
   elem%matno = ip      !store section position in array
   elem%lnods(1:nnode) = INT (param(3:2+nnode))  !connectivities
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i8,i5,3i7)",ERR=9999) elem%numel,isec,elem%lnods
   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   IF( plast )THEN   !reserve space for Gauss variables
     ALLOCATE( elem%ehist(5,ngaus), elem%strap(nstre,ngaus) )
     elem%ehist = 0d0
     elem%strap = 0d0
   ELSE
     NULLIFY( elem%ehist, elem%strap)
   END IF
   CALL add_ele09e (elem, heade, taile) !add element to list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmda9
