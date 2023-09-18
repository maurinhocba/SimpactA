 SUBROUTINE elmd19(nelem, heade, taile, iwrit, eulrf)
 !******************************************************************
 !
 !*** READ element DATA for 2-D triangular element (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4),INTENT(IN):: iwrit
   INTEGER(kind=4),INTENT(INOUT):: nelem
   LOGICAL, INTENT(IN) :: eulrf
   TYPE(ele19),POINTER:: heade, taile !pointers to first and last elements
   !--- Local variables
   LOGICAL:: found,plast,mises
   INTEGER(kind=4):: i,ip,chnode,isec,nvare,mtype,osec
   TYPE(ele19),POINTER:: elem      !auxiliar pointer
   TYPE(section),POINTER:: sec     !auxiliar pointer


   ! *** READ the element nodal connections, and the property numbers.

   IF(iwrit == 1) WRITE (lures,"(//6X,'ELEMENT     MAT       NODE NUMBERS ', &
                               & 10x,'angle')",ERR=9999)

   osec = -1
   DO
     ! loop over elements

     CALL listen('ELMD19')          !read a card
     IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

     nelem = nelem + 1         !increase element counter
     CALL new_ele19e(elem)     !reserve memory for data

     elem%numel = INT (param(1))    !element label
     isec = INT (param(2))         !associated material
     IF( isec /= osec )THEN        !checks
       CALL sect_search (isec,found,sec,ip)
       IF (.NOT.found) CALL runend('ELMD19: Invalid Section Number used')
       mtype = sec%secty        !check section type for element
       IF( mtype /= 0 ) CALL runend('ELMD19: Invalid sect. type for SOLID')
       mtype = sec%mtbas%mtype       !check material type for element
       IF( mtype /= 4 )CALL runend('ELTTR2D: Only Material 4 can be used')
       plast = sec%mtbas%matdef(3) > 1  !plasticity
       mises = sec%mtbas%matdef(3) == 2 !Mises plasticity
       ! visco = sec%mtbas%matdef(6) == 3 !Viscoplasticity ... not yet
       nvare = 1
       IF( plast )THEN
         nvare = 5   !number of internal variables per integration point (plastic)
       END IF
       IF( mtype == 4 .AND. mises ) nvare = nvare + 1 !von Mises ... check this WBC
       ! IF( mtype == 4 .AND. visco ) nvare = nvare + 1 !add efective plastic strain rate ... not yet
       osec = isec
     END IF
     elem%matno = ip        !store section position in array

     elem%lnods(1:nnode) = INT (param(3:5))  !connectivities
     elem%angle = param(6)   ! angle between Default Ort_1 and local Ort_1
     !    echo according to flag
     IF(iwrit == 1) WRITE(lures,"(i13,i8,3i7,F10.4)",ERR=9999) elem%numel,isec, &
                                                     elem%lnods(1:3),elem%angle
     !     original nodes are changed
     DO i=1,nnode
       elem%lnods(i) = chnode(elem%lnods(i))
     END DO

     IF( plast )THEN
       ALLOCATE (elem%gausv(nvare),elem%strpl(5)) !,ngaus))   !reserve space for Gauss variables
       IF( mises .AND. mtype == 4 )THEN
         elem%strpl = (/ 0d0, 0d0, 0d0, 0d0, 0d0 /) 
         !IF (eulrf .AND. visco) THEN !initializes elastic Finger for spatial viscoplastic formulation
         !  elem%gausv = (/ 1d0, 1d0, 0d0, 1d0, 0d0, 0d0, 0d0 /)  ... not yet

         IF (eulrf) THEN !initializes elastic Finger tensor for Spatial Formulation
           elem%gausv = (/ 1d0, 1d0, 0d0, 1d0, 0d0, 0d0 /)
         ELSE  !initializes plastic gradient tensor for Lagrangian Formulation
           elem%gausv = (/ 1d0, 0d0, 0d0, 1d0, 1d0, 0d0 /)
         END IF

       END IF
     END IF

     CALL add_ele19e (elem, heade, taile) !add element to list

   END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd19
