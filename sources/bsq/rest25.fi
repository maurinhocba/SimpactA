 SUBROUTINE rest25 (nelem, nreqs, heade, taile, ngrqs, nbs, bhead, btail, stint)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 nelem, & !number of elements in the set
                                 nbs      !number of branching sides in the set
 INTEGER (kind=4), POINTER  :: ngrqs(:)   !elements for requested output
 REAL(kind=8), POINTER :: stint(:,:)         !forces, moments and shear
 TYPE (ele25), POINTER  :: heade,taile    !pointers to first and last element
 TYPE (sideb), POINTER  :: bhead,btail    !pointers to first and last side

 INTEGER (kind=4) :: i,j,nlayr,nvare,ielem,nn
 TYPE (ele25), POINTER :: e
 TYPE (sideb), POINTER :: sb

 IF( nbs > 0 )THEN
   DO i=1,nbs
     ! read  variables associated to the side
     READ (51) nn
     ALLOCATE(sb)
     ALLOCATE( sb%lnods(-1:2*nn), & ! j,k, i=1,2*nn (-) inverse order in conns
               sb%alph0(nn),    & ! i=1,nn-1 angles between elems. nn = 1/ls
               sb%fc(nn,nn),    & ! (1:nn,1:nn-1) matrix F  (1:nn,nn) relative stiffness
               sb%t0(3,nn),  &
               sb%tb(3,nn),  &
               sb%c(4,nn,nn),      & ! normal derivative of shape functions
               sb%gamma(nn),    & ! distortion at each element
               sb%bbr(3,-1:2*nn) )   ! average rotation B - matrix
     READ (51) sb%lnods,sb%alph0,sb%fc,sb%c
     !READ (51) sb%gamma,sb%bb  !unnecessary, rebuild in RESV13
     sb%nn = nn
     IF(i == 1)THEN   !for first branching side
       bhead => sb      !keep first pointer
       btail => sb      !auxiliar last pointer
     ELSE
       btail%next => sb !updates pointer to new branching side
       btail => sb      !updates pointer
     END IF
     NULLIFY (btail%next)  !unnecessary
     NULLIFY (sb)          !unnecessary
   END DO
 END IF

 DO ielem=1,nelem    !for each element

   CALL new_ele25e(e)    !get memory for an element
   READ (51) e%numel,e%matno,e%lnods,e%lside,nvare,nlayr !,e%ldv2sd
   READ (51) e%area,e%lb,e%angle,e%cartd,e%normd,e%ns,e%a0,e%ci,e%gamma, &
             e%stra0,e%stra1,stint(:,ielem)
   IF( nvare > 0 )THEN
     ALLOCATE (e%gausv(nvare,nlayr))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,nvare),j=1,nlayr) !read internal variables
   ELSE
     NULLIFY( e%gausv )
   END IF
   IF( ANY(e%lnods < 0 ))THEN !if a branching side exist
     DO i=1,4
       j = fn(i)              !first node of associated side element (local)
       IF( e%lnods(j) >= 0 )CYCLE
       j = -e%lnods(j)        !position in side connectivities
       sb => bhead
       DO
         IF( sb%nn*2 >= j )THEN  !check number of nodes
           IF( ABS(sb%lnods(j)) == e%lnods(i)) THEN !node found
            !  transfer information
             e%si(i)%p => sb      !keep a pointer to the branching side
             EXIT
           END IF
         END IF
         sb => sb%next
       END DO
     END DO
   END IF
   CALL add_ele25e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 ) THEN
   ALLOCATE( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF
 RETURN
 END SUBROUTINE rest25
