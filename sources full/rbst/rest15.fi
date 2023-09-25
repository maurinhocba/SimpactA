 SUBROUTINE rest15 (nelem, nreqs, heade, taile, ngrqs, nrf, rhead, rtail, stint)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 nelem, & !number of elements in the set
                                 nrf      !number of reinforced sides in the set
 INTEGER (kind=4), POINTER  :: ngrqs(:)   !elements for requested output
 REAL(kind=8), POINTER :: stint(:,:)         !forces, moments and shear
 TYPE (ele15), POINTER  :: heade,taile    !pointers to first and last element
 TYPE (reinf), POINTER  :: rhead,rtail    !pointers to first and last side

 INTEGER (kind=4) :: i,j,nlayr,nvare,ielem,k,l
 TYPE (ele15), POINTER :: e
 TYPE (reinf), POINTER :: sb

 DO i=1,nrf
   ! read  variables associated to the side
   ALLOCATE(sb)
   READ (51) sb%numel,sb%lnods,sb%secno
   READ (51) sb%l0,sb%curv,sb%fc,sb%dc,sb%stint
   IF(i == 1)THEN   !for first reinforced side
     rhead => sb      !keep first pointer
   ELSE
     rtail%next => sb !updates pointer to new branching side
   END IF
   rtail => sb      !auxiliar last pointer
   !NULLIFY (btail%next)  !unnecessary
   !NULLIFY (sb)          !unnecessary
 END DO

 DO ielem=1,nelem    !for each element

   ALLOCATE (e)      !get memory for an element
   READ (51) e%numel,e%matno,e%lnods,e%lside,nvare,nlayr
   READ (51) e%area1,e%lb,e%angle,e%a,e%b,e%c,e%cd,e%a0,e%ci,e%gamma, &
             e%stra1,stint(:,ielem),e%mems
   IF( nvare > 0 )THEN
     ALLOCATE (e%gausv(nvare,nlayr))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,nvare),j=1,nlayr) !read internal variables
   ELSE
     NULLIFY( e%gausv )
   END IF
   CALL add_ele15e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 ) THEN
   ALLOCATE( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 sb => rhead             !rather long loop
 DO i=1,nrf
   DO j=1,2
     k = sb%lnods(j+2)
     IF( k == 0 )CYCLE
     l = MOD(k,10)
     k = k/10
     e => heade
     DO ielem=1,k-1
       e => e%next
     END DO
     e%si(l)%p => sb
   END DO
   sb => sb%next      !auxiliar last pointer
 END DO

 RETURN
 END SUBROUTINE rest15
