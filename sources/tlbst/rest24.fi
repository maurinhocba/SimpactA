 SUBROUTINE rest24 (nelem, nreqs, heade, taile, shearc, ngrqs, nbs, &
                    bhead, btail, stint, factors, ninv, moments)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element TLBST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 nelem, & !number of elements in the set
                                 nbs,   & !number of branching sides in the set
                                 shearc   !flag for shear forces computation
 INTEGER (kind=4), POINTER  :: ngrqs(:)   !elements for requested output
 REAL(kind=8), POINTER :: stint(:,:)         !forces, moments and shear
 TYPE (ele24), POINTER  :: heade,taile    !pointers to first and last element
 TYPE (sideb), POINTER  :: bhead,btail    !pointers to first and last side
 REAL(kind=8), POINTER :: factors(:)      !factors to smooth moments
 REAL(kind=8), POINTER :: moments(:,:)    !smoothed nodal moments
 INTEGER(kind=4), POINTER :: ninv(:)      !inverse nodal relation


 INTEGER (kind=4) :: i,j,nlayr,nvare,ielem,nn
 TYPE (ele24), POINTER :: e
 TYPE (sideb), POINTER :: sb

 IF( nbs > 0 )THEN
   DO i=1,nbs
     ! read  variables associated to the side
     READ (51) nn
     ALLOCATE(sb)
     ALLOCATE( sb%lnods(-1:nn), & ! j,k, i=1,nn (-) inverse order in conns
               sb%alph0(nn),    & ! i=1,nn-1 angles between elems. nn = 1/ls
               sb%fc(nn,nn),    & ! (1:nn,1:nn-1) matrix F  (1:nn,nn) relative stiffness
               sb%c(3,nn),      & ! normal derivative of shape functions
               sb%gamma(nn),    & ! distortion at each element
               sb%bb(3,-1:nn) )   ! average rotation B - matrix
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

   CALL new_ele24e(e)      !get memory for an element
   READ (51) e%numel,e%matno,e%lnods,e%lside,nvare,nlayr
   READ (51) e%area1,e%lb,e%angle,e%a,e%b,e%c,e%cd,e%a0,e%ci,e%gamma, &
             e%stra1,stint(:,ielem),e%mems
   IF( nvare > 0 )THEN
     ALLOCATE (e%gausv(nvare,nlayr))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,nvare),j=1,nlayr) !read internal variables
   ELSE
     NULLIFY( e%gausv )
   END IF
   IF( ANY(e%lnods < 0 ))THEN !if a branching side exist
     DO i=1,3
       IF( e%lnods(i+3) >= 0 )CYCLE
       j = -e%lnods(i+3)   !position in side connectivities
       sb => bhead
       DO
         IF( sb%nn >= j )THEN  !check number of nodes
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
   CALL add_ele24e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 ) THEN
   ALLOCATE( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 IF (shearc < 0 )THEN
   READ (51) nn
   ALLOCATE(factors(nn),ninv(npoin),moments(3,nn))
   READ (51) (factors(i),i=1,nn)
   READ (51) (ninv(i),i=1,npoin)
 END IF

 RETURN
 END SUBROUTINE rest24
