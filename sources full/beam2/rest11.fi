 SUBROUTINE rest11(nelem,nreqs,ntype,ngaus,nbn,heade,taile,ngrqs, &
                   nhead,stint)
 !******************************************************************
 !
 !*** restore an element set at restart
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) nelem, nreqs, ntype, ngaus, nbn
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL (kind=8), POINTER :: stint(:,:)
 TYPE(ele11), POINTER :: heade,taile
 TYPE(nodei), POINTER :: nhead


 ! local variables
 TYPE(ele11), POINTER :: e
 TYPE(nodei), POINTER :: ni,ntail
 INTEGER (kind=4) :: ielem,i,nvare,nlayr
 LOGICAL :: fcur

 DO ielem=1,nelem
   CALL new_ele11e(e)     !reserve memory for data
   READ(51)e%numel, e%matno, e%lnods, e%l0, e%cab,  &
           e%lambd, stint(:,ielem), fcur
   IF( ntype == 3 )THEN
     ALLOCATE ( e%ro(2,ngaus) )
     READ(51) e%ro
   END IF

   READ (51)nvare,nlayr
   IF( nvare > 0 )THEN
     ALLOCATE( e%gausv(nvare,nlayr,ngaus) )
     READ (51)e%gausv
   END IF

   IF( fcur) THEN   !reserve space for curvatures
     ALLOCATE( e%cur(2) )
     READ (51)e%cur
   END IF
   CALL add_ele11e (e, heade, taile) !add element to list
 END DO

 IF ( nreqs > 0 )THEN
   ALLOCATE ( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)
 END IF

 IF ( nbn > 0 )THEN
   DO i=1,nbn
     ALLOCATE ( ni )
     READ (51) ni%nn,ni%ba
     ALLOCATE( ni%lnods(0:ni%nn), ni%alph0(1:ni%nn),ni%bb(2,0:ni%nn))
     READ (51) ni%lnods,ni%alph0
     ni%bb = 0d0
     IF( i == 1 )THEN
        nhead => ni
        ntail => ni
     END IF
     ntail%next => ni
     ntail => ntail%next
   END DO
 END IF

 RETURN
 END SUBROUTINE rest11
