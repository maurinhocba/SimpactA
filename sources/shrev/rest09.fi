 SUBROUTINE rest09(nelem, nreqs,nnode,ngaus,nstre,axesc, heade,ngrqs, &
                   posgp,shape,deriv,weigh)
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) nelem, nreqs, nnode,ngaus,nstre,axesc
 INTEGER (kind=4), POINTER :: ngrqs(:)
 TYPE(ele09), POINTER :: heade
 REAL (kind=8) :: posgp(:),shape(:,:),deriv(:,:),weigh(:)

 ! local variables
 TYPE(ele09), POINTER :: e
 INTEGER (kind=4) :: ielem,i,isec,osec
 LOGICAL :: plast
 TYPE (section), POINTER :: sec      !auxiliar pointer

 READ (51) posgp, shape, deriv, weigh
 e => heade
 osec = -1
 DO ielem=1,nelem
   ALLOCATE (e)           !reserve memory for data
   ALLOCATE (e%lnods(nnode),       &
             e%phil(MAX(axesc,1)), &
             e%jac(ngaus),         &
             e%r0(ngaus),          &
             e%stra0(nstre,ngaus), &
             e%stres(nstre,ngaus), &
             e%sede0(ngaus),       &
             e%sedef(ngaus)        )
   READ(51)e%numel, e%matno, e%lnods, e%phil, e%jac, e%r0, e%stra0, &
           e%stres, e%sede0, e%sedef
   isec = e%matno          !associated section
   !     check if associated material exists
   IF( isec /= osec )THEN
     sec = psecs(isec)%p
     plast = sec%mtbas%matdef(3) > 1   !plasticiy
     osec = isec
   END IF

   IF( plast) THEN   !reserve space for Gauss variables
     ALLOCATE( e%ehist(5,ngaus), e%strap(nstre,ngaus) )
     READ (51)e%ehist, e%strap
   END IF
   e => e%next
 END DO

 IF ( nreqs > 0 )THEN
   ALLOCATE ( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)
 END IF

 RETURN
 END SUBROUTINE rest09
