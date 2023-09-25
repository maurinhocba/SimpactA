 SUBROUTINE resvp6(heade,coora,euler,resid,energ,istop,bottom, &
                   top,coorb,coort,ifact,nstre,zigzag)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma) for shell elements
 !
 !********************************************************************
 USE npo_db, ONLY : psia
 IMPLICIT NONE

 LOGICAL, INTENT(IN) :: bottom,top  !flags for top and bottom surface computation
 INTEGER (kind=4), INTENT(IN) :: nstre  !
 INTEGER (kind=4), INTENT(OUT) :: istop,    & !flag for error detection
                               ifact(:)       !weigths for surface computation
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             euler(:,:)       !present total time

 REAL(kind=8),INTENT (IN OUT) :: resid(:,:),& !internal nodal forces
                                 energ(:),  & !bottom surface coordinates
                                 coorb(:,:),& !bottom surface coordinates
                                 coort(:,:)   !top surface coordinates
 TYPE (ele06), POINTER  :: heade              !pointer to first element
 LOGICAL, INTENT(IN) :: zigzag

 INTEGER (kind=4) ielem,nvarl,isec,l,n,matty,nlayr,osec,np,ndofe,secty
 LOGICAL :: plast

 REAL (kind=8), ALLOCATABLE :: bmsig(:,:)
 REAL (kind=8) thick,thico,x(3,nnode),t(3,nnode),lb(3,2,nnode),psi(2,nnode),sf(3)
 TYPE(ele06), POINTER :: e
 TYPE(section), POINTER :: sec
 TYPE(mater), POINTER :: mat


 !***  loop over all the elements

 osec = -1
 e => heade
 DO
   IF( .NOT.ASSOCIATED(e) )EXIT
   !       zero the internal forces array
   ielem = e%numel
   isec = e%matno
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     secty = sec%secty
     matty = sec%secdef(3)
     IF( sec%mabas > 0 )THEN
       mat => sec%mtbas
       plast = matty <= 2 .AND. mat%matdef(3) > 1
       nlayr = sec%iprop(1)
     ELSE     !for sectyp = 15
       plast = .FALSE.
       nlayr = 1          !integrated
     END IF
     thico = sec%rprop(1)
     IF(matty == 0)THEN  !layered isotropic
       np = 12
     ELSE                !in therms of forces and moments
       np = 21
       IF(matty == 5)THEN
         IF( sec%rprop(4) /= 0d0 )THEN
           sf = sec%rprop(4)
         ELSE
           sf = sec%rprop(71:73)
         END IF
       END IF
     END IF
     np = MIN(np,sec%secdef(2)) !avoids problems in debug
     osec = isec
     ndofe = 5
     IF( secty == 15 .AND. zigzag )ndofe = 7
     nvarl = nnode*ndofe
     IF( ALLOCATED(bmsig) )DEALLOCATE (bmsig)
     ALLOCATE( bmsig(ndofe,nnode))
   END IF

   !     gathers the actual coordinates and directors
   x  = coora(:,e%lnods)
   lb(:,1,:) = euler(1:3,e%lnods)
   lb(:,2,:) = euler(4:6,e%lnods)
   t  = euler(7:9,e%lnods)
   IF( zigzag )THEN
     DO n=1,nnode
       l = e%lnods(n)
       psi(1,n) = psia(1,l)*e%nangl(1,n)+psia(2,l)*e%nangl(2,n)
       psi(2,n) = psia(2,l)*e%nangl(1,n)-psia(1,l)*e%nangl(2,n)
     END DO
   END IF
   thick = thico
   CALL intrf6(nvarl,plast,ndofe,nstre,                  &
               matty,nlayr,sec%rprop(6:),sec%rprop(np:), &
               sf,x,t,lb,psi,e%nangl,e%jas,              &
               e%dvolu,e%cartd,e%jacin,                  &
               e%stra0,e%gamm0,e%strsg,                  &
               e%ambda,e%ehist,e%strap,e%stres,          &
               bmsig,energ,thick,istop,ielem)
   DO n=1,nnode
     l = e%lnods(n)
     resid(1:5,l) = resid(1:5,l) + bmsig(1:5,n)
     IF(zigzag)resid(7:8,l) = resid(7:8,l) + bmsig(6:7,n)
   END DO
   IF(top .OR. bottom)THEN
     DO n=1,nnode
       l = e%lnods(n)
       IF( top  )coort(:,l) = coort(:,l) + euler(7:9,l)*thick
       IF(bottom)coorb(:,l) = coorb(:,l) - euler(7:9,l)*thick
       ifact(l) = ifact(l) + 2
     END DO
   END IF
   e => e%next
 END DO
 RETURN
 END SUBROUTINE resvp6
