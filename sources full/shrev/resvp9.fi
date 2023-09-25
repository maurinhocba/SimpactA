 SUBROUTINE resvp9(ntype,nelem,nnode,ngaus,axesc,nstre,heade, &
                   coora,euler,resid,weigh,shape,deriv,       &
                   energ,istop,bottom,top,coorb,coort,ifact)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma) for beam elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ntype,nelem,nnode,ngaus,axesc,nstre
 INTEGER (kind=4), INTENT(IN OUT) :: istop,ifact(:)
 REAL (kind=8), INTENT(IN) :: coora(:,:),euler(:,:),  &
                             weigh(:),shape(:,:),deriv(:,:)
 REAL (kind=8), INTENT(OUT) :: resid(:,:),coort(:,:),coorb(:,:),energ(:)
 LOGICAL, INTENT(IN) :: bottom,top
 TYPE(ele09), POINTER :: heade


 INTEGER (kind=4) ielem,l,n,isec,np,nd,osec,matty
 REAL (kind=8) phig(nnode),x(2,nnode),bts(3,nnode),tn(2,nnode),props(35)
 TYPE(ele09), POINTER :: e
 TYPE(section), POINTER :: sec
 LOGICAL :: plast

 osec = -1
 IF( ntype == 1 )THEN
   nd = 11
   np = 14
 ELSE
   nd = 21
   np = 26
 END IF
 !***  loop over all the elements
 e => heade
 DO ielem=1,nelem
   !   gathers the actual coordinates and local coordinate transf.
   x    = coora(:,e%lnods)
   phig = euler(1,e%lnods)
   isec  = e%matno
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     matty = sec%secdef(3)
     plast = sec%mtbas%matdef(3) > 1
     props = sec%rprop
     osec  = isec
   END IF
   CALL intrf9(nstre,ntype,plast,nnode,ngaus,matty,axesc,     &
              weigh,shape,deriv,props(nd:),props(5:),props(np:), &
              e%jac,e%r0,e%stra0,e%stres,e%sede0,e%sedef,     &
              e%ehist,e%strap,bts,x,phig,e%phil,              &
              energ,istop,tn,props(3))
   DO n=1,nnode
     l = e%lnods(n)
     resid(1:ndofe,l) = resid(1:ndofe,l) + bts(1:ndofe,n)
   END DO
   IF(top .OR. bottom)THEN
     DO n=1,nnode
       l = e%lnods(n)
       IF( top  )coort(1:2,l) = coort(1:2,l) + tn(1:2,n)
       IF(bottom)coorb(1:2,l) = coorb(1:2,l) - tn(1:2,n)
       ifact(l) = ifact(l) + 2
     END DO
   END IF
   e => e%next
 END DO
 RETURN
 END SUBROUTINE resvp9
