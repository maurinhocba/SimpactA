 SUBROUTINE resvp2(ndime,nelem,coora,resid,heade)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma) for truss elements
 !
 !********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ndime,nelem
 REAL (kind=8) coora(:,:),resid(:,:)
 TYPE(ele02), POINTER :: heade

 INTEGER (kind=4) ielem,l,n,imat,gtype,omat
 REAL (kind=8) x(ndime,nnode),bts(ndime,nnode),es,xl,dx(ndime),facto
 TYPE(ele02), POINTER :: e
 TYPE(section), POINTER :: sec
 LOGICAL :: newm


 !***  loop over all the elements
 omat = -1               !initializes
 e => heade              !point to first element
 DO ielem = 1,nelem      !loop over the set of elements
   imat  = e%matno          !section of present element
   newm = imat /= omat      !compare with previous section
   IF ( newm )THEN            !if a new section
     sec => psecs(imat)%p                !search and point
     omat = imat                         !remember section number
     gtype = sec%iprop(1)                !type of stress-strain relation
   END IF
   !       gathers the actual coordinates
   x = coora(:,e%lnods)                  !present coordinates
   dx = x(:,nnode) - x(:,1)   !distance vector
   CALL stra02(ndime,gtype,e%gausv(1),es,dx,xl) !computes strain measure

   CALL stre02 (sec, es, e%gausv(:), .TRUE., xl ) !computes stress

   facto = e%gausv(4) / e%gausv(1)/xl

   bts(1:ndime,1) = - facto*dx
   bts(1:ndime,2) =   facto*dx

   DO n = 1,nnode
     l = e%lnods(n)
     resid(1:ndime,l) = resid(1:ndime,l) + bts(1:ndime,n)
   END DO
   e => e%next
 END DO
 RETURN
 END SUBROUTINE resvp2
