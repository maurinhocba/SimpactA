 SUBROUTINE resvp8(ndime,nelem,nnode,ngaus,axesc,coora,  &
              euler,velnp,resid,weigh,shape,deriv,heade,istop)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma) for beam elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: ndime,nelem,nnode,ngaus,axesc,   &
                    istop
 REAL (kind=8) :: coora(:,:),euler(9,*),resid(:,:),weigh(:),    &
                  shape(:,:),deriv(:,:),velnp(:,:)
 TYPE (ele08), POINTER :: heade

 INTEGER (kind=4), PARAMETER :: ndofe = 6, nstre = 6
 REAL (kind=8) x(ndime,nnode),v(ndofe,nnode),lbd(3,3,nnode),          &
               bmsig(ndofe,nnode)

 INTEGER (kind=4) ielem,l,n,isec,osec
 TYPE(section), POINTER :: sec
 TYPE (ele08), POINTER :: e
 LOGICAL :: isotro,visco,plast

 !***   loop over all the elements

 osec = -1
 e => heade
 DO ielem=1,nelem
   isec = e%matno
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     isotro = sec%iprop(1) == 0
     visco = sec%mtbas%mtype == 8
     plast = .FALSE.
     IF( isotro ) plast = sec%rprop(21) > 0d0
     osec  = isec
   END IF
   ! gathers the actual coordinates and local coordinate transf.
   DO n=1,nnode      !for each element node
     l = e%lnods(n)  !global node
     x(1:ndime,n) = coora(1:ndime,l)  !actual coordinates
     lbd(:,:,n)   = RESHAPE(euler(1:9,l),(/3,3/))      !actual local system
     v(:,n) = velnp(:,l)            !recover velocity
   END DO

   CALL intrf8(ndime,ndofe,ngaus,nstre,nnode,axesc,     &
               sec%rprop(:),weigh,shape,deriv,          &
               e%jac,e%stra0,e%stran,   &
               e%stres,e%epdef,         &
               e%sedef,bmsig,x,v,lbd,           &
               e%llbd,istop,e%numel,isotro,visco,plast)
   IF(istop == 1) RETURN
   DO n=1,nnode
     l = e%lnods(n)
     resid(1:ndofe,l) = resid(1:ndofe,l)+bmsig(1:ndofe,n)
   END DO
   e => e%next
 END DO
 RETURN
 END SUBROUTINE resvp8
