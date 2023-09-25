 SUBROUTINE resvp7(heade,coora,euler,resid,energ,istop,ap1,ngamm,shape,posgp,bottom, &
                   top,coorb,coort,ifact,stype,ngaus,nnode,stabq,nstre,zigzag)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma) for shell elements
 !
 !********************************************************************
 USE npo_db, ONLY : psia
 USE ctrl_db, ONLY : npoin
 IMPLICIT NONE

 LOGICAL, INTENT(IN) :: bottom,top  !flags for top and bottom surface computation
 INTEGER (kind=4), INTENT(IN) ::  stype, &    !formulation type
                                  ngamm, &    !number of shear values
                                  nnode, &    !number of nodes in LNODS
                                  nstre, &    !number of stress components
                                  ngaus       !number of integration points
 INTEGER (kind=4), INTENT(OUT) :: istop,    & !flag for error detection
                               ifact(:)       !weigths for surface computation
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             euler(:,:),    & !present total time
                             ap1(:,:,:),    & !
                             shape(:,:),    & !
                             posgp(:,:),    & !
                             stabq            !stabilization factor

 REAL(kind=8),INTENT (IN OUT) :: resid(:,:),& !internal nodal forces
                                 energ(:),  & !bottom surface coordinates
                                 coorb(:,:),& !bottom surface coordinates
                                 coort(:,:)   !top surface coordinates
 TYPE (ele07), POINTER  :: heade              !pointer to first element
 LOGICAL, INTENT(IN) :: zigzag

 INTEGER (kind=4) ielem,nvarl,isec,l,n,matty,nlayr,osec,np,mtype,i,m,ni,ndofe,secty
 LOGICAL ::  plast,quad,mems(3)

 REAL(kind=8) :: thick,thico,tgp(3),ds,sf(3)
 REAL(kind=8), ALLOCATABLE ::  bmsig(:,:), x(:,:),t(:,:),lb(:,:,:),psi(:,:)
 TYPE(ele07), POINTER :: e
 TYPE(section), POINTER :: sec
 TYPE(mater), POINTER :: mat

 IF( ngaus /= 1 )THEN
   ALLOCATE( x(3,nnode),t(3,nnode),lb(3,2,nnode))
   IF(zigzag) ALLOCATE( psi(2,6))
 ELSE
   ALLOCATE( t(3,3),lb(3,2,3))
   IF( nnode == 6 )THEN
     quad = .FALSE.
     ALLOCATE( x(3,3))
   ELSE
     quad = .TRUE.
     ALLOCATE( x(3,6))
   END IF
   IF(zigzag) ALLOCATE( psi(2,3))
 END IF

 !***  loop over all the elements

 osec = -1
 e => heade
 mtype = 4-stype
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
     IF(matty == 0)THEN   !layered isotropic
       np = 12
       ds = sec%rprop(9)*thico  !shear stabilization modulus
     ELSE                !in therms of forces and moments
       np = 21
       ds = sec%rprop(11)       !shear stabilization modulus
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
     IF( ngaus /= 1 )THEN
       nvarl = nnode*ndofe
       IF( ALLOCATED(bmsig) )DEALLOCATE (bmsig)
       ALLOCATE( bmsig(ndofe,nnode))
     ELSE
       IF( quad ) THEN
         ALLOCATE( bmsig(ndofe,6))
       ELSE
         ALLOCATE( bmsig(ndofe,3))
       END IF
     END IF
   END IF
   !     gathers the actual coordinates and directors
   thick = thico
   bmsig = 0d0
   IF( ngaus /= 1 )THEN
     x  = coora(:,e%lnods)
     lb(:,1,:) = euler(1:3,e%lnods)
     lb(:,2,:) = euler(4:6,e%lnods)
     t  = euler(7:9,e%lnods)
     IF( ASSOCIATED(e%nangl) )THEN
       DO n=1,nnode
         l = e%lnods(n)
         psi(1,n) = psia(1,l)*e%nangl(1,n)+psia(2,l)*e%nangl(2,n)
         psi(2,n) = psia(2,l)*e%nangl(1,n)-psia(1,l)*e%nangl(2,n)
       END DO
     END IF
     CALL intrf7(nvarl,plast,ndofe,nstre,shape,posgp,ngaus,ngamm,nnode,        &
                 stype,matty,nlayr,sec%rprop(6:),sec%rprop(np:),sf,            &
                 x,t,lb,psi,e%nangl,e%jas,                                     &
                 e%dvolu,e%cartd,e%jacin,                          &
                 e%stra0,e%gamm0,e%strsg,e%qbar,e%beta,            &
                 e%ambda,e%ehist,e%strap,e%stres,                  &
                 bmsig,                                            &
                 dn(:,:,mtype:4),ap1,energ,thick,istop,ielem)
     DO n=1,nnode
       l = e%lnods(n)
       resid(1:5,l) = resid(1:5,l) + bmsig(1:5,n)
       IF(zigzag)resid(7:8,l) = resid(7:8,l) + bmsig(6:7,n)
     END DO
     !WRITE(58,"(7e12.5)")bmsig
     !WRITE(58,"(7e12.5)")
     IF(top .OR. bottom)THEN
       DO n=1,nnode
         l = e%lnods(n)
         IF( top  )coort(:,l) = coort(:,l) + euler(7:9,l)*thick
         IF(bottom)coorb(:,l) = coorb(:,l) - euler(7:9,l)*thick
         ifact(l) = ifact(l) + 2
       END DO
     END IF
   ELSE
     x  = coora(:,e%lnods(1:3))
     lb(:,1,:) = euler(1:3,e%lnods(4:6))
     lb(:,2,:) = euler(4:6,e%lnods(4:6))
     t  = euler(7:9,e%lnods(4:6))
     IF( ASSOCIATED(e%nangl) )THEN
       DO n=1,3
         l = e%lnods(n)
         psi(1,n) = psia(1,l)*e%nangl(1,n)+psia(2,l)*e%nangl(2,n)
         psi(2,n) = psia(2,l)*e%nangl(1,n)-psia(1,l)*e%nangl(2,n)
       END DO
     END IF
     IF(quad)THEN
       DO m =1,3                       !for the side nodes in the patch
         i = m+3                       !associated node
         ni = e%lnods(i+3)               !node number
         IF(ni > 0) THEN
           x(1:3,i) = coora(1:3,ni)    !get coordinates of adjacent node
           mems(m) = .TRUE.               !side element exists
         ELSE
           mems(m) = .FALSE.              !side element does not exist
           x(1:3,i) = 0d0            !null values (unnecessary perhaps)
         END IF
       END DO
     END IF
     CALL intr07(ndofe,nstre,ngaus,ngamm,ap1(:,:,1),plast,                 &
                 matty,nlayr,sec%rprop(6:),sec%rprop(np:),sf,              &
                 x,t,lb,tgp,e%dvolu(1),e%cartd(:,:,1),e%jacin(:,:,1),e%cd, &
                 mems,e%stra0(:,1),e%gamm0,e%strsg(:,1),e%qbar,psi,e%nangl,&
                 e%jas,e%ambda(:,1),e%ehist,e%strap,e%stres,               &
                 bmsig,energ,thick,istop,ielem,stabq,quad,ds)
     DO n=1,3
       l = e%lnods(n)
       m = e%lnods(n+3)
       resid(1:3,l) = resid(1:3,l) + bmsig(1:3,n)
       resid(4:5,m) = resid(4:5,m) + bmsig(4:5,n)
       IF(quad)THEN
         ni = e%lnods(n+6)               !node number
         IF(ni > 0) resid(1:3,ni) = resid(1:3,ni) + bmsig(1:3,n+3)
       END IF
       IF(zigzag)resid(7:8,l) = resid(7:8,l) + bmsig(6:7,n)
     END DO
     !WRITE(58,"(10e12.4)")e%strsg(1:6,1),e%strsg(9:12,1)
     !WRITE(58,"(4e12.4)")e%strsg(7:8,1),e%strsg(13:14,1)
     !WRITE(58,"(7e12.5)")bmsig
     IF(top .OR. bottom)THEN
       DO n=1,3
         l = e%lnods(n)
         IF( top  )coort(:,l) = coort(:,l) + tgp*thick
         IF(bottom)coorb(:,l) = coorb(:,l) - tgp*thick
         ifact(l) = ifact(l) + 2
       END DO
     END IF
   END IF
   e => e%next
 END DO
 DEALLOCATE( bmsig,x,t,lb)
 !WRITE(58,"(8e12.5)")resid
 RETURN
 END SUBROUTINE resvp7

