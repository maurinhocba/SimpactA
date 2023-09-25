 SUBROUTINE tcsurf ( )

 ! apply convection and radiation contributions

 USE heat_db
 USE ctrl_db, ONLY : ndime,ntype,ndoft,istep,ttime,tref
 USE npo_db, ONLY : iftmp,coora,tempe,tresi
 USE curv_db
 ! local variables
 IMPLICIT NONE
 TYPE (srf_con), POINTER :: s
 INTEGER(kind=4) :: i,ii,j,l,n,m,mdof,ndc,ndt,ndb
 REAL(kind=8) :: t1(ndime),t2(ndime),t3(ndime),hconc,hradc,hcont,hconb, &
                 hradt,hradb,thc,trc,tht,trt,thb,trb
 LOGICAL :: hc,rc,ht,rt,hb,rb
 REAL (kind=8) :: functs,k,f,avt,avtt,avtb,fac,dlum,dcon,rad
 REAL (kind=8),ALLOCATABLE :: x(:,:),rl(:,:) !x(ndime,nn),rl(nn,ndoft)
 REAL (kind=8),PARAMETER   :: twopi=6.283185307209586

 !___________________________________

 s => headsc !point to first surface

 DO m=1,ncsur  ! loop over each surface
   s%carea = s%carea .OR. MOD(istep,s%afreq) == 1  !compute area?
   ALLOCATE(x(ndime,s%nn))
   IF( s%carea )THEN !if TRUE
     DO i=1,s%nelem  !loop over each segment
       x(:,1:s%nn) = coora(:,s%lnods(:,i))         !coordinates
       IF( ndime == 2 )THEN
         t1 = x(:,2) - x(:,1)
         s%area(i) = SQRT(t1(1)*t1(1)+t1(2)*t1(2)) !segment length
       ELSE
         t1 = x(:,2) - x(:,1)                      !first side
         t2 = x(:,3) - x(:,1)                      !second side
         CALL vecpro(t1,t2,t3)                     !element normal
         s%area(i) = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))/2d0
       END IF
     END DO
     s%carea = .FALSE.  !modify flag
   END IF
   !***  get Convection/radiation factors and ambient/radiation temperatures ***
   !
   ! for SOLID problems (or shells with only one central thermal DOF)
   ndc = 1  !Associated dof with central surface
   hconc = s%hf(1)    !convection factor
   hc = hconc > 0d0 !convection on the surface?
   IF(hc)THEN       !if TRUE
     thc = s%tf(1)    !constant ambient temperature
     IF( thc == 0d0 .AND. s%tcur(1) > 0 )thc = functs (s%tcur(1),ttime)
   END IF
   hradc = s%hf(4)    !radiation factor
   rc = hradc > 0d0 !radiation on this surface?
   IF(rc)THEN       !if TRUE
     trc = s%tf(4)    !constant radiation temperature
     IF( trc == 0d0 .AND. s%tcur(4) > 0 )trc = functs (s%tcur(4),ttime)
     trc = (trc+tref)**4
   END IF
   ndt = ndc !
   tht = thc !for convenience
   trt = trc !
   mdof = 1
   ! for SHELLS with 2 or 3 thermal DOFs
   IF(s%mdof > 1)THEN
     ndt = s%mdof      !Associated dof with top surface
     ndb = s%mdof-1    !Associated dof with bottom surface
     hcont = s%hf(2)     !convection factor for Top surface
     ht = hcont > 0    !convection on top surface?
     IF(ht)THEN        !if TRUE
       tht = s%tf(2)     !constant ambient temperature
       IF( tht == 0d0 .AND. s%tcur(2) > 0 )tht = functs (s%tcur(2),ttime)
     END IF
     hconb = s%hf(3)     !convection factor for bottom surface
     hb = hconb > 0    !convection on bottom surface?
     IF(hb)THEN        !if TRUE
       thb = s%tf(3)     !constant ambient temperature
       IF( thb == 0d0 .AND. s%tcur(3) > 0 )thb = functs (s%tcur(3),ttime)
     END IF
     hradt = s%hf(5)     !radiation factor on top surface
     rt = hradt > 0    !radiation on top surface?
     IF(rt)THEN        !if TRUE
       trt = s%tf(5)     !constant radiation temperature
       IF( trt == 0d0 .AND. s%tcur(5) > 0 )trt = functs (s%tcur(5),ttime)
       trt = (trt+tref)**4
     END IF
     hradb = s%hf(6)     !radiation factor on bottom surface
     rb = hradb > 0    !radiation on bottom surface?
     IF(rb)THEN        !if TRUE
       trb = s%tf(6)     !constant radiation temperature
       IF( trb == 0d0 .AND. s%tcur(6) > 0 )trb = functs (s%tcur(6),ttime)
       trb = (trb+tref)**4
     END IF
     mdof = ndb
   END IF
   !
   !  loop over each segment on the surface
   !  below code is valid ONLY for NDIME=3, NN=3 (triangles)
   ALLOCATE(rl(s%nn,ndoft))
   DO i=1,s%nelem
     rl = 0d0      !initializes residual contribution

     IF(s%solid)THEN !for SOLID problems
       t1  = tempe(ndc,s%lnods(:,i)) !get temperatures on the surface
       avt = SUM(t1)                 !average temperature
       fac  = 1d0 
       IF( ndime == 2 )THEN !for segments       
         IF( ntype == 3 )THEN
           x(:,1:s%nn) = coora(:,s%lnods(:,i)) !coordinates of segment
           rad = (x(1,1) + x(1,2))/2d0         !mid segment radius         
           fac = twopi*rad !factor for axilsymmetric problem (2 pi R)
         END IF  
         dcon =  6d0
         dlum =  2d0
       ELSE                 !for triangles
         dcon = 12d0
         dlum =  3d0
       END IF
       ! convection contribution (h)
       IF( hc )THEN                !if convection is defined
         k  = s%area(i)*fac*hcont/dcon      !2D: hLt/6 or 3D: hA/12
         t3 = (avt + t1(:))*k               !
         f  = -s%area(i)*fac*hcont/dlum*thc !2D: -hLt/2 or 3D: -hA/3 Tenv
         rl(:,ndc) = rl(:,ndc) + t3 + f
       END IF
       ! radiation contribution (r)
       IF( rc )THEN                !if radiation is defined
         k  = s%area(i)*fac*hradt/dlum   !2D: rLt/3 or 3D: hA/3
         f  = ((avt/dlum)**4 - trc)*k    !
         rl(:,ndc) = rl(:,ndc) + f
       END IF

     ELSE !for SHELL problems
       IF( hc .OR. rc .OR. ht .OR. rt )THEN !top or central surface
         t1 = tempe(ndt,s%lnods(:,i)) !get temperatures on the surface
         avtt = SUM(t1)                !average temperature
       END IF
       IF( hb .OR. rb ) THEN       !bottom surface
         t2 = tempe(ndb,s%lnods(:,i)) !get temperatures on the surface
         avtb = SUM(t2)                !average temperature
       END IF
       ! convection contribution (h)
       IF( hc .OR. ht )THEN         !top surface
         k  = s%area(i)*hcont/12d0         ! hA/12
         t3 = (avtt + t1(:))*k             !
         f  = -s%area(i)*hcont/3d0*tht     ! -hA/3 Tenv
         rl(:,ndt) = rl(:,ndt) + t3 + f
       END IF
       IF( hb )THEN                 !bottom surface
         k  = s%area(i)*hconb/12d0         ! hA/12
         t3 = (avtb + t2(:))*k             !
         f  = -s%area(i)*hconb/3d0*thb     ! -hA/3 Tenv
         rl(:,ndb) = rl(:,ndb) + t3 + f
       END IF
       ! radiation contribution (r)
       IF( rc .OR. rt )THEN         !top surface
         k  = s%area(i)*hradt/3d0          ! hA/3
         f  = ((avtt/3d0)**4 - trt)*k
         rl(:,ndt) = rl(:,ndt) + f
       END IF
       IF( rb )THEN                 !bottom surface
         k  = s%area(i)*hradb/3d0           ! hA/3
         f  = ((avtb/3d0)**4 - trb)*k
         rl(:,ndb) = rl(:,ndb) + f
       END IF

     END IF

     ! sum contributions on residual
     DO ii=1,s%nn          !for each node in the element
       n = s%lnods(ii,i)   !node
       DO j=mdof,ndoft     !
         l = iftmp(j,n)
         IF( l > 0 ) tresi(l) = tresi(l) + rl(ii,j) !sum internal force
       END DO
     END DO
   END DO
   s => s%next
 END DO
 RETURN
 END SUBROUTINE tcsurf
