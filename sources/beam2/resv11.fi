 SUBROUTINE resv11(nele0,nbn,nstre,ntype,ngaus,head0,nhead,sh,iffix,    &
                   coora,resid,ierr,ttime,bottom,top,coorb,coort,ifact, &
                   stabs,stint)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma) for beam elements
 !
 !********************************************************************
 !$ USE omp_lib
 !USE npoi_db, ONLY : tempe, itemp
 IMPLICIT NONE
 LOGICAL, PARAMETER :: logst = .TRUE.  !include as data latter
 INTEGER (kind=4), INTENT(IN) :: nele0,nbn,nstre,ntype,ngaus,iffix(:)
 INTEGER (kind=4), INTENT(IN OUT) :: ierr,ifact(:)
 REAL (kind=8), INTENT(IN) :: coora(:,:), & !present coordinates
                              sh(2,2),    & !shape functions
                              stabs,      & !stabilization factor
                              ttime         !present total time
 REAL (kind=8), INTENT(IN OUT) :: resid(:,:),coort(:,:),coorb(:,:),stint(:,:)
 LOGICAL, INTENT(IN) :: bottom,top
 TYPE(ele11), POINTER :: head0
 TYPE(nodei), POINTER :: nhead

 REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0

 INTEGER (kind=4) iel0,  & !first element of the loop - 1
                  nelem, & !number of elements at each processor
                  ielem, & !element number
                  isec,  & !associated material
                  i,j,k,l,g,nm,& !different indexes
                  ns,    & !number of strains and stresses
                  nlayr, & !number of layers
                  mtype, & !associated material type
                  fixit(2), ln(4), pp(5,2), & !case
                  numpt, & !number of points in curve
                  osec     !associated material of previous element
 REAL (kind=8) :: y0(2),yi(2), & !nodal coordinates
                  t(2),a0,     & !tangent vector, initial angle
                  li,dvolu,    & !element length
                  n(2),        & !element normal
                  beta,        & !average rotation
                  x(2,nnode),  & !element nodal coordinates
                  pff(2,2),    & !first fundamental form
                  curv(2,2),   & !curvatures
                  ro(2,2),     & !initial radius and curvatures
                  stres(2),    & !stresses (different measures)
                  css(2),      & !strains
                  stran(2),    & !strains
                  force(2),    & !integrated forces
                  momen(2),    & !moments
                  bmem(4,2),   & !membrane B matrix
                  bfle(2,4,2), & !bending B matrix (standard)
                  lb(3),lc(3), & !eigenvalues
                  thick,width, & !thickness and width (original)
                  thnew,       & !present thickness
                  zk,          & !distance to mid surface
                  aux,auxi(6), & !auxiliar value
                  c(3),cs,gs,  & !Elastic constitutive matrix for plane stress
                  propi(7),    & !Plastic Material properties
                  chi(12),     & !Hill 48 coefficients or rubber material
                  rl(2,4),     & !nodal internal equivalent forces
                  deatht,      & !end time for plasticity
                  min_tr,      & !Minimum Allowable Thickness ratio
                  max_tr         !Maximum Allowable Thickness ratio
 TYPE(ele11), POINTER :: e,    &  !pointer to present element
                         heade    !pointer to first element in the loop
 TYPE(section), POINTER :: sec   !pointer to section
 TYPE(mater), POINTER :: mat     !pointer to material
 TYPE(nodei), POINTER :: ni,nd   !pointers to left and rigth branching nodes
 LOGICAL :: fix(2),  &! TRUE if side is clamped or symmetry line
            newmt,   &! TRUE if material constant computation necessary
            natst,   &! TRUE for large strain analysis
            elast,   &! TRUE if material is strictly elastic
            shell,   &! TRUE if bending included
            simpl,   &! FALSE if one of the nodes is branched
            siden(2),&! TRUE if the node is on a side
            kinhd,   &! TRUE if kinematic hardening is to be considered
            plast,   &! TRUE if plasticity is to be considered
            pflag     ! TRUE if plastic flow


 !REAL (kind=8), ALLOCATABLE :: cm(:,:), prop(:,:), volfr(:), rr(:,:)
 REAL (kind=8), POINTER :: val(:,:)

 CHARACTER(len=1 ) letra  !for console input

 ! Gauss points throught the thickness
 !REAL (kind=8), ALLOCATABLE :: thf(:),wei(:)   !Gass point position and weigth (TTT)
 REAL (kind=8) :: thf(21),wei(21)   !Gass point position and weigth (TTT)

 INTERFACE
   INCLUDE 'gaussq.h'
   INCLUDE 'core06ps.h'
   INCLUDE 'rubberps.h'
   INCLUDE 'stre02_1.h'
   INCLUDE 'stre02_6.h'
 END INTERFACE

 !Auxiliar variables for parallel code
 !!$ INTEGER(kind=4), PARAMETER :: maxel=34816,minel=32768  !8 threads
 !!$ INTEGER(kind=4), PARAMETER :: maxel=18432,minel=16384  !4 threads
 !$ INTEGER(kind=4), PARAMETER :: maxel=10240,minel= 8192   !2 threads
 !$ INTEGER(kind=4):: nloop,iloop,irest,ii    ! number of loops, present loop, and last group
 INTEGER(kind=4):: lim_inf,  & ! first element for each processor
                   lim_sup,  & ! last element for each processor
                   num_chk     ! number of elements chunks

 !$ INTEGER(kind=4),PARAMETER:: MAXLOCK=2048              !number of keys
 !$ INTEGER(kind=4),SAVE     :: lock1(0:MAXLOCK-1)        !keys
 !$ LOGICAL,SAVE             :: inilock=.FALSE.           !initialization flag

 !Initialize exclusion variables for parallel code (only first time)
 !$ IF (.NOT.inilock) THEN                   !first time only
 !$   DO ii=0,MAXLOCK-1                     !for each key
 !$     CALL omp_init_lock(lock1(ii))          !unlock
 !$     PRINT *,ii,lock1(ii)
 !$   END DO
 !$   inilock = .TRUE.                         !change flag
 !$   print *,'first time ok'
 !$ END IF
 !---------------- FIRST deal with branching nodes (non-parelellized)
 ni => nhead     !point to the first branching node
 DO i=1,nbn           !for each branching node
   y0 = coora(:,ni%lnods(0))          !branching node coordinates
   beta = 0d0                         !initializes average rotation
   ni%bb(:,0) = 0d0                   !initializes B matrix
   DO j=1,ni%nn                       !for each connected node
     yi = coora(:,ni%lnods(j))        !coordinates of the node
     t  = yi - y0                     !oriented vector 0-i
     li = SQRT(DOT_PRODUCT(t,t))      !segment length
     t = t/li                         !unit tangent vector
     n = (/ -t(2),t(1) /)*ni%fc(j)/li      !unit normal vector, divided by L, times Ri
     beta = beta + ni%fc(j)*ATAN3(t(2),t(1),ni%alph0(j)) !average rotation
     ni%bb(:,0) = ni%bb(:,0) - n      !average rotation B matrix : at node
     ni%bb(:,j) = n                   !                          : at opposite node
   END DO
   ni%ba = beta              !average rotation value
   ni => ni%next             !point to next node
 END DO
 !----------------

 ns   = nstre/2  !number of stress components
 k = 1  !generate array of pointers to STINT
 DO i=1,ngaus
   pp(1,i) = k
   pp(2,i) = pp(1,i)+ns-1
   pp(3,i) = pp(2,i)+1
   pp(4,i) = pp(2,i)+ns
   pp(5,i) = pp(4,i)+1
   k = pp(5,i) + 1
 END DO

 heade => head0
 nelem = nele0
 iel0 = 0
 !  check maximum number of elements to be processed
 !$ IF( nelem > maxel )THEN
 !$   nloop = nelem/minel
 !$   irest = MOD(nelem,minel)
 !$  IF( irest >= maxlock ) nloop = nloop + 1
 !$  nelem = minel
 !$ ELSE
 !$  nloop = 1
 !$ END IF
 !$ DO iloop=1,nloop

 !initialize variables
 num_chk =  1    !only one processor or thread
 lim_inf =  1    !first element
 lim_sup = nelem !last element

 !$OMP  PARALLEL                                                                           &
 !$OMP& SHARED (num_chk, lock1, psecs, nelem, nndpd,                                       &
 !$OMP&         nstre,ntype,ngaus,iffix, ierr, ifact, coora, sh, stabs, ttime, ns,         &
 !$OMP&         bottom,top, nhead, resid, coort, coorb, stint, letra, iel0, pp, heade)     &
 !$OMP& PRIVATE(lim_inf, lim_sup,ii,                                                       &
 !$OMP&         ielem, isec, i,j,k,l,g,nm, nlayr, mtype, fixit, ln,  numpt, osec,     &
 !$OMP&         y0,yi, t,a0, li,dvolu, n, x, pff, curv, ro, stres, css,               &
 !$OMP&         stran, force, momen, bmem, bfle, lb,lc, thick,width, thnew, zk,       &
 !$OMP&         aux,auxi, c,cs,gs, propi, chi, rl, deatht, min_tr,                    &
 !$OMP&         max_tr, e, sec, mat, ni,nd, fix, newmt, natst, elast, shell,               &
 !$OMP&         simpl, siden, kinhd, plast, pflag, val, thf,wei)
 !$ IF( omp_get_num_threads() > 1 )THEN
   !$ num_chk = (nelem+omp_get_num_threads()-1)/omp_get_num_threads()
   !$ lim_inf = MIN(nelem,1+omp_get_thread_num()*num_chk) + iel0    !Find the first element in the list
   !$ lim_sup = MIN(nelem,num_chk*(1+omp_get_thread_num())) + iel0  !Find the last element in the list
 !$ END IF

 e => heade                          !point to head of the list
 DO ielem=iel0+1,lim_inf-1             !Loop over the list to find the first element in the list
   e => e%next                       !Point the first element of the list
 END DO

 !***  loop over all the elements
 auxi = 0d0
 osec = -1       !initializes

 ELEMS : DO ielem=lim_inf,lim_sup                   !loop over each element
   isec  = e%matno                  !section
   newmt = osec /= isec             !same section than previous ?
   IF(newmt)THEN                    !if new section/material, get properties
     sec => psecs(isec)%p           !point to new section
     mat  => sec%mtbas              !point to associated material
     osec  = isec                   !keep section lable
     nlayr = sec%iprop(1)           !number of layers
     !IF( ALLOCATED(thf)) DEALLOCATE(thf,wei)
     !ALLOCATE(thf(nlayr),wei(nlayr))
     shell = nlayr > 1              !membrane or shell
     mtype = mat%mtype              !type of base material
     thick = sec%rprop(1)           !original thickness
     width = sec%rprop(10)          !original width
     min_tr= sec%rprop(2)           !Minimum Thickness ratio
     max_tr= sec%rprop(3)           !Maximum Thickness ratio
     natst = logst .AND. mtype /= 6  !use log strains
     elast = mat%matdef(3) == 1      !elastic?
     CALL gaussq(nlayr,thf(1:nlayr),wei(1:nlayr)) !integration points through the thickness
     thf(1:nlayr) = thf(1:nlayr)/2d0  !positions
     wei(1:nlayr) = wei(1:nlayr)/2d0  !weights
     !alpha = mat%prope(6)             !Thermical Dilatation Coeff
     !IF( .NOT.itemp ) alpha = 0d0

     SELECT CASE(mtype)
     CASE (1)                     !elastic isotropic - plastic orthotropic
       IF( .NOT.elast )THEN
         propi(1:7) = mat%propp(1:7)       ! isotropic & kinematic hardening parameters
         propi(5) = REAL( mat%matdef(4),8) ! isotropic hardening model
         chi    = mat%propp(16:27)         ! hill coefficients
         IF( mat%matdef(8) == 1  ) chi(1:6) = chi(7:12)  !associative
         kinhd =  mat%matdef(5) > 1 .AND. ntype == 1 !kinematic hardening
         deatht = mat%propp(5)             !end of plasticity
         IF (mat%matdef(4) == 5 )THEN
           val => mat%chead%val
           numpt = mat%chead%np
         ELSE
           NULLIFY (val)
           numpt = 0
         END IF
         plast = deatht > ttime .AND. propi(1) > 0  !consider plasticity ?
       ELSE
         plast = .FALSE.
       END IF
       c(1:3) = mat%prope(7:9)           ! plane stress elasticity matrix
       gs     = sec%rprop(13)
     CASE (5)                     !elastic orthotropic
       ! plane stress elasticity matrix
       c(1) = mat%props(1)
       c(2) = mat%props(2)
       c(3) = mat%props(5)
       elast = .TRUE.
       plast = .FALSE.                        !consider plasticity ?
     CASE (6)                     !hyperelastic isotropic
       ! plane stress elasticity matrix
       chi = mat%prope(7:18)      !elastic coefficients
       IF( .NOT.elast )THEN
         propi(1:2) = mat%propp(1:2)       ! isotropic hardening parameters
         plast = .TRUE.
       ELSE
         plast = .FALSE.                        !consider plasticity ?
       END IF
     END SELECT
     IF( shell )THEN
       cs = sec%rprop(12)/12d0*stabs     !3*2*2 stabilization
     ELSE
       cs = 0d0
     END IF
     newmt = .FALSE.                        !
   END IF

   ! gathers the initial coordinates and fixities

   simpl = .TRUE.
   siden = .FALSE.
   NULLIFY( ni,nd )          !nullify pointers to left and rigth branching nodes
   fixit = 0                 !Both sides exist and are smooth
   ln = e%lnods(1:4)         !copy nodes
   fix = .FALSE.             !initializes
   DO j=1,nnode              !nnode = 4
     k = ln(j)          !node value
     SELECT CASE ( k )       !
     CASE (:-1)        !find pointers to intersection node
       simpl = .FALSE.
       IF( j == 1 )THEN!first node
         ni => nhead   !point to first node
         DO i=1,-k-1   !loop
           ni => ni%next
         END DO
         fixit(1) = 3
       ELSE            !second node
         nd => nhead   !point to first node
         DO i=1,-k-1   !loop
           nd => nd%next
         END DO
         fixit(2) = 3  !branching node
       END IF

     CASE (0)          !end node
       IF( j == 1 )THEN      !first node
         fix(1) = iffix(ln(2)) /= 0
         IF(iffix(ln(2)) < 0) THEN  !node on a side
           simpl = .FALSE.
           a0 = e%cab(1,1)
           nm = -iffix(ln(2))  !compute coordinates of the right node
           x(:,1) =  (/ coora(2,nndpd(3,nm))-coora(2,nndpd(2,nm)), &
                        coora(1,nndpd(2,nm))-coora(1,nndpd(3,nm)) /)
           li = DOT_PRODUCT(x(:,1),x(:,1))
           x(:,1) = x(:,1)/li
           e%cab(1,1) = ATAN2(x(2,1),x(1,1)) + a0       !a0 initial angle of solid
           fixit(1) = 2      !clamped
           siden(1) = .TRUE.
         ELSE
           IF( fix(1) )THEN
             fixit(1) = 2      !clamped
           ELSE
             fixit(1) = 1      !free
           END IF
         END IF
       ELSE                  !second node
         fix(2) = iffix(ln(3)) /= 0
         IF(iffix(ln(3)) < 0) THEN  !node on a side
           simpl = .FALSE.
           a0 = e%cab(1,2)
           nm = -iffix(ln(3))  !compute coordinates of the right node
           x(:,4) = (/ coora(2,nndpd(3,nm))-coora(2,nndpd(2,nm)), &
                       coora(1,nndpd(2,nm))-coora(1,nndpd(3,nm)) /)
           li = DOT_PRODUCT(x(:,4),x(:,4))
           x(:,4) = x(:,4)/li
           e%cab(1,2) = ATAN2(x(2,4),x(1,4)) + a0       !a0 initial angle of solid
           fixit(2) =  2     !clamped
           siden(2) = .TRUE.
         ELSE
           IF( fix(2) )THEN
             fixit(2) =  2     !clamped
           ELSE
             fixit(2) =  1     !free
           END IF
         END IF
       END IF

     CASE (1: )        !inner node
       x(:,j) = coora(:,k)   !get coordinates
     END SELECT

   END DO
   ! compute shear force
   IF( ngaus == 1 )THEN
     SELECT CASE (e%lnods(5))    !left element
     CASE (0)      !simple supported
       bmem(1,1) =  e%l0(2)
       bmem(1,2) = -stint(pp(3,1),ielem)
     CASE (-1)     !clamped
       bmem(1,1) =  e%l0(2)
       bmem(1,2) =  stint(pp(3,1),ielem)
     CASE (-2)     !branching node
       bmem(1,1) =  e%l0(3)
       bmem(1,2) =  stint(pp(3,1),ielem)
     CASE DEFAULT
       bmem(1,1) = e%l0(1)
       bmem(1,2) = stint(pp(3,1),e%lnods(5))
     END SELECT
     bmem(2,1) = e%l0(2)
     bmem(2,2) = stint(pp(3,1),ielem)
     SELECT CASE (e%lnods(6))
     CASE (0)      !simple supported
       bmem(3,1) =  e%l0(2)
       bmem(3,2) = -stint(pp(3,1),ielem)
     CASE (-1)     !clamped
       bmem(3,1) =  e%l0(2)
       bmem(3,2) =  stint(pp(3,1),ielem)
     CASE (-2)     !branching node
       bmem(3,1) =  e%l0(1)
       bmem(3,2) =  stint(pp(3,1),ielem)
     CASE DEFAULT
       bmem(3,1) =  e%l0(3)
       bmem(3,2) =  stint(pp(3,1),e%lnods(6))
     END SELECT

     lb(1) = -1d0/(bmem(1,1)+bmem(2,1))
     lb(3) =  1d0/(bmem(2,1)+bmem(3,1))
     lb(2) = -lb(1)-lb(3)
     stint(pp(5,1),ielem) = lb(1)*bmem(1,2)+lb(2)*bmem(2,2)+lb(3)*bmem(3,2)
   ELSE
     aux  =  1.732050808d0/ e%l0(2)
     stint(pp(5,1:2),ielem) = (stint(pp(3,2),ielem) - stint(pp(3,1),ielem))*aux
   END IF
   IF( ntype == 3) ro(1:2,1:ngaus) = e%ro(1:2,1:ngaus)    !for axilsymmetric get initial R
   CALL stra11(ntype,ngaus,x,ni,nd,fix,ln,   &   !compute stretchings and curvatures
               e%l0,e%cab,ro,e%lambd,pff,curv,t,n,sh)
   IF( ASSOCIATED( e%cur ) )THEN
     IF( ngaus == 1 )THEN
       e%cur = (/ curv(1,1)-curv(1,2), curv(1,1)+curv(1,2) /)/2d0
     ELSE
       e%cur = curv(1,:) /2d0
     END IF
   END IF
   !check maximum and miminum thickness ratio
   IF( e%lambd(2) < min_tr .OR. e%lambd(2) > max_tr)THEN      !check
     WRITE(*,"(' thickness ratio limit exceeded, TH_RATIO: ',f5.3, &
           &  /,' Do you want to continue [Y/N]')")e%lambd(2)
     DO                            !pause until answer
       READ(*,'(a1)')letra         !read the answer
       SELECT CASE (letra)         !process answer
       CASE ( 'y', 'Y')         !continue
         IF( e%lambd(2) < min_tr )THEN    !minimum exceeded
           WRITE(*,"(' New value for minimum thickness ratio')")
           READ(*,*) min_tr        !read value
           sec%rprop(2) = min_tr   !Store Minimum Thickness ratio
         ELSE                      !maximum exceeded
           WRITE(*,"(' New value for maximum thickness ratio')")
           READ(*,*) max_tr        !read value
           sec%rprop(3) = max_tr   !Store Maximum Thickness ratio
         END IF
       CASE ( 'n', 'N')         !stop process
         ierr = 1                  !modify flag
       CASE DEFAULT
         CYCLE                     !answer not valid, ask again
       END SELECT
       EXIT
     END DO
   END IF
   ! consider thermal strains
   !IF( alpha > 0 )THEN     ! compute
   !  t0 = 1d0+alpha*SUM(tempe(e%lnods(2:3),1))
   !  e%lamb = e%lamb*t0        ! modify thickness ratio
   !END IF

   thnew = thick*e%lambd(2)            !thickness at new configuration

   !IF( alpha > 0 )t1 = alpha*SUM(tempe(e%lnods(2:3),2))*2d0/thnew

   IF( bottom .OR. top)THEN       !bottom or top surfaces required
     DO i=1,2                   !for each side
       j = ln(i+1)         !node
       !$ CALL omp_set_lock(lock1(MOD(j,MAXLOCK)))    !Start blocking acess in shared variables
       IF( fix(i) )THEN         !if side clamped
         yi = (/ -SIN(e%cab(1,i)),COS(e%cab(1,i)) /) !normal at plane of symmetry
         IF(bottom) coorb(1:2,j) = coorb(1:2,j) - yi(1:2)*thnew !bottom surface
         IF(top   ) coort(1:2,j) = coort(1:2,j) + yi(1:2)*thnew !top surface
       ELSE
         IF(bottom) coorb(1:2,j) = coorb(1:2,j) - n(1:2)*thnew !bottom surface
         IF(top   ) coort(1:2,j) = coort(1:2,j) + n(1:2)*thnew !top surface
       END IF
       !$ CALL omp_unset_lock(lock1(MOD(j,MAXLOCK)))  !End blocking acess in shared variables
       ifact(j) = ifact(j) + 2  !increase weigthing factor
     END DO
   END IF

   rl = 0d0                             !initializes internal forces (4-nodes only)
   DO g=1,ngaus

     force = 0d0                        !initializes Integrated forces
     momen = 0d0                        !initializes Integrated moments
     !Trought the thickness (t.t.t.) integration loop

     DO l=1,nlayr                        !for each layer
       zk = thf(l)*thnew                    !Z coordinate
       css(1:ns) = pff(1:ns,g)+curv(1:ns,g)*zk  !layer U^2
       !IF( alpha > 0 )THEN                 !consider thermal strains
       !  j0 = (t0 + zk*t1)**(2d0/3d0)      !modifies jacobian due to temperature
       !  stran = stran/j0                  !modifies C
       !END IF
       IF( mtype /= 6 )THEN
         IF( natst )THEN !Hencky strains
           DO i=1,ns
             stran(i) = LOG(css(i))/2d0
           END DO
         ELSE            !Green Lagrange Strains
           DO i=1,ns
             stran(i) = (css(i)-1d0)/2d0
           END DO
         END IF
       END IF
       SELECT CASE ( mtype )
       CASE (1)  !one component material
         IF( ntype == 1 )THEN
           !!(1:6)stre0,stres,?,ep,p,eps
           IF(plast)THEN
             auxi(4) = e%gausv(1,l,g)   !plastic strain
             auxi(6) = e%gausv(2,l,g)   !equivalent plastic strain
             IF( kinhd) auxi(5) = e%gausv(3,l,g)   !back stress
           END IF
           CALL stre02_1(mat,stran(1),auxi,.TRUE.,aux)
           stres(1) = auxi(2)
           IF(plast) THEN
             e%gausv(1,l,g) = auxi(4)  !plastic strain
             e%gausv(2,l,g) = auxi(6)  !equivalent plastic strain
             IF( kinhd) e%gausv(3,l,g) = auxi(5)  !back stress
           END IF
         ELSE
           CALL streb2(stran,stres,c,propi,chi,e%gausv(:,l,g), &
                     ierr,plast,elast,numpt,val)
         END IF
         IF (ierr == 1) CYCLE ELEMS       ! if there is an error CYCLE TO END (RETURN or EXIT are not allowed in parallel)
       CASE (5)
         !elastic (trial) stresses (SHELLs Only)
         stres(1)= c(1)*stran(1)+c(2)*stran(2)
         stres(2)= c(2)*stran(1)+c(3)*stran(2)

       CASE (6)
         lb(1) = SQRT(css(1))           !stretch
         IF( ntype == 1 ) THEN
           IF(plast)THEN
             auxi(6) = e%gausv(2,l,g)   !equivalent plastic strain
             CALL stre02_6(mat,lb(1),auxi,aux)
             stres(1) = auxi(2)/css(1)  !second Piol Kirchhoff
             e%gausv(2,l,g) = auxi(6)   !equivalent plastic strain
           ELSE
             CALL rubber1d(chi(1),lb(1),.FALSE.,mat%matdef(8),stres(1),aux)
           END IF
         ELSE
           lb(2) = SQRT(css(2))
           lb(3) = 1d0/(lb(1)*lb(2))    !unnecessary
           IF( plast )THEN
             lc(1) = lb(1)/EXP(e%gausv(1,l,g))  !elastic stretchs
             lc(2) = lb(2)/EXP(e%gausv(2,l,g))
             aux = e%gausv(3,l,g)       ! equivalent plastic strains
             CALL core06ps(chi,propi(1:2),mat%matdef(8),pflag,lc,stres,aux,ierr)
             IF (ierr == 1) CYCLE ELEMS       ! if there is an error CYCLE TO END (RETURN or EXIT are not allowed in parallel)
             IF( pflag )THEN  ! if plastic flow in the step
               e%gausv(1,l,g) = LOG(lb(1)/lc(1))   !plastic  log strains
               e%gausv(2,l,g) = LOG(lb(2)/lc(2))   !plastic  log strains
               e%gausv(3,l,g) = aux                !effective plastic strain
             END IF
             stres(1) = stres(1)/css(1)
             stres(2) = stres(2)/css(2)
           ELSE
             CALL rubberps(chi,lb,mat%matdef(8),stres)
           END IF
         END IF
       END SELECT

       IF( natst )THEN
         ! Computes 2nd P-K stress from Kirchhoff stresses
         stres(1) = stres(1)/css(1)
         IF( ns == 2 )stres(2) = stres(2)/css(2)
       END IF

      !***   compute Int(B**t*sigma) on element level

       force(1:ns) = force(1:ns)+stres(1:ns)*wei(l)     !t.t.t. integrated forces
       IF(shell)momen(1:ns) = momen(1:ns)+stres(1:ns)*zk*wei(l)  !integrated moments
     END DO
     IF( ntype == 1 )THEN
       force = force*thick*width            !Original thick because it's TLF
       momen = momen*thick*width            !Original thick because it's TLF
     ELSE
       force = force*thick                  !Original thick because it's TLF
       momen = momen*thick                  !Original thick because it's TLF
     END IF
     stint(pp(1,g):pp(2,g),ielem) = force(1:ns)             !store N for post-processing
     stint(pp(3,g):pp(4,g),ielem) = momen(1:ns)       !store M for post-processing
     dvolu = e%l0(2)/ngaus                        !Original length used (TLF)
     IF( ntype == 3 ) dvolu = dvolu*ro(1,g)*twopi
     force(1:ns) = force(1:ns)*dvolu
     momen(1:ns) = momen(1:ns)*dvolu

     !       actual configuration B Matrices
     ! compute B-matrix
     CALL bmem11(ntype,e%l0(2),ro(1,g),e%lambd,t,bmem,sh(:,g))
     CALL proma3(rl(1,2),bmem(1,1),force(1),4,ns,.FALSE.)  !membrane contributions
     ! compute standard bending B matrix

     IF(shell)THEN
       CALL bfle11(ntype,fixit,e%l0,ro(1,g),e%lambd,t,n,x,e%cab,sh(:,g),bfle,g)
       CALL proma3(rl(1,1),bfle(1,1,1),momen(1),8,ns,.FALSE.) !bending contributions
       ! add  bending equivalent forces for nodes on a side
       IF( siden(1) )THEN !first node on a side
         y0 = e%cab(2,1)*momen(1)*x(:,1)*sh(1,g)
         nm = -iffix(ln(2))  !position in array
         l  = nndpd(2,nm)
         !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
         resid(1:2,l) = resid(1:2,l) + y0
         !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
         l  = nndpd(3,nm)
         !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
         resid(1:2,l) = resid(1:2,l) - y0
         !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
       END IF
       IF( siden(2) )THEN !second node on a side
         y0 = e%cab(2,2)*momen(1)*x(:,4)*sh(2,g)
         nm = -iffix(ln(3))  !position in array
         l  = nndpd(2,nm)
         !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
         resid(1:2,l) = resid(1:2,l)  - y0
         !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
         l  = nndpd(3,nm)
         !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
         resid(1:2,l) = resid(1:2,l)  + y0
         !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
       END IF
     END IF

   END DO  !gauss point loop

   IF( ngaus == 1 .AND. cs > 0d0 .AND. simpl )THEN
     !one gauss point + stabilization for non-branching elements only
     aux = cs*curv(1,2)*dvolu   !stabilization moment
     CALL bfle11(ntype,fixit,e%l0,ro(1,g),e%lambd,t,n,x,e%cab,sh(:,g),bfle,0)
     CALL proma3(rl(1,1),bfle(1,1,1),aux,8,1,.FALSE.) !bending contributions
   END IF

   IF( siden(1) ) e%cab(1,1) = a0     !restore initial distortion
   IF( siden(2) ) e%cab(1,2) = a0     !restore initial distortion

   DO i=1,nnode
     l = e%lnods(i)
     IF(l > 0)THEN
       !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
       resid(1:2,l) = resid(1:2,l) + rl(1:2,i)
       !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
     END IF
   END DO
   ! include average rotation of the branching node
   IF( ASSOCIATED(ni) )THEN
     aux = 0d0                                         !initializes
     DO g=1,ngaus                                      !for each Gauss pont
       dvolu = sh(1,g)/ngaus                           !factor
       IF( ntype == 3 ) dvolu = dvolu*ro(1,g)*twopi    !axilsymmetric problems
       aux = aux + stint(pp(3,g),ielem) * dvolu        !Mss * dvolu
     END DO
     aux = 2d0*aux*e%lambd(1)                          !2*ls/l0* Int(Mss) *l0
     DO i=0,ni%nn                                      !for each node
       l = ni%lnods(i)                                 !global node
       !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
       resid(1:2,l) = resid(1:2,l) + aux*ni%bb(1:2,i)  !sum
       !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
     END DO
   END IF
   IF ( ASSOCIATED(nd) )THEN
     aux = 0d0                                         !initializes
     DO g=1,ngaus                                      !for each Gauss pont
       dvolu = sh(2,g)/ngaus                           !factor
       IF( ntype == 3 ) dvolu = dvolu*ro(1,g)*twopi    !axilsymmetric problems
       aux = aux + stint(pp(3,g),ielem) *dvolu         !Mss * dvolu
     END DO
     aux = 2d0*aux*e%lambd(1)                          !2*ls/l0* Int(Mss)*l0
     DO i=0,nd%nn                                      !for each node
       l = nd%lnods(i)                                 !global node
       !$ CALL omp_set_lock(lock1(MOD(l,MAXLOCK)))    !Start blocking acess in shared variables
       resid(1:2,l) = resid(1:2,l) - aux*nd%bb(1:2,i)  !sum
       !$ CALL omp_unset_lock(lock1(MOD(l,MAXLOCK)))  !End blocking acess in shared variables
     END DO
   END IF

   e => e%next
 END DO ELEMS
 !$OMP END PARALLEL
 !$ iel0 = iel0 + nelem
 !$ IF( iloop < nloop )THEN
 !$   DO i=1,nelem
 !$     heade => heade%next
 !$   END DO
 !$   IF( iloop + 1 == nloop ) THEN
 !$     IF( irest < maxlock)THEN
 !$       nelem = nelem + irest
 !$     ELSE
 !$       nelem = irest
 !$     END IF
 !$   END IF
 !$ END IF
 !$ END DO
 RETURN
 END SUBROUTINE resv11
