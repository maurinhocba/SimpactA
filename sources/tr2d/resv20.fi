 SUBROUTINE resv20( nele0, head0, ntype, coora, resid, ierr, ttime,  &
                    angdf, coorc, eulrf)
 !------------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT T-L 2-D TRIANGLE
 !
 !------------------------------------------------------------------------
 !$ USE omp_lib
 IMPLICIT NONE

        ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nele0,     & !total elements in set
                                 ntype        !problem type
 INTEGER (kind=4), INTENT(OUT) :: ierr        !flag for error detection
 REAL(kind=8),INTENT (IN) :: coorc(:,:),    & !previous coordinates
                             coora(:,:),    & !present  coordinates
                             ttime,         & !present total time
                             angdf            !default angle X_1 - Ort_1
 LOGICAL,INTENT (IN) :: eulrf                 !user finger tensor
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
 TYPE (ele20), POINTER  :: head0              !pointer to first element

 ! local variables
 LOGICAL :: newmt,   &! TRUE if material constant computation necessary
            elast,   &! TRUE if material is strictly elastic
            plast,   &! TRUE if plasticity is to be considered
            visco,   &! TRUE if viscoplasticity is to be considered
            mises,   &! TRUE if von Mises plasticity with multiplicative dec.
            sides(3)  ! TRUE if side element exist
 INTEGER (kind=4) :: osec,  & !associated material of previous element
                     numpt, & !number of points in curve
                     mtype, & !associated material type
                     is       !isotropic hardening model
 REAL (kind=8) :: fr(7),      & !Flow rule coefficients
                  km,gm,      & !bulk modulus and shear modulus
                  thick,      & !thickness for plane stress
                  propi(7),   & !Material properties
                  propv(3),   & !Viscoplastic material properties
                  chi(12),chib(4),    & !Hill-48 coeff
                  deatht        !end time for plasticity

 INTEGER (kind=4) ielem,nelem, & !element number
                  isec,        & !associated material
                  i,n            !different indexes

 REAL (kind=8) stres(4),   & !stresses (different measures)
               sigma(4),   & !stresses (different measures)
               stran(4),   & !C=U^2  also Log strains
               dstpl(4),   & !plastic strain increments
               efpst,      & !effective plastic strain
               strpl(4),   & !plastic strain (previous step)
               bmat(12,4), & !B matrix
               x(2,6),     & !nodal coordinates of the element
               xa(2,6),    & !auxiliar nodal coordinates of the element
               t(2,2),     & !deformation gradient at sides
               r1,r2,      & !eigevector components in local system
               lb(3),      & !eigenvalues
               aux,        & !auxiliar value
               c(4),       & !Elastic constitutive matrix for plane stress
               delta,      & !volumetric strain/3 and Jacobian
               press,      & !first stress invariant (mean pressure)
               rl(2,6)       !nodal internal equivalent forces

 REAL (kind=8), POINTER :: val(:,:)

 LOGICAL :: pflag
 TYPE (ele20), POINTER :: e, & !pointer to an element data
                          heade
 TYPE (section), POINTER :: sec    !pointer to section data
 TYPE (mater), POINTER :: mat  !pointer to an material data

 ! variables for error capture in strains evaluations
 REAL (kind=8), ALLOCATABLE :: conec(:),  & ! element nodes (conectivities)
                               estra(:),  & ! C=U^2  also Log strains
                               coorx(:,:)   ! nodal coordinates of the element

 INTERFACE
   INCLUDE 'stre14.h'
   INCLUDE 'ud_shell_2.h'
   INCLUDE 'ud_solid_3.h'
   INCLUDE 'corr17.h'
   INCLUDE 'corr20.h'
   INCLUDE 'rubberps.h'
   INCLUDE 'rubber2d.h'
 END INTERFACE

 !Auxiliar variables for parallel code
 !!$ INTEGER(kind=4), PARAMETER :: maxel=34816,minel=32768  !8 threads
 !!$ INTEGER(kind=4), PARAMETER :: maxel=18432,minel=16384  !4 threads
 !$ INTEGER(kind=4), PARAMETER :: maxel=10240,minel= 8192    !2 threads
 !$ INTEGER(kind=4):: nloop,iloop,irest
 INTEGER(kind=4):: piter,   & !
                   lim_inf, & ! first element for each processor
                   lim_sup, & ! last element for each processor
                   num_chk    ! number of elements chunks

 !$ INTEGER(kind=4),PARAMETER:: MAXLOCK=2048              !number of keys... this must be checked WBC
 !$ INTEGER(kind=4),SAVE     :: lock1(0:MAXLOCK-1)        !keys
 !$ LOGICAL,SAVE             :: inilock=.FALSE.           !initialization flag

 heade => head0
 nelem = nele0

 !Initialize exclusion variables for parallel code (only first time)
 !$ IF (.NOT.inilock) THEN                   !first time only
 !$   DO piter=0,MAXLOCK-1                     !for each key
 !$     CALL omp_init_lock(lock1(piter))          !unlock
 !$   END DO
 !$   inilock = .TRUE.                         !change flag
 !$ END IF
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

 !initializes parallel variables (one-processor/one-thread compatibility)
 lim_inf =  1    !first element
 lim_sup = nelem !last element

 !$OMP  PARALLEL                                                                           &
 !$OMP& SHARED  (num_chk,lock1,nelem,ntype,ierr,coorc,coora,ttime,angdf,eulrf,resid, &
 !$OMP&          heade,psecs,conec,estra,coorx)                                            &
 !$OMP& PRIVATE (piter,lim_inf,lim_sup,osec,e,newmt,elast,plast,visco,mises,sides,numpt,   &
 !$OMP&          mtype,is,fr,km,gm,thick,propi,propv,chi,deatht,ielem,isec,i,n,stran,dstpl,&
 !$OMP&          efpst,strpl,x,xa,t,r1,r2,lb,aux,c,delta,press,rl,val,pflag,sec,mat,stres, &
 !$OMP&          sigma,bmat,chib)

 !$ IF( omp_get_num_threads() > 1 )THEN
   !$ num_chk = (nelem+omp_get_num_threads()-1)/omp_get_num_threads() !multi-processors
   !$ lim_inf = MIN(nelem,1+omp_get_thread_num()*num_chk)     !Find the first element in the list
   !$ lim_sup = MIN(nelem,num_chk*(1+omp_get_thread_num()))   !Find the last element in the list
 !$ END IF

 !initializes variables
 sigma(4)  = 0d0  !for plane stress
 stres(4)  = 0d0  !for plane stress
 bmat(:,4) = 0d0  !for plane stress or strain

 e => heade                          !point to first element in the list
 DO piter=1,lim_inf-1                !Loop over the list to find the first element for EACH THREAD
   e => e%next                       !Point the first element of the list
 ENDDO

 !***  Computes nodal internal forces
 osec = -1        !no previous material

 ELEMS: DO piter=lim_inf,lim_sup     !loop over all the elements
   IF (ierr == 1) CYCLE ELEMS        ! if there is an error CYCLE TO END (RETURN or EXIT are not allowed in parallel)
   ielem = e%numel                   !element number
   isec = e%matno                   !associated material
   !     SET MATERIAL CONSTANTS
   newmt = osec /= isec             !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     sec => psecs(isec)%p
     osec = isec                    !keep present material
     mat => sec%mtbas                !point to associated material
     mtype = mat%mtype               !type of base material
     is    = mat%matdef(4)           !isotropic hardening model
     mises = .FALSE.                 !initializes
     ! ONLY MTYPES 1, 5, 6 and 30) ARE ALLOWED
     IF( mtype == 1)THEN              !Isotropic Elastic - Anisotropic Plastic
       IF( mat%matdef(3) /= 1 )THEN
         mises = mat%matdef(3) == 2     !if von Mises yield function
         elast = .FALSE.                !not elastic
         deatht = mat%propp(5)          !end of plasticity
         plast = deatht > ttime         !consider plasticity ?
         propi(1:7) = mat%propp(1:7)    !isotropic hardening parameters
         propi(5) = REAL(is,8)          !isotropic hardening model
         IF ( is == 5 )THEN             !points defined yield curve
           val => mat%chead%val
           numpt = mat%chead%np
         ELSE
           NULLIFY (val)
           numpt = 0
         END IF
         IF (mat%matdef(6) == 3) THEN !elastic-viscoplastic material
           visco = .TRUE.
           propv(1:3) = mat%propp(28:30) !viscoplastic parameters
         ELSE
           visco = .FALSE.
           propv = 0d0
         END IF
         chi = mat%propp(16:27)         !Hill yield function coefficients
         SELECT CASE (ntype)
         CASE ( 1 )  !plane stress
           SELECT CASE (mat%matdef(3))
           CASE (2,3) !Mises Hill48
             ! D matrix, derivative of yield function
             fr(1) = chi(2)+chi(3)  !g+h
             fr(2) = -chi(3)        !-h
             fr(3) = chi(1)+chi(3)  !f+h
             fr(4) = 2d0*chi(6)     !2n
             fr(5) = 2d0*chi(5)     !2m
             fr(6) = 2d0*chi(4)     !2l
             ! B matrix, flow rule matrix
             chib(1) = chi(8)+chi(9)  !g+h
             chib(2) = -chi(9)        !-h
             chib(3) = chi(7)+chi(9)  !f+h
             chib(4) = 2d0*chi(12)    !2n
           CASE (4)  !Hill-79
             fr(1:3) = chi(1:3)
           CASE (5)  !Hill-90
             fr(1:4) = chi(1:4)
             chib(1:4) = chi(5:8)
           END SELECT
         CASE ( 2 )  !plane strain
         ! A matrix, derivative of yield function (associative only)
           fr(1) =  chi(8) + chi(9) + 0.5d0  !0.5d0 added to make it diagonal
           fr(2) = -chi(8)          + 0.5d0  ! in case of von Mises
           fr(3) =  chi(7) + chi(8) + 0.5d0  ! it is not essential
           fr(4) = 2d0*chi(11)
           fr(5) = -chi(9)          + 0.5d0
           fr(6) = -chi(7)          + 0.5d0
           fr(7) =  chi(7) + chi(9) + 0.5d0
         CASE ( 3 ) !axilsymmetric
           aux = (chi(7) + chi(8) )/2d0      ! (F+G)/2
           fr(1) =  aux    + chi(9) + 0.5d0  !0.5d0 added to make it diagonal
           fr(2) = -aux             + 0.5d0  ! in case of von Mises
           fr(3) =  aux    + aux    + 0.5d0  ! it is not essential
           fr(4) =  chi(10) + chi(11)        ! L + M
           fr(5) = -chi(9)          + 0.5d0
           fr(6) = -aux             + 0.5d0
           fr(7) =  aux    + chi(9) + 0.5d0
         END SELECT
       ELSE
         elast = .TRUE.
         plast = .FALSE.
       END IF
       !  elasticity constants
       IF( ntype == 1 )THEN           !plane stress
         c(1:4) = mat%prope(7:10)   ! plane stress elasticity matrix
         thick  = sec%rprop(1)      ! thickness for plane stress
       ELSE                           ! plane strain or axilsymmetric
         km = mat%prope(4)*3d0      ! bulk modulus * 3
         gm = mat%prope(3)          ! shear modulus
       END IF
     ELSE IF( mtype == 5 )THEN                !orthotropic elastic
       IF( ntype == 1 ) THEN           !plane stress
         c(1) = mat%prope(16)   !(1,1)  11-11
         c(2) = mat%prope(17)   !(2,1)  11-22
         c(3) = mat%prope(18)   !(2,2)  22-22
         c(4) = mat%prope(19)   !(3,3)  12-12
       ELSE
         propi(1) = mat%prope(29)   !(1,1)  11-11
         propi(2) = mat%prope(30)   !(2,1)  11-22
         propi(3) = mat%prope(32)   !(4,1)  11-33
         propi(4) = mat%prope(34)   !(2,2)  22-22
         propi(5) = mat%prope(36)   !(4,2)  22-33
         propi(6) = mat%prope(39)   !(3,3)  12-12
         propi(7) = mat%prope(44)   !(4,4)  33-33
       END IF
       elast = .TRUE.
       plast = .FALSE.
     ELSE IF( mtype == 6 )THEN            ! Elastic Rubber
       !km = mat%prope(4)*3d0               ! bulk modulus * 3
       !gm = mat%prope(3)                   ! initial shear modulus
       propi(1:2) = mat%propp(1:2)         ! plastic properties
       plast = propi(2) > 0 ! .AND. props(5) > ttime
     ELSE IF( mtype == 30 )THEN            ! User defined
       plast = .TRUE.                    !consider plasticity ?
     END IF
     newmt = .FALSE.
   END IF

   IF(eulrf)THEN ! <--- spatial formulation
     !*** computes configuration in previous and present steps
     x(:,1:3) = coorc(:,e%lnods(1:3))  !previous coordinates of vertex nodes
     xa(:,1:3) = coora(:,e%lnods(1:3))  !present coordinates of vertex nodes
     DO i =4,6                          !for the rest of nodes in the patch
       n = e%lnods(i)                   !node number
       IF(n /= 0) THEN                  !if node exist
         x(:,i) = coorc(:,n)        !get original coordinates
         xa(:,i) = coora(:,n)       !get present coordinates
         sides(i-3) = .TRUE.            !side element exists
       ELSE
         sides(i-3) = .FALSE.       !side element does not exists
         x(:,i) = 0d0                !nothing
         xa(:,i) = 0d0                  !nothing
       END IF
     END DO
     !*** computes actual derivates from shape functions and present configuration
     CALL cdac20(e%numel,e%lnods,xa,e%cd,e%area1,ierr,thick,angdf,e%angle,ntype)

   ELSE         ! <--- lagrangian formulation
     !*** computes configuration in present step
     x(:,1:3) = coora(:,e%lnods(1:3))  !present coordinates of vertex nodes
     DO i =4,6                       !for the rest of nodes in the patch
       n = e%lnods(i)                !node number
       IF(n /= 0)THEN               !if node exist
         x(:,i) = coora(:,n)         !get present coordinates
         sides(i-3) = .TRUE.           !side element exists
       ELSE
         sides(i-3) = .FALSE.          !side element does not exists
         x(:,i) = 0d0                !nothing
       END IF
     END DO
   END IF

   !computes deformation gradient (total or updated)
   t(:,:) = MATMUL(x,e%cd(:,1:2))  ! in-plane derivatives
   IF( ntype == 2 ) lb(3) = 1d0                            ! hoop stretch
   IF( ntype == 3 ) lb(3) = DOT_PRODUCT(x(1,:),e%cd(:,3))  ! hoop stretch

   IF( mises .AND. mtype == 1 .AND. ntype /= 1)THEN
     CALL corr20(eulrf,t,lb(3),e%gausv,is,propi,propv,gm,km, &
                    e%stint,sigma,numpt,val,ierr)

   ELSE
     CALL stra20(t,ntype,lb(3),stran)  !metric tensor at center
     !computes Hencky (logarithmic) in-plane strains
     CALL lgst14(stran,r1,r2,lb,'RESV20',ierr)
     IF( ierr == 1 )THEN
       IF( ALLOCATED(conec) ) CYCLE ELEMS ! if was allocated by previous thread then cycle
       ALLOCATE(conec(6),estra(4),coorx(2,6))
       conec = e%lnods(1:6) !save data for error capture
       estra = stran        !
       coorx = x            !
       CYCLE ELEMS
     END IF

     IF( ntype == 1 )THEN            !plane stress
       IF( mtype == 1 )THEN          !one component material
         stres(1:2) = 0d0
         CALL stre14(stran,stres,c,propi,chib,fr,e%gausv,ierr,4,plast,elast,&
                     val,numpt,aux,mat%matdef(3),pflag)
       ELSE IF ( mtype == 5 )THEN ! Orthotropic elastic
         stres(1) = c(1)*stran(1)+c(2)*stran(2)
         stres(2) = c(2)*stran(1)+c(3)*stran(2)
         stres(3) = c(4)*stran(3)
       ELSE IF ( mtype == 6 )THEN ! Hyperelastic
         !IF( plast )THEN
           !!       elastic-plastic behavior (NEW)
           !efpst = e%gausv(5)                !effect plastic strain
           !IF( efpst > 0 )THEN
           !  !elastic (trial) strains
           !  strpl(1:4) = e%gausv(1:4)            !previous (twice) plastic strains
           !  dstpl(1:4) = stran(1:4) - strpl(1:4)   !trial Elastic log strains (Twice)
           !  ! compute eigen-decomposition of elastic part
           !  CALL eige17(dstpl(1),s1,s2,lc(1))
           !  lc(3) = dstpl(4)
           !  lc = EXP(lc/2d0)  !lc are twice the principal shear strains
           !ELSE
           !  lc = lb
           !  s1 = r1
           !  s2 = r2
           !END IF
           !!
           !!   compute deviatoric principal stresses and check plasticity
           !!
           !CALL ogdene(bet(1),lc(1),eprop(1,1),nto,.FALSE.) !although .TRUE.
           !IF( plast ) THEN
           !  ss(1) = s1
           !  ss(2) = s2
           !  CALL core10 (eprop(1,1), propi(1), newmt, pflag, lc(1), bet(1), &
           !               ss ,efpst, dstpl(1), ielem, ierr, 4, nto)
           !  IF( pflag )THEN  ! if plastic flow in the step
           !    e%gausv(5) = efpst      !effective plastic strain
           !    e%gausv(1:4) = stran(1:4) - dstpl(1:4)   !plastic  log strains
           !  END IF
           !END IF
           !
           !! compute (Total) stress tensor
           !stres(1) = s1*bet(1)*s1 + s2*bet(2)*s2 - bet(3)  !(1,1)
           !stres(2) = s2*bet(1)*s2 + s1*bet(2)*s1 - bet(3)  !(2,2)
           !stres(3) = s1*bet(1)*s2 - s2*bet(2)*s1           !(1,2)
           !stres(4) = 0d0                                   !(3,3)
         !ELSE
           CALL rubberps(chi,lb,mat%matdef(8),sigma,r1=r1,r2=r2)
         !END IF
       ELSE IF( mtype == 30 )THEN     ! User defined material
         CALL ud_shell_2(stran,stres,e%gausv,ttime,ielem,ierr, &
                        newmt,mat%matdef(8),mat%props,mat%chead)

       END IF
     ELSE    !(ntype == 2 or 3) Plane Strain or Axisymmetric
       stran(4) = LOG(lb(3))
       IF( mtype == 1 )THEN          !one component material
         ! compute volumetric and deviatoric parts
         delta = ( stran(1) + stran(2) + stran(4) )/3d0
         stran(1) = 2d0*(stran(1) - delta)  !twice the deviatoric shear strain
         stran(2) = 2d0*(stran(2) - delta)
         stran(4) = 2d0*(stran(4) - delta)
         !elastic (trial) strains
         IF( .NOT.elast )THEN
           strpl(1:4) = e%gausv(1:4)               !previous (twice) plastic strains
           stran(1:4) = stran(1:4) - strpl(1:4)    !trial Elastic log strains
         END IF
         stres= gm*stran                         !Trial elastic shear stresses
         IF( plast ) THEN            !Consider plasticity
           efpst = e%gausv(5)        !effect plastic strain
           CALL corr17(stres(1),stres(2),stres(3),stres(4),efpst,gm, &
                       propi,fr,ierr,dstpl,pflag,is,ielem,numpt,val)
           IF(ierr == 1) CYCLE ELEMS             !no convergence in plasticity
           IF( pflag )THEN                   !if plastic flow
             e%gausv(1:4) = e%gausv(1:4) + dstpl(1:4)  !total plastic shear strains
             e%gausv(5)   = e%gausv(5)   + efpst       !Total effect. plastic strain
           END IF
         END IF
       ELSE IF ( mtype == 5 )THEN ! Orthotropic elastic
         stres(1) = propi(1)*stran(1)+propi(2)*stran(2)+propi(3)*stran(3)
         stres(2) = propi(2)*stran(1)+propi(4)*stran(2)+propi(5)*stran(3)
         stres(3) = propi(6)*stran(3)
         stres(4) = propi(3)*stran(1)+propi(5)*stran(2)+propi(7)*stran(3)
       ELSE IF ( mtype == 6 )THEN ! Hyperelastic
         !IF(plast)THEN
           !efpst = e%gausv(5)                !effect plastic strain
           !IF( efpst > 0 )THEN
           !  !elastic (trial) strains
           !  strpl(1:4) = e%gausv(1:4)            !previous (twice) plastic strains
           !  dstpl(1:4) = stran(1:4) - strpl(1:4)   !trial Elastic log strains (Twice)
           !  ! compute eigen-decomposition of elastic part
           !  CALL eige17(dstpl(1),s1,s2,lc(1))
           !  lc(3) = dstpl(4)
           !  lc = EXP(lc/2d0)  !lc are twice the principal shear strains
           !ELSE
           !  aux = EXP(delta)  !aux = jacobian (delta = ln(J)/3
           !  lc = lb/aux
           !  s1 = r1
           !  s2 = r2
           !END IF
           !!
           !!   compute deviatoric principal stresses and check plasticity
           !!
           !CALL ogdene(bet(1),lc(1),eprop(1,1),nto,.FALSE.)
           !IF( plast ) THEN
           !  ss(1) = s1
           !  ss(2) = s2
           !  CALL core10 (eprop(1,1), propi(1), newmt, pflag, lc(1), bet(1), &
           !               ss ,efpst, dstpl(1), ielem, ierr, 4, nto)
           !  IF( pflag )THEN  ! if plastic flow in the step
           !    e%gausv(5) = efpst      !effective plastic strain
           !    e%gausv(1:4) = stran(1:4) - dstpl(1:4)   !plastic  log strains
           !  END IF
           !END IF
           !
           !! compute (deviatoric) stress tensor
           !stres(1) = s1*bet(1)*s1 + s2*bet(2)*s2    !(1,1)
           !stres(2) = s2*bet(1)*s2 + s1*bet(2)*s1    !(2,2)
           !stres(4) = bet(3)                         !(3,3)
           !stres(3) = s1*bet(1)*s2 - s2*bet(2)*s1    !(1,2)
         !ELSE
           CALL rubber2d(chi,lb,mat%matdef(8),sigma,r1,r2)
         !END IF
       ELSE IF( mtype == 30 )THEN     ! User defined material
         CALL ud_solid_3(stran,delta,stres,e%gausv,ttime,ielem,ierr, &
                         newmt,mat%matdef(8),mat%props,mat%chead)
       END IF
     END IF

     IF( ntype /= 1 )THEN  !only for one-component material
       IF( mtype == 1 )THEN
         press = km*delta
         stres(1) = stres(1) + press
         stres(2) = stres(2) + press
         stres(4) = stres(4) + press
       ELSE
         ! nothing yet
       END IF
     END IF

     IF( mtype == 6 )THEN
       e%stint = sigma  !2nd Piola Kirchhoff stresses
     ELSE
       ! Computes Hencky stress on the natural Frame
       sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2+2d0*stres(3)*r1*r2
       sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1-2d0*stres(3)*r1*r2
       sigma(3) =(stres(2)-stres(1))*r1*r2+stres(3)*(r1*r1-r2*r2)
       ! Computes 2nd P-K stress on the natural Frame
       stres(1) = sigma(1)/lb(1)**2
       stres(2) = sigma(2)/lb(2)**2
       IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
         stres(3) = sigma(3)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
       ELSE                                !lb(1) = lb(2)
         stres(3) = sigma(3)/lb(1)/lb(2)
       END IF
      ! Computes 2nd P-K on the Lagrangian Frame
       sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2-2d0*stres(3)*r1*r2
       sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1+2d0*stres(3)*r1*r2
       sigma(3) =(stres(1)-stres(2))*r1*r2+stres(3)*(r1*r1-r2*r2)
       IF( ntype /= 1 )sigma(4) = stres(4)/lb(3)**2
       ! compute Kirchhoff stress tensor for post-process
       e%stint(1) = sigma(1)*t(1,1)*t(1,1)+sigma(2)*t(1,2)*t(1,2)+2d0*sigma(3)*t(1,2)*t(1,1)
       e%stint(2) = sigma(1)*t(2,1)*t(2,1)+sigma(2)*t(2,2)*t(2,2)+2d0*sigma(3)*t(2,1)*t(2,2)
       e%stint(3) = sigma(1)*t(1,1)*t(2,1)+sigma(2)*t(1,2)*t(2,2)+sigma(3)*(t(1,1)*t(2,2)+t(1,2)*t(2,1))
       e%stint(4) = stres(4)
     END IF
   END IF
   !***   compute Int(B**t*sigma) on element level

   sigma = sigma*e%area1         !Original dvol used (TLF)

   !       actual configuration B Matrx
   CALL bmat20(e%cd,bmat,t,sides,ntype,lb(3))

   CALL proma3(rl(1,1),bmat(1,1),sigma(1),12,4,.TRUE.)  !nodal contributions

   DO i=1,6             !for each node in the element
     n = e%lnods(i)     !node
     IF( n > 0) THEN
       !$ CALL omp_set_lock(  lock1(MOD(n,MAXLOCK)) )  !Start blocking acess in shared variables
       resid(1:2,n) = resid(1:2,n)+rl(1:2,i) !sum internal force
       !$ CALL omp_unset_lock(lock1(MOD(n,MAXLOCK)) )  !End blocking acess in shared variables
     END IF
   END DO
   e => e%next          !point to next element in the list
 END DO ELEMS
 !$OMP END PARALLEL
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

 !print out the error data and stop the program
 IF( ALLOCATED(conec) )THEN
   WRITE(55,"(6i5)",ERR=9999) conec
   WRITE(55,"(4e12.4)",ERR=9999) estra
   WRITE(55,"(2e15.4)",ERR=9999) coorx
   DEALLOCATE(conec,estra,coorx)
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE resv20
