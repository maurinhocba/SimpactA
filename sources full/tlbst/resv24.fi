 SUBROUTINE resv24( nele0,head0, iffix, coora, resid, logst, ierr,       &
                    bottom, top, coorb, coort, ifact, nbs, bhead, stint, &
                    shearc,moments,factors,ninv,velnp)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT TL BBST++
 !
 !-------------------------------------------------------------------
 !$ USE omp_lib
 USE gvar_db, ONLY : updiv,newiv,elastic
 USE kinc_db, ONLY : ndepd
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN) :: logst       !flags for small or large strain analysis
 LOGICAL, INTENT(IN) :: bottom,top  !flags for top and bottom surface computation
 INTEGER (kind=4), INTENT(IN) :: nele0,     & !number of elements in the set
                                 iffix(:),  & !nodal fixities
                                 nbs          !number of branching sides
 INTEGER (kind=4), INTENT(OUT) :: ierr,     & !flag for error detection
                               ifact(:)       !weigths for surface computation
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             velnp(:,:)       !present velocities (not used yet)

 REAL(kind=8),INTENT (IN OUT) :: resid(:,:),& !internal nodal forces
                                 coorb(:,:),& !bottom surface coordinates
                                 coort(:,:)   !top surface coordinates
 TYPE (ele24), POINTER  :: head0              !pointer to first element
 TYPE (sideb), POINTER  :: bhead              !pointer to first element
 REAL(kind=8), POINTER :: stint(:,:)          !moments and shear
 INTEGER(kind=4), INTENT (IN) :: shearc       !/=0 if shear are to be computed
 REAL(kind=8), POINTER :: moments(:,:),factors(:)  !bending moments and factors for its computation
 INTEGER(kind=4), POINTER :: ninv(:)          !node inverse order

 ! local variables
 REAL (kind=8), PARAMETER :: fac(4) = (/ -1d0, 1d0, 1d0, -1d0 /)  !quad positions

 TYPE (ele24), POINTER :: e       !pointer to an element data
 TYPE (ele24), POINTER  :: heade  !pointer to first element
 TYPE (sideb), POINTER  :: si     !pointer to branching side

 LOGICAL :: isidf(3),& ! TRUE if side is clamped or symmetry line
            branch,  & ! TRUE if a branching side exists at the element
            mti        ! TRUE if more than one solid side

 INTEGER (kind=4) sides(3),               & ! 1 if side element exist
                  nbn,nbo,lnb(18),lna(4)    !nodal connectivities for branching element
 INTEGER (kind=4)  :: iel,i3,jj,nelem,iel0  !element order

 REAL (kind=8) stra1(6),   & !first and second fundamental forms
               bmat(3,6,6),& !membrane and bending B matrix
               x(3,6),     & !nodal coordinates of the patch
               v(3,6),     & !nodal velocities of the patch
               t(3,15),    & !Local triad + local in-plane vector at side elms
               rl(3,6),    & !nodal internal equivalent forces
               lh,sh(4),tt(3),k0

 REAL (kind=8), ALLOCATABLE :: bb(:,:,:),rb(:,:)                     !branching sides

 LOGICAL :: large                 ! TRUE if log strains must be used
 LOGICAL, ALLOCATABLE :: pflag(:) ! TRUE if plastic flow in the step

 LOGICAL :: newmt,   &! TRUE if material constant computation necessary
            natst,   &! TRUE for large strain analysis
            elast,   &! TRUE if material is strictly elastic
            shell,   &! TRUE if bending included
            tplast,  &! TRUE if plastic properties are temperature dependant
            plast,   &! TRUE if plasticity is to be considered
            shear     ! TRUE if transverse shear stresses are considered in the constitutive relation


 INTEGER (kind=4) ielem,     & !element number
                  i,j,k,n,l, & !different indexes
                  error        !error flag

 INTEGER (kind=4) isec,  & !associated material
                  nlayr, & !number of layers
                  numpt, & !number of points in curve
                  nyoung,nnu,ncons0,nuniax,nconsn,nalpha, & !number of points in curve
                  osec,  & !associated material of previous element
                  nvar     !number of internal variables per layer

 REAL (kind=8) stres(3),  & !stresses (different measures)
               sigma(3),  & !stresses (different measures)
               stran(3),  & !C=U^2  also Log strains
               stine(10), & !t.t.t integrated stresses (forces and moments)
               r1,r2,     & !eigevenctor components in local system
               lb(3),     & !eigenvalues, also used as aux for temperature comp.
               lc(3),     & !point loads leading to torsor moments
               thnew,     & !present thickness
               zk,        & !distance to mid surface
               pw(3),     & !integrated plastic work of the step
               young,nu,  & !elastic modulus and poisson ratio
               vm,        & !auxiliar value (Von Mises stress)
               aux,area1, & !auxiliar value and element area
               t0,t1,t2,j0,  & !thermical dilatation coeff
               thick,     & !thickness (original)
               c(4),      & !Elastic constitutive matrix for plane stress
               alpha,     & !thermical dilatation coeff
               propi(13), & !Plastic Material properties
               chi(12),   & !Hill 48 coefficients
               chib(4),chid(6)  !yield and potential coefficients
 ! Gauss points throught the thickness
 REAL (kind=8), ALLOCATABLE :: thf(:),wei(:),shf(:)

 REAL (kind=8), POINTER :: val(:,:)
 REAL (kind=8), ALLOCATABLE :: varin(:,:)                      !internal variables

 TYPE (section), POINTER :: sec         !pointer to a section data
 TYPE (mater), POINTER :: mat           !pointer to a material data
 TYPE (curve), POINTER :: cyoung,cnu,ccons0,cuniax,cconsn,calpha  !pointers to different curves

 INTERFACE
   INCLUDE 'stre14.h'
 END INTERFACE


 !Auxiliar variables for parallel code
 !!$ INTEGER(kind=4), PARAMETER :: maxel=34816,minel=32768  !8 threads
 !!$ INTEGER(kind=4), PARAMETER :: maxel=18432,minel=16384  !4 threads
 !$ INTEGER(kind=4), PARAMETER :: maxel=10240,minel= 8192    !2 threads
 !$ INTEGER(kind=4):: nloop,iloop,irest    !
 INTEGER(kind=4):: lim_inf,  & ! first element for each processor
                   lim_sup,  & ! last element for each processor
                   num_chk     ! number of elements chunks

 !$ INTEGER(kind=4),PARAMETER:: MAXLOCK=2048              !number of keys
 !$ INTEGER(kind=4),SAVE     :: lock1(0:MAXLOCK-1)        !keys
 !$ LOGICAL,SAVE             :: inilock=.FALSE.           !initialization flag

 NULLIFY (val)                     !no point-defined isotropic hardening (for stre14 call)
 heade => head0
 nelem = nele0
 iel0 = 0

 !***  compute data associated to branching sides  only if branching sides exists
 IF( nbs > 0 )   CALL bran24(nbs,bhead,coora)         !compute angles and B matrices

 !Initialize exclusion variables for parallel code (only first time)
 !$ IF (.NOT.inilock) THEN                   !first time only
 !$   DO i=0,MAXLOCK-1                     !for each key
 !$     CALL omp_init_lock(lock1(i))          !unlock
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

 !initialize variables
 num_chk =  1    !only one processor or thread
 lim_inf =  1    !first element
 lim_sup = nelem !last element

 !$OMP  PARALLEL                                                                           &
 !$OMP& SHARED (num_chk, lock1, psecs, pmats, ndepd,                                       &
 !$OMP&         updiv, logst, bottom, top, coort, coorb, iffix, coora, velnp,              &
 !$OMP&         shearc, ninv, factors,  val,                                               &
 !$OMP&         heade, ierr, newiv, resid, stint,                                          &
 !$OMP&         nbs, iel0                 )                                                &
 !$OMP& PRIVATE(lim_inf, lim_sup, e, si, isidf, branch, mti, sides, nbn, nbo, lnb, lna,    &
 !$OMP&         iel, i3, jj, stra1, bmat, x, v, t, rl, lh, sh, tt, k0,                     &
 !$OMP&         bb, rb,                                                                    &
 !$OMP&         large, pflag,                                                              &
 !$OMP&         newmt,  natst, elast, shell, tplast, plast, shear,                   &
 !$OMP&         ielem, i, j, k, n, l, error,                                               &
 !$OMP&         isec, nlayr, numpt, osec, nvar,                                            &
 !$OMP&         nyoung, nnu, ncons0, nuniax, nconsn, nalpha,                               &
 !$OMP&         stres, sigma, stran, stine, r1,r2, lb, lc, thnew, zk,                      &
 !$OMP&         aux, area1, pw, young, nu, vm, t0, t1, t2, j0,                             &
 !$OMP&         thick, c, alpha, propi, chi, chib, chid,                                   &
 !$OMP&         thf, wei, shf, varin,                                                      &
 !$OMP&         cyoung, cnu, ccons0, cuniax, cconsn, calpha,                               &
 !$OMP&         sec, mat)
 !$ IF( omp_get_num_threads() > 1 )THEN
   !$ num_chk = (nelem+omp_get_num_threads()-1)/omp_get_num_threads()
   !$ lim_inf = MIN(nelem,1+omp_get_thread_num()*num_chk) + iel0    !Find the first element in the list
   !$ lim_sup = MIN(nelem,num_chk*(1+omp_get_thread_num())) + iel0  !Find the last element in the list
 !$ END IF

 IF( nbs > 0 ) THEN                     !only if branching sides exists
   nbo = 6                              !initializes
   ALLOCATE( bb(3,nbo,6),rb(3,nbo) )    !get memory for auxiliar matrices
 END IF

 e => heade                          !point to head of the list
 DO iel=iel0+1,lim_inf-1             !Loop over the list to find the first element in the list
   e => e%next                       !Point the first element of the list
 END DO

 !***  Computes nodal internal forces, integrating trought the thickness

 osec = -1                           !no previous material
 error = 0

 ELEMS: DO iel=lim_inf,lim_sup    !loop over the elements ()
   IF (ierr == 1) CYCLE ELEMS       ! if there is an error CYCLE TO END (RETURN or EXIT are not allowed in parallel)
   ielem = e%numel                   !element number
   isec = e%matno                    !associated material
   newmt = osec /= isec              !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties

     sec => psecs(isec)%p             !point to section
     nlayr = sec%iprop(1)             !number of layers
     nvar  = sec%iprop(2)             !number of internal variables per layer
     thick = sec%rprop(1)             !original thickness
     osec = isec                      !keep present section
     shear = sec%iprop(4) == 1        !consider transverse shear stresses in constitutive relation
     !visco = .FALSE.                  !presently viscous properties are not available for MAT 4
     !***** Allocate vectors thf, wei and shf *****************************************
     IF (ALLOCATED(thf)) DEALLOCATE(thf, wei, shf, pflag)
     ALLOCATE(thf(nlayr), wei(nlayr), shf(nlayr), pflag(nlayr))
     !**************************SHELL2*************************************************
     mat => sec%mtbas                 !point to associated material
     shell = nlayr > 1                !membrane or shell
     IF( .NOT.shell )rl(:,4:6) = 0d0
     natst = logst      !use log strains
     elast = mat%matdef(3) == 1       !elastic
     CALL gaussq(nlayr,thf(1),wei(1)) !integration points through the thickness
     thf(1:nlayr) = thf(1:nlayr)/2d0  !positions
     wei(1:nlayr) = wei(1:nlayr)/2d0  !weights
     young = mat%prope(1)             !young modulus
     nu = mat%prope(2)                !poisson ratio
     alpha = mat%prope(6)/3d0         !Thermical Dilatation Coeff (3 because used with 3 nodes)
     nyoung= INT(sec%mtbas%props(1))  !associated curve for Young modulus
     nnu   = INT(sec%mtbas%props(2))  !associated curve for Young modulus
     nalpha= INT(sec%mtbas%props(6))  !associated curve for alpha
     IF( nyoung > 0 ) CALL cur_point( nyoung,sec%mtbas%chead,cyoung) !point to the curve
     IF( nnu    > 0 ) CALL cur_point( nnu   ,sec%mtbas%chead,cnu   ) !point to the curve
     IF( nalpha > 0 ) CALL cur_point( nalpha,sec%mtbas%chead,calpha) !point to the curve

     IF( elast )THEN
       plast = .FALSE.
     ELSE   !for elasto-plastic mats
       propi(5) = REAL( mat%matdef(4),8) ! isotropic hardening model 0-1-2 only
       propi(1:3) = mat%propp(1:3)       ! isotropic hardening parameters
       IF( propi(5) == 3 ) propi(2) = (propi(2)/propi(1))**(1d0/propi(3)) !bar epsilon IHLUDW
       chi      = mat%propp(16:27)       ! hill coefficients
       plast = propi(1) > 0  !consider plasticity ?
       IF(plast)THEN                     !plasticity present
         SELECT CASE (mat%matdef(3))     !according to yield surface
         CASE (2,3) !Mises Hill48
           ! D matrix, derivative of yield function
           chid(1) = chi(2)+chi(3)  !g+h
           chid(2) = -chi(3)        !-h
           chid(3) = chi(1)+chi(3)  !f+h
           chid(4) = 2d0*chi(6)     !2n
           chid(5) = 2d0*chi(5)     !2m
           chid(6) = 2d0*chi(4)     !2l
           ! B matrix, flow rule matrix
           chib(1) = chi(8)+chi(9)  !g+h
           chib(2) = -chi(9)        !-h
           chib(3) = chi(7)+chi(9)  !f+h
           chib(4) = 2d0*chi(12)    !2n
         CASE (4)  !Hill-79
           chid(1:3) = chi(1:3)
         CASE (5)  !Hill-90
           chid(1:4) = chi(1:4)
           chib(1:4) = chi(5:8)
         END SELECT
         ncons0= INT(sec%mtbas%props(3))        !associated curve for CONS0
         nuniax= INT(sec%mtbas%props(4))        !associated curve for UNIAX
         nconsn= INT(sec%mtbas%props(5))        !associated curve for CONSN or HARDS
         IF( ncons0 > 0 ) CALL cur_point( ncons0,sec%mtbas%chead,ccons0) !point to the curve
         IF( nuniax > 0 ) CALL cur_point( nuniax,sec%mtbas%chead,cuniax) !point to the curve
         IF( nconsn > 0 ) CALL cur_point( nconsn,sec%mtbas%chead,cconsn) !point to the curve
       ELSE
         ncons0= 0
         nuniax= 0
         nconsn= 0
       END IF
     END IF
     c(1:4) = mat%prope(10:13)          ! plane stress elasticity matrix (initial)
     numpt = nyoung + nnu + ncons0 + nuniax + nconsn + nalpha
     tplast = plast .AND. ((ncons0+nuniax+nconsn) > 0)

     !- Obtain shear factors from the rprop array ---------------------------------
     shf(1:nlayr) = sec%rprop(13:12+nlayr)  !shear factors (isotropic material)

     IF( plast .OR. logst )THEN
       IF(ALLOCATED( varin ) )THEN
          IF( SIZE(varin,1) /= nvar .OR. SIZE(varin,2) /= nlayr ) DEALLOCATE(varin)
       END IF
       IF( .NOT.ALLOCATED(varin)) ALLOCATE( varin(nvar,nlayr) )
       IF(plast) plast = .NOT.elastic       ! for elastic behavior with previous internal variables
     END IF
     newmt = .FALSE.                        !same material than previous ?
   END IF

   !*** computes actual configuration
   lb(1) = SUM(tempe(1,e%lnods(1:3)))                  !average temperature (center or bottom)
   IF(ndoft > 1) lb(2) = SUM(tempe(2,e%lnods(1:3)))    !average temperature (bottom or top)
   IF( lb(1) /= 0d0 ) THEN
     IF( numpt > 0 )THEN !temperature dependant properties
       ! temperature at element center (middle surface)
       IF( ndoft == 2)THEN
         aux = (lb(1)+lb(2))/6d0   !average temperature
       ELSE
         aux = lb(1)/3d0   !average temperature
       END IF
       i = 1
       IF( nyoung > 0 .OR. nnu > 0 ) THEN
         IF( nyoung > 0 ) young = inte_cr(cyoung%val,nyoung,aux,i)
         IF( nnu    > 0 ) nu    = inte_cr(cnu%val   ,nnu   ,aux,i)
         ! plane stress elasticity matrix
         vm = 1d0-nu*nu
         c(1) = young/vm
         c(2) = c(1)*nu
         c(3) = c(1)
         c(4) = young/(2d0*(1d0+nu))
       END IF
       IF( nalpha > 0 ) alpha = inte_cr(calpha%val,nalpha,aux,i)/3d0
       IF( tplast )THEN       !termo-plasticity
         SELECT CASE (INT(propi(5)))           !according to IH properties
         CASE (1)           !perfect plasticity
           IF( nuniax > 0 ) propi(1) = inte_cr(cuniax%val,nuniax,aux,i)
         CASE (2)           !linear hardening
           IF( nuniax > 0 ) propi(1) = inte_cr(cuniax%val,nuniax,aux,i)
           IF( nconsn > 0 ) propi(2) = inte_cr(cconsn%val,nconsn,aux,i)
         CASE (3)           !Ludwig-Nadai approximation
           IF( ncons0 > 0 ) propi(1) = inte_cr(ccons0%val,ncons0,aux,i)
           IF( nuniax > 0 ) propi(2) = inte_cr(cuniax%val,nuniax,aux,i)
           IF( nconsn > 0 ) propi(3) = inte_cr(cconsn%val,nconsn,aux,i)
           IF( nuniax > 0 )THEN
             aux = (propi(2)/propi(1))**(1d0/propi(3))
           ELSE
             aux = (mat%propp(2)/propi(1))**(1d0/propi(3))
           END IF
         END SELECT
       END IF
     END IF
   END IF

   branch = ANY(e%lnods < 0 )  !a branching exists for this element
   IF( branch )THEN
     nbn = 3                   !first three nodes of the main triangle
     lnb(1:3) = e%lnods(1:3)
   END IF
   isidf = .FALSE.                  !initializes
   x(:,1:3) = coora(:,e%lnods(1:3)) ! get coordinates of the main triangle
   !IF( visco)THEN
   !  v(:,1:3) = velnp(:,e%lnods(1:3)) ! get velocities of the main triangle
   !  v(:,4:6) = 0d0                 ! initializes extra node velocities
   !END IF
   DO i =1,3                        ! for each extra-node in the patch
     i3= i + 3                      ! position of extra node
     n = e%lnods(i3)                ! node number
     IF(n > 0) THEN                 ! adjacent element exist ?
       sides(i) = 1                   ! side element exist
       x(:,i3) = coora(:,n)           ! get coordinates of the side node
       !IF( visco) v(:,i3) = velnp(:,n) ! get velocities of the extra node
       IF(branch)THEN                 ! if a branch exist
         nbn = nbn + 1                  ! update number of nodes
         lnb(nbn) = n                   ! keep node
       END IF
     ELSE IF( n == 0 )THEN            ! boundary node
       j = e%lnods(hh(3,i))             ! first node of the side
       k = e%lnods(hh(2,i))             ! second node of the side
       IF( iffix(j) < 0 .AND. iffix(k) < 0 )THEN !node on a side?
         l = -iffix(j)  !compute position of the i+1 node
         lna(1:2) = nndpd(2:3,l)          ! add first two nodes of the solid quadrilateral
         l = -iffix(k)                    ! compute position of the i+2 node
         lna(3:4) = (/ nndpd(3,l), nndpd(2,l) /)  ! add third and fourth nodes of the solid quadrilateral
         t(:,1:4) =  coora(:,lna(:))      ! quadrilateral coordinates
         t(:,5) = -t(:,1)+t(:,2)+t(:,3)-t(:,4) !phi,xita derivative
         t(:,6) = x(:,hh(2,i))-x(:,hh(3,i))    !side vector
         CALL vecpro(t(1,5),t(1,6),x(1,i3))    !outward normal vector
         CALL vecuni(3,x(1,i3),aux)       !unit normal to the quadrilateral
         sides(i) = 2                     ! (2) ==> side element is solid
       ELSE
         x(1:3,i3) = 0d0                  ! null values (unnecessary perhaps)
         sides(i) = 0                     ! (0) ==> side does NOT exists
       END IF
       isidf(i) = ( iffix(j) /= 0 .AND. iffix(k) /= 0 )  !clamped side ?
     ELSE !IF( n == -1 )THEN ! branching side
       sides(i) = -1                    !(-1) => side element is a branching side
       si => e%si(i)%p                  !point SI to the corresponding branching side
       ! assign present distorsion at side
       j = -n                           !position
       IF( si%lnods(j) > 0 )THEN        !if order
         e%gamma(i) = si%gamma(j)         !(+)
       ELSE                               !reversed order
         e%gamma(i) = -si%gamma(j)        !(-)
       END IF
       ! generate connectivities for enlarged B matrix
       k = nbn + 1                      !next position
       nbn = nbn + si%nn-1              !update number of nodes
       lnb(k:nbn) = (/ ABS(si%lnods(1:j-1)), ABS(si%lnods(j+1:si%nn)) /) !+/-
     END IF
   END DO

   ! consider thermal strains
   k0 = 0d0
   IF( alpha > 0 )THEN              ! compute coefficients for J(temp) evaluation
     SELECT CASE (ndoft)
     CASE (1)                         !constant
       t0 = 1d0+alpha*lb(1)              ! T0
       t1 = 0d0
       t2 = 0d0
     CASE (2)                         !linear
       t0 = 1d0+alpha*( lb(1)+lb(2))/2d0 ! T0
       t1 =     alpha*(-lb(1)+lb(2))     ! T1
       t2 = 0d0
     CASE (3)                         !quadratic
       lb(3) = SUM(tempe(3,e%lnods(1:3)))
       t0 = 1d0+alpha*lb(1)                         ! T0
       t1 =     alpha*(          -lb(2)+lb(3))      ! T1
       t2 =     alpha*(-2d0*lb(1)+lb(2)+lb(3))*2d0  ! T2
     END SELECT
     k0 = t1/thick
   END IF
   ! array LB released after this to be used as eigenvalues
   ! compute first and second fundamental forms
   CALL stra24(e%a,e%b,e%c,e%cd,e%ci,e%a0,e%area1,t,e%lb,x,sides,isidf, &
               e%mems,.FALSE.,stra1,e%gamma,k0)
   e%stra1 = stra1                  ! keep Mid Surface fundamental forms
   IF( alpha > 0d0 ) e%lb = e%lb*t0**3     ! modify thickness ratio including thermal expansion

   thnew = thick*e%lb                  !real thickness at new configuration
   IF( alpha > 0 )THEN                 !divide 2nd and 3rd coefficient by the thickness
     t1 = t1/thnew
     t2 = t2/thnew/thnew
   END IF

   ! compute external contact surfaces
   IF( bottom .OR. top) THEN      !bottom or top surfaces required
     ! determine element normal to compute top-bottom surfaces
     IF( ALL(.NOT.isidf))THEN     !no clamped sides
       tt = t(:,3)                !element normal
     ELSE                         !one or more clamped sides
       n = 0                      !initializes number of clamped sides
       DO i=1,3                   !loop over each side
         IF( isidf(i) )THEN       !if side clamped
           n = n+1                !increase number of clamped sides
           rl(:,n) = e%c(:,i,2)/2d0  !normal to plane of symmetry
         END IF
       END DO
       IF( n == 1)THEN            !for one clamped side (or symmetry plane)
         CALL vecuni(3,rl(1,1),aux)
         aux = DOT_PRODUCT(rl(:,1),t(:,3))    !proyect normal over plane
         rl(:,3) = t(:,3) - aux*rl(:,1)       !proyect normal over plane
       ELSE                       !mfor two clamped side (or symmetry planes)
         CALL vecpro(rl(1,1),rl(1,2),tt) !normal to both planes
       END IF
       CALL vecuni(3,tt,aux)             !unit vector
     END IF
     tt = tt*thnew          !
   END IF
   area1 = e%area1       !element original area
   stine = 0d0           !Integrated forces, moments, shears and equivalent stresses
   ! compute shear forces based on Moments at previous step
   IF( shearc < 0 )THEN       ! compute shear forces from smoothed moments
     DO i=1,3
       j = ninv(e%lnods(i))                 ! element node
       IF(factors(j) > 0d0) THEN            ! inner or boundary node
         stine(7) = stine(7) - moments(1,j)*e%b(i) + moments(3,j)*e%a(i)
         stine(8) = stine(8) - moments(3,j)*e%b(i) + moments(2,j)*e%a(i)
       ELSE                                 ! node on branching side
         ! Here it should be modified
         stine(7) = stine(7) - stint(4,iel)*e%b(i) + stint(6,iel)*e%a(i)
         stine(8) = stine(8) - stint(6,iel)*e%b(i) + stint(5,iel)*e%a(i)
       END IF
     END DO
   ELSE IF( shearc > 0 )THEN  ! compute shear forces from Gauss-points moments
     DO i=1,3
       n = e%lside(i)                 ! element side
       IF(n > 0) THEN                 ! element exist ?
         stres(1:3) = stint(4:6,n)*e%ci(i) + stint(4:6,iel)*(1-e%ci(i))
       ELSE
         stres(1:3) = stint(4:6,iel)
       END IF
       stine(7) = stine(7) + stres(1) * e%b(i) - stres(3) * e%a(i)
       stine(8) = stine(8) + stres(3) * e%b(i) - stres(2) * e%a(i)
     END DO
     stine(7:8) = 2d0*stine(7:8)
   END IF

   IF( ASSOCIATED(e%gausv) )THEN  !check if internal variables already exist
     varin = e%gausv
   ELSE
     varin = 0d0
   END IF
   pflag = .FALSE.        !initializes flag

   !Trought the thickness (t.t.t.) integration loop
   pw = 0d0                            !initializes plastic work
   DO l=1,nlayr                        !for each layer
     large = natst                     !use log strains
     IF( plast ) large = large .OR. ASSOCIATED(e%gausv)  !use log strains
     zk = thf(l)*thnew                     !Z coordinate
     stran = stra1(1:3)+stra1(4:6)*zk      !layer U^2
     IF( alpha > 0 )THEN                   !consider thermal strains
       !j0 = (EXP((t0-1d0 + zk*t1 + zk*zk*t2)))**(2d0/3d0)   !modifies jacobian due to temperature
       j0 = (t0 + zk*t1 + zk*zk*t2)**2d0     !modifies fiber length (squared) due to temperature
       stran = stran/j0                      !modifies C excluding thermal expansion
     END IF
     IF( large )THEN
       CALL lgst14(stran,r1,r2,lb,'RESVPL',error)    !Hencky (logarithmic) strains
       IF( error == 1 )THEN
         WRITE(55,"(i7,3i6,7e12.4)") e%numel,e%lnods(1:3),stra1(1:6),zk ! print to .rep file
         WRITE(55,"(3e15.4)") x                                         ! print to .rep file
         ierr = 1                       ! flag to exit program
         error = 0                      ! element flag
       END IF
     ELSE
       stran(1:2) = (stran(1:2) - 1d0 )*0.5d0
     END IF

     IF( shear )THEN           ! to include shear influence
       stres(1:2) = stint(7:8,iel)*shf(l)
     ELSE                      ! do not include shear influence
       stres(1:2) = 0d0
     END IF
     ! integrate elastic-plastic equations
     CALL stre14(stran,stres,c,propi,chib,chid,varin(:,l),ierr,3, &
                 plast,elast,val,numpt,aux,mat%matdef(3),pflag(l))

     IF( l == 1     ) stine( 9) = aux  !keep equivalent stress at bottom layer
     IF( l == nlayr ) stine(10) = aux  !and top layer

     IF( pflag(l) )THEN                !keep plastic disipated work
       aux = (varin(4,l) - aux)*vm         !plastic work (per volume unit)
       pw(1) = pw(1) + wei(l)*aux          !integrated t t thickness
       pw(2) = pw(2) + wei(l)*aux*zk       !first moment
       pw(3) = pw(3) + wei(l)*aux*zk*zk    !second moment
     END IF

     IF( large )THEN  !Compute 2nd Piola-Kirchhoff stresses
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
     ELSE
       sigma = stres  !small strain
     END IF

     !***   compute Int(B**t*sigma) on element level

     stine(1:3) = stine(1:3)+sigma*wei(l)     !t.t.t. integrated forces
     IF(shell)stine(4:6) = stine(4:6)+sigma*zk*wei(l)  !integrated moments
   END DO
   stine(1:6) = stine(1:6)*thick                  !Original thick because it's TLF

   e%pwork = pw                         !store pw for heat transfer part

   IF(ANY(pflag(1:nlayr)))THEN
     IF( updiv )THEN
       IF(.NOT.ASSOCIATED(e%gausv))  ALLOCATE(e%gausv(nvar,nlayr))
       e%gausv = varin
     ELSE
       newiv = .TRUE.
     END IF
   END IF

  !       average configuration B Matrices
   CALL bmem24(e%a,e%b,e%cd,e%a0,e%gamma,bmat(:,:,1:3),t,e%mems)
   IF(shell) CALL bfle24(e%a,e%b,e%c,e%ci,sides,isidf,bmat(:,:,4:6),t(:,1:9))

   !IF( visco )THEN
   !  CALL proma1(stran(1),v(1,1),bmat(1,1,1),1,3,18)
   !  stine(1:3) = r1*stran + stine(1:3)
   !  IF(shell) THEN
   !    CALL proma1(stran(1),v(1,1),bmat(1,1,4),1,3,18)
   !    stine(4:6) = r2*stran + stine(4:6)
   !  END IF
   !END IF

   stint(1:10,iel) = stine                      !store N, M, Q & vMS  for post-processing

   stine(1:3) = stine(1:3)*area1/3d0  ! /3 to use with average matrix
   stine(4:6) = stine(4:6)*area1      !

   IF ( branch ) THEN
     IF( nbn /= nbo )THEN
       DEALLOCATE( bb,rb )
       ALLOCATE( bb(3,nbn,6),rb(3,nbn) )
       nbo = nbn
     END IF
     CALL bbra24(bmat,e%lb,sides,nbn,e%si,bb,e%lnods(4:6))
     CALL proma3(rb(1,1),bb(1,1,1),stine(1),3*nbn,6,.TRUE.) !membrane + bending contributions

     mti = .FALSE.
     DO i=3,1,-1
       IF( sides(i) == 2 )THEN
         j = hh(3,i)         ! first node of the side
         k = hh(2,i)         ! second node of the side
         jj = 3+2*i          ! position of normal vector t at the side
         t(:,10) = x(:,k)-x(:,j)         !side J-K
         CALL vecuni(3,t(1,10),lh)       !unit side vector
         lh = e%a(i)**2+e%b(i)**2        !(ls/2A)^2 = (1/h)^2
         ! bending moment times original side length / 2 / present thickness
         aux = -(stine(4)*e%b(i)**2+stine(5)*e%a(i)**2-2d0*stine(6)*e%b(i)*e%a(i)) &
               /SQRT(lh)/thnew     !times the side length, divided by the thickness /2
         stres = x(:,i+3)*aux            !forces normal to the face
         ! torsor moment  = forces applied on nodes J and K
         aux = ((stine(4)-stine(5))*e%a(i)*e%b(i) + stine(6)*(e%b(i)**2-e%a(i)**2))/lh/area1
         lc = aux*t(:,3)                 !point loads leading to torsor moment
         rl(:,j) = rl(:,j) - lc          !modify node J forces
         rl(:,k) = rl(:,k) + lc          !modify node K forces
         !"S" forces = torsor moment * original length / 2 / present thickness
         aux = aux*SQRT(lh)*area1/thnew
         stres = stres + aux*t(:,10)     !add to node forces
         ! update shell element external nodes
         IF( mti )THEN !recompute LNA
           n = -iffix(e%lnods(j))  !compute position of the i+1 node
           lna(1:2) = nndpd(2:3,n)
           n = -iffix(e%lnods(k))  !compute position of the i+2 node
           lna(3:4) = (/ nndpd(3,n),nndpd(2,n) /)
         END IF
         t(:,1:4) = coora(:,lna(:))            !master node positions
         t(:,11) = t(:,jj)*thnew/2d0           !distance along normal
         t(:,12)  = x(:,j) - t(:,11)  !under j == 1
         t(:,13)  = x(:,j) + t(:,11)  !above j == 2
         t(:,14)  = x(:,k) + t(:,11)  !above k == 3
         t(:,15)  = x(:,k) - t(:,11)  !under k == 4
         DO n=1,4
           ! recompute proyections
           CALL projq4(t(1,n+11),t(1,1),e%cd(n,1,i),e%cd(n,2,i),sh(1),.FALSE.)
           lb = fac(n)*stres     !equivalent force
           DO j=1,4
             i3 = lna(i)          !node
             !$ CALL omp_set_lock(lock1(MOD(i3,MAXLOCK)))    !Start blocking acess in shared variables
             resid(:,i3) = resid(:,i3) + sh(j)*lb
             !$ CALL omp_unset_lock(lock1(MOD(i3,MAXLOCK)))  !End blocking acess in shared variables
           END DO
         END DO
       END IF
       mti = .TRUE.
     END DO

     DO i=1,nbn            !for each node in the patch
       n = lnb(i)          !node
       !$ CALL omp_set_lock(lock1(MOD(n,MAXLOCK)))    !Start blocking acess in shared variables
       resid(1:3,n) = resid(1:3,n)+rb(1:3,i) !sum internal force
       IF( i < 4 )THEN
         IF( bottom .OR. top )THEN
           IF(bottom) coorb(:,n) = coorb(:,n) - tt !bottom surface
           IF(top   ) coort(:,n) = coort(:,n) + tt !top surface
           ifact(n) = ifact(n) + 2  !increase weigthing factor
         END IF
       END IF
       !$ CALL omp_unset_lock(lock1(MOD(n,MAXLOCK)))  !End blocking acess in shared variables
     END DO
   ELSE
     CALL proma3(rl(1,1),bmat(1,1,1),stine(1),18,3,.TRUE.)  !membrane contributions
     ! compute bending B matrix
     IF(shell) CALL proma3(rl(1,1),bmat(1,1,4),stine(4),18,3,.FALSE.) !bending contributions
     mti = .FALSE.
     DO i=3,1,-1
       IF( sides(i) == 2 )THEN
         j = hh(3,i)         ! first node of the side
         k = hh(2,i)         ! second node of the side
         jj = 3+2*i          ! position of normal vector t at the side
         t(:,10) = x(:,k)-x(:,j)         !side J-K
         CALL vecuni(3,t(1,10),lh)       !unit side vector
         lh = e%a(i)**2+e%b(i)**2        !(ls/2A)^2 = (1/h)^2
         ! bending moment times original side length / 2 / present thickness
         aux = -(stine(4)*e%b(i)**2+stine(5)*e%a(i)**2-2d0*stine(6)*e%b(i)*e%a(i)) &
               /SQRT(lh)/thnew     !times the side length, divided by the thickness /2
         stres = x(:,i+3)*aux            !forces normal to the face
      !   ! torsor moment  = forces applied on nodes J and K
      !   aux = ((stine(4)-stine(5))*e%a(i)*e%b(i) + stine(6)*(e%b(i)**2-e%a(i)**2))/lh/area1
      !   lc = aux*t(:,3)                 !point loads leading to torsor moment
      !   rl(:,j) = rl(:,j) - lc          !modify node J forces
      !   rl(:,k) = rl(:,k) + lc          !modify node K forces
      !   !"S" forces = torsor moment * original length / 2 / present thickness
      !   aux = aux*SQRT(lh)*area1/thnew
      !   stres = stres + aux*t(:,10)     !add to node forces
         ! update shell element external nodes
         IF( mti )THEN !recompute LNA and Y
           n = -iffix(e%lnods(j))  !compute position of the i+1 node
           lna(1:2) = nndpd(2:3,n)
           n = -iffix(e%lnods(k))  !compute position of the i+2 node
           lna(3:4) = (/ nndpd(3,n),nndpd(2,n) /)
         END IF
         t(:,1:4) = coora(:,lna(:))            !master node positions
         t(:,11) = t(:,jj)*thnew/2d0           !distance along normal
         t(:,12)  = x(:,j) - t(:,11)  !under j == 1
         t(:,13)  = x(:,j) + t(:,11)  !above j == 2
         t(:,14)  = x(:,k) + t(:,11)  !above k == 3
         t(:,15)  = x(:,k) - t(:,11)  !under k == 4
         DO n=1,4
           ! recompute proyections
           CALL projq4(t(1,n+11),t(1,1),e%cd(n,1,i),e%cd(n,2,i),stra1(1),.FALSE.)
           lb = fac(n)*stres     !equivalent force
           DO j=1,4
             i3 = lna(i)          !node
             !$ CALL omp_set_lock(lock1(MOD(i3,MAXLOCK)))    !Start blocking acess in shared variables
             resid(:,i3) = resid(:,i3) + stra1(j)*lb
             !$ CALL omp_unset_lock(lock1(MOD(i3,MAXLOCK)))  !End blocking acess in shared variables
           END DO
         END DO
       END IF
       mti = .TRUE.
     END DO
     DO i=1,6             !for each node in the patch
       n = e%lnods(i)     !node
       IF(n > 0) THEN
         !$ CALL omp_set_lock(lock1(MOD(n,MAXLOCK)))    !Start blocking acess in shared variables
         resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
         IF( i < 4 )THEN
           IF( bottom .OR. top )THEN
             IF(bottom) coorb(:,n) = coorb(:,n) - tt !bottom surface
             IF(top   ) coort(:,n) = coort(:,n) + tt !top surface
             ifact(n) = ifact(n) + 2  !increase weigthing factor
           END IF
         END IF
         !$ CALL omp_unset_lock(lock1(MOD(n,MAXLOCK)))  !End blocking acess in shared variables
       END IF
     END DO
   END IF
   e => e%next          !point to next node in the list
 END DO ELEMS
 IF( nbs > 0 ) DEALLOCATE( bb,rb )
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
 IF( shearc < 0 ) THEN
   ! compute smoothed moments
   moments = 0d0
   e => heade
   iel = 1
   DO
     DO i=1,3
       j = ninv(e%lnods(i))
       moments(:,j) = moments(:,j) + stint(4:6,iel)/e%area1*factors(j)
     END DO
     e => e%next
     IF( .NOT.ASSOCIATED(e) )EXIT
     iel = iel+1
   END DO
 END IF
 RETURN
 END SUBROUTINE resv24
