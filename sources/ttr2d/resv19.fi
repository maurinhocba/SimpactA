 SUBROUTINE resv19( heade, ntype, coora, resid, ierr, ttime,  &
                    angdf, coorc, eulrf)
 !------------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT 2-D LINEAR TRIANGLE
 !                      WITH THERMO-MECHANICAL COUPLING
 !
 !------------------------------------------------------------------------
 IMPLICIT NONE

        ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ntype        !problem type
 INTEGER (kind=4), INTENT(OUT) :: ierr        !flag for error detection
 REAL(kind=8),INTENT (IN) :: coorc(:,:),    & !previous coordinates
                             coora(:,:),    & !present  coordinates
                             ttime,         & !present total time
                             angdf            !default angle X_1 - Ort_1
 LOGICAL,INTENT (IN) :: eulrf                 !user finger tensor
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
 TYPE (ele19), POINTER  :: heade              !pointer to first element

 ! local variables
 LOGICAL :: newmt,    &! TRUE if material constant computation necessary
            elast,    &! TRUE if material is strictly elastic
            plast,    &! TRUE if plasticity is to be considered
            visco,    &! TRUE if viscoplasticity is to be considered
            mises,    &! TRUE if von Mises plasticity with multiplicative dec.
            sides(3), &! TRUE if side element exist
            tplast     ! TRUE if some plastic properties is thermal dependant

 INTEGER (kind=4) :: osec,  & !associated material of previous element
                     numpt, & !number of points in curve
                     mtype, & !associated material type
                     is       !isotropic hardening model

 REAL (kind=8) :: young,nu,   & !elastic young modulus and poisson ratio
                  alpha,      & !thermal expansion coeficient
                  km,gm,      & !bulk modulus and shear modulus
                  delas(3),   & !elastic properties derivatives with respect temperature
                  fac,        & !factor for thermal strains
                  thick,      & !thickness for plane stress
                  propi(7),   & !plastic Material properties
                  propv(3),   & !Viscoplastic material properties
                  chi(12),    & !Hill-48 coeff
                  deatht        !end time for plasticity

 INTEGER (kind=4) ielem, & !element number
                  isec,  & !associated material
                  ic,    & !sector of curve for properties interpolation
                  i,n      !different indexes

 INTEGER (kind=4) ncurv(6), &    !associated curves for thermal properties
                  numct          !total number of thermal curves

 REAL (kind=8) eltmp,      & !temperature in gauss point
               dtmpe,      & !increment temperature in gauss point
               stres(4),   & !stresses (different measures)
               sigma(4),   & !stresses (different measures)
               bmat(12,4), & !B matrix
               x(2,6),     & !nodal coordinates of the element
               xa(2,6),    & !auxiliar nodal coordinates of the element
               t(2,2),     & !deformation gradient at sides
               lb(3),      & !eigenvalues
               aux,        & !auxiliar value
               c(4),       & !Elastic constitutive matrix for plane stress
               rl(2,6)       !nodal internal equivalent forces

 REAL (kind=8), POINTER :: val(:,:)

 TYPE (ele19), POINTER :: e    !pointer to an element data
 TYPE (section), POINTER :: sec    !pointer to section data
 TYPE (mater), POINTER :: mat  !pointer to an material data
 TYPE (curve), POINTER :: crvy,crvn,crvc,crvu,crvh,crva !pointers to thermal properties curves

 INTERFACE
   !INCLUDE 'stre14.h'
   !INCLUDE 'corr17.h'
   INCLUDE 'corr19.h'
 END INTERFACE

 !***  Computes nodal internal forces

 sigma(4)  = 0d0                     !initializes for plane stress
 stres(4)  = 0d0                     !initializes for plane stress
 bmat(:,4) = 0d0                     !initializes for plane stress or strain
 delas(3)  = 0d0                     !initializes
 osec = -1                           !no previous material
 e => heade                          !point to first element in the list
 DO                                  !loop over all the elements
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   ielem = e%numel                   !element number
   isec  = e%matno                   !associated material
   !     SET MATERIAL CONSTANTS
   newmt = osec /= isec             !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     sec => psecs(isec)%p
     osec = isec                    !keep present material
     mat => sec%mtbas                !point to associated material
     mtype = mat%mtype               !type of base material(ONLY MTYPE = 4, ARE ALLOWED)
     is    = mat%matdef(4)           !isotropic hardening model (1,2,3)
     mises = .FALSE.                 !initializes
     !isoch = .TRUE.                  !material is isochoric as default
     ncurv = 0                       !initialize
     ! evaluate elastic properties
     young = mat%prope(1)             !young modulus
     nu    = mat%prope(2)             !poisson ratio
     alpha = mat%prope(6)!            !Thermical Dilatation Coeff
     ! ASOCIATED CURVES - ELASTIC VALUES
     ncurv(1) = INT(sec%mtbas%props(1))        !associated curve for Young modulus
     ncurv(2) = INT(sec%mtbas%props(2))        !       "             Poisson ration
     ncurv(3) = INT(sec%mtbas%props(6))        !       "             alpha
     IF( ncurv(1) > 0 ) CALL cur_point(ncurv(1),sec%mtbas%chead,crvy) !pointer to Young curve
     IF( ncurv(2) > 0 ) CALL cur_point(ncurv(2),sec%mtbas%chead,crvn) !pointer to Poisson curve
     IF( ncurv(3) > 0 ) CALL cur_point(ncurv(3),sec%mtbas%chead,crva) !pointer to Alpha curve
     ! evaluates plastic properties
     IF( mat%matdef(3) /= 1 )THEN
       mises = mat%matdef(3) == 2     !if von Mises yield function
       elast = .FALSE.                !not elastic
       deatht = mat%propp(4)          !end of plasticity
       plast = deatht > ttime         !consider plasticity?
       IF( is == 3 )THEN    !Ludwick-Nadai
         propi(1:3) = mat%propp(1:3)               !isotropic hardening parameters
         aux = (propi(2)/propi(1))**(1d0/propi(3))
         propi(2) = aux !evaluate EFREF
       ELSEIF( is == 2)THEN !Linear
         propi(1:2) = mat%propp(1:2)               !isotropic hardening parameters
       ELSE                 !None
         propi(1) = mat%propp(1)    !isotropic hardening parameters
       END IF
       propi(5) = REAL(is,8)          !isotropic hardening model
       NULLIFY (val)                  !no point-defined isotropic hardening
       !IF (mat%matdef(6) == 3) THEN !elastic-viscoplastic material   ... not yet
       !  visco = .TRUE.
       !  propv(1:3) = mat%propp(28:30) !viscoplastic parameters
       !ELSE
       visco = .FALSE.
       propv = 0d0
       !END IF
       chi = mat%propp(16:27)         !Hill yield function coefficients
       ! ASOCIATED CURVES - PLASTIC VALUES
       ncurv(4) = INT(sec%mtbas%props(3))        !associated curve for CONS0
       ncurv(5) = INT(sec%mtbas%props(4))        !         "           UNIAX
       ncurv(6) = INT(sec%mtbas%props(5))        !         "           CONSN (or HARDS)
       IF( ncurv(4) > 0 ) CALL cur_point(ncurv(4),sec%mtbas%chead,crvc) !pointer to CONS0 curve
       IF( ncurv(5) > 0 ) CALL cur_point(ncurv(5),sec%mtbas%chead,crvu) !pointer to UNIAX curve
       IF( ncurv(6) > 0 ) CALL cur_point(ncurv(6),sec%mtbas%chead,crvh) !pointer to CONSN (or HARDS) curve
     ELSE
       elast = .TRUE.
       plast = .FALSE.
     END IF
     !  elasticity constants
     IF( ntype == 1 )THEN      !plane stress
       c(1:4) = mat%prope(10:13)  ! plane stress elasticity matrix
       thick  = sec%rprop(1)      ! thickness for plane stress
     ELSE                      ! plane strain or axilsymmetric
       km = mat%prope(4)*3d0      ! bulk modulus * 3
       gm = mat%prope(3)          ! shear modulus
     END IF
     numct  = SUM(ncurv(1:3)) !elastic temperature dependat properties
     tplast = plast .AND. (SUM(ncurv(4:6)) > 0) !plastic temp. dependant prop.
   END IF

   ! if problem is thermal-elastic/plastic then evaluates properties
   eltmp = SUM(tempe(1,e%lnods(1:3)))/3d0 !average temperature
   dtmpe = SUM(dtemp(1,e%lnods(1:3)))/3d0 !average incremental temperature
   IF( numct > 0 )THEN !temperature dependant elastic properties
     ic = 1
     ! elastic properties
     IF( ncurv(1) > 0 .OR. ncurv(2) > 0 ) THEN
       IF( ncurv(1) > 0 )THEN
         young    = inte_cr(crvy%val,crvy%np,eltmp,ic)
         delas(1) = crvy%val(3,ic) !dyoung/dT
       END IF
       IF( ncurv(2) > 0 )THEN
         nu       = inte_cr(crvn%val,crvn%np,eltmp,ic)
         delas(2) = crvn%val(3,ic) !dnu/dT
       END IF
       IF( ntype == 1 )THEN     !plane stress
         c(1) = young/(1d0-nu*nu)    ! plane stress elasticity matrix
         c(2) = c(1)*nu              !
         c(3) = c(1)                 !
         c(4) = young/(2d0*(1d0+nu)) !
       ELSE                     ! plane strain or axilsymmetric
         gm = young/(1d0+nu)/2d0      ! shear modulus
         km = young/(1d0-2d0*nu)      ! bulk modulus * 3
         delas(3) = (delas(1)*(1d0+nu)-young*delas(2)) &
                    / (1d0+nu)**2 / 2d0 !dgm/dT
       END IF
     END IF
     IF( ncurv(3) > 0 ) alpha = inte_cr(crva%val,crva%np,eltmp,ic)!/3d0
   END IF
   ! if temperature dependant plastic properties
   IF( tplast )THEN
     ic = 1
     IF( is == 3 )THEN    !Ludwick-Nadai
       IF( ncurv(4) > 0 ) propi(1) = inte_cr(crvc%val,crvc%np,eltmp,ic)
       IF( ncurv(5) > 0 ) propi(2) = inte_cr(crvu%val,crvu%np,eltmp,ic)
       IF( ncurv(6) > 0 ) propi(3) = inte_cr(crvh%val,crvh%np,eltmp,ic)
       aux = (propi(2)/propi(1))**(1d0/propi(3))
       propi(2) = aux !evaluate EFREF
     ELSEIF( is == 2)THEN !Linear
       IF( ncurv(5) > 0 ) propi(1) = inte_cr(crvu%val,crvu%np,eltmp,ic)
       IF( ncurv(6) > 0 ) propi(2) = inte_cr(crvh%val,crvh%np,eltmp,ic)
     ELSE                 !None
       IF( ncurv(5) > 0 ) propi(1) = inte_cr(crvu%val,crvu%np,eltmp,ic)
     END IF
   END IF

   IF(eulrf)THEN ! <--- spatial formulation
     !*** computes configuration in previous and present steps
     x(:,1:3) = coorc(:,e%lnods(1:3))  !previous coordinates of vertex nodes
     xa(:,1:3) = coora(:,e%lnods(1:3)) !present coordinates of vertex nodes
     DO i =4,6                         !for the rest of nodes in the patch
       n = e%lnods(i)                  !node number
       IF(n /= 0) THEN                 !if node exist
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
     CALL cdac19(e%numel,e%lnods,xa,e%cd,e%area1,ierr,thick,angdf,e%angle,ntype)
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
   fac = 1d0 !factor for thermal strains for ntype = 1
   IF( ntype > 1 ) fac   = 1d0 + nu       ! factor for thermal plane strains and axylsimmetric
   t(:,:) = MATMUL(x,e%cd(:,1:2))  ! in-plane derivatives
   IF( ntype == 2 ) lb(3) = 1d0 ! transversal stretch
   IF( ntype == 3 ) lb(3) = DOT_PRODUCT(x(1,:),e%cd(:,3))  ! hoop stretch

   IF( mises .AND. ntype /= 1 )THEN
     CALL corr19(t,lb(3),e%gausv,is,propi,propv,gm,km,e%stint,sigma, &
                 numpt,val,ierr,e%pwork,e%strpl,delas(3),fac,alpha,dtmpe)

   !ELSE
   !
   ! CALL stra19(t,ntype,lb(3),stran)  !metric tensor at center
   ! !computes Hencky (logarithmic) in-plane strains
   ! CALL lgst14(stran,r1,r2,lb,'RESV19',ierr)
   ! stran(4) = LOG(lb(3))
   ! IF( ierr == 1 )THEN
   !   WRITE(55,"(4i5,4e12.4)",ERR=9999) e%lnods(1:4),stran(1:4)
   !   WRITE(55,"(2e15.4)",ERR=9999) x
   !   CALL runen3('RESV19: negative eigenvalues of U^2')
   ! END IF
   !
   ! IF( ntype == 1 )THEN  !plane stress
   !
   !   IF( mtype == 1 )THEN          !one component material
   !     stres(1:2) = 0d0
   !     CALL stre14(stran,stres,c,propi,chi,e%gausv,ierr,4,newmt,plast,elast,&
   !                 val,numpt,aux,mat%matdef(3),pflag)
   !   END IF
   !
   ! ELSE    !(ntype == 2 or 3) Plane Strain or Axisymmetric
   !
   !   stran(4) = LOG(lb(3))
   !   IF( mtype == 1 )THEN          !one component material
   !     ! compute volumetric and deviatoric parts
   !     delta = ( stran(1) + stran(2) + stran(4) )/3d0
   !     stran(1) = 2d0*(stran(1) - delta)  !twice the deviatoric shear strain
   !     stran(2) = 2d0*(stran(2) - delta)
   !     stran(4) = 2d0*(stran(4) - delta)
   !     !elastic (trial) strains
   !     IF( .NOT.elast )THEN
   !       strpl(1:4) = e%gausv(1:4)               !previous (twice) plastic strains
   !       stran(1:4) = stran(1:4) - strpl(1:4)    !trial Elastic log strains
   !     END IF
   !     stres= gm*stran                         !Trial elastic shear stresses
   !     IF( plast ) THEN            !Consider plasticity
   !       efpst = e%gausv(5)        !effect plastic strain
   !       CALL corr17(stres(1),stres(2),stres(3),stres(4),efpst,gm, &
   !                   propi,fr,ierr,dstpl,pflag,is,newmt,numpt,val)
   !       IF(ierr == 1) RETURN              !no convergence in plasticity
   !       IF( pflag )THEN                   !if plastic flow
   !         e%gausv(1:4) = e%gausv(1:4) + dstpl(1:4)  !total plastic shear strains
   !         e%gausv(5)   = e%gausv(5)   + efpst       !Total effect. plastic strain
   !       END IF
   !     END IF
   !
   !   END IF
   ! END IF
   !
   ! IF( ntype /= 1 )THEN  !only for one-component material
   !   IF( mtype == 1 )THEN
   !     press = km*delta
   !     stres(1) = stres(1) + press
   !     stres(2) = stres(2) + press
   !     stres(4) = stres(4) + press
   !   END IF
   ! END IF
   !
   ! ! Computes Hencky stress on the natural Frame
   ! sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2+2d0*stres(3)*r1*r2
   ! sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1-2d0*stres(3)*r1*r2
   ! sigma(3) =(stres(2)-stres(1))*r1*r2+stres(3)*(r1*r1-r2*r2)
   ! ! Computes 2nd P-K stress on the natural Frame
   ! stres(1) = sigma(1)/lb(1)**2
   ! stres(2) = sigma(2)/lb(2)**2
   ! IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
   !   stres(3) = sigma(3)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
   ! ELSE                                !lb(1) = lb(2)
   !   stres(3) = sigma(3)/lb(1)/lb(2)
   ! END IF
   ! ! Computes 2nd P-K on the Lagrangian Frame
   ! sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2-2d0*stres(3)*r1*r2
   ! sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1+2d0*stres(3)*r1*r2
   ! sigma(3) =(stres(1)-stres(2))*r1*r2+stres(3)*(r1*r1-r2*r2)
   ! IF( ntype /= 1 )sigma(4) = stres(4)/lb(3)**2
   ! ! compute Kirchhoff stress tensor for post-process
   ! e%stint(1) = sigma(1)*t(1,1)*t(1,1)+sigma(2)*t(1,2)*t(1,2)+2d0*sigma(3)*t(1,2)*t(1,1)
   ! e%stint(2) = sigma(1)*t(2,1)*t(2,1)+sigma(2)*t(2,2)*t(2,2)+2d0*sigma(3)*t(2,1)*t(2,2)
   ! e%stint(3) = sigma(1)*t(1,1)*t(2,1)+sigma(2)*t(1,2)*t(2,2)+sigma(3)*(t(1,1)*t(2,2)+t(1,2)*t(2,1))
   ! e%stint(4) = stres(4)

   END IF

   !***   compute Int(B**t*sigma) on element level
   sigma = sigma*e%area1         !Original dvol used (TLF) or updated (ULF)

   !       actual configuration B Matrx
   CALL bmat19(e%cd,bmat,t,sides,ntype,lb(3))

   CALL proma3(rl(1,1),bmat(1,1),sigma(1),12,4,.TRUE.)  !nodal contributions

   DO i=1,6             !for each node in the element
     n = e%lnods(i)     !node
     IF( n > 0 ) resid(1:2,n) = resid(1:2,n)+rl(1:2,i) !sum internal force
   END DO
   e => e%next          !point to next element in the list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE resv19
