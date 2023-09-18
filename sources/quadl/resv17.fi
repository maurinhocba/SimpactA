 SUBROUTINE resv17( heade, ntype, coora, resid, ierr, ttime, coord, ngaus)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT T-L 2-D QUAD
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ntype, &     !problem type
                                 ngaus        !number of gauss points
 INTEGER (kind=4), INTENT(OUT) :: ierr        !flag for error detection
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             coord(:,:),    & !original coordinates
                             ttime            !present total time
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
 TYPE (ele17), POINTER  :: heade              !pointer to first element

 ! local variables
 LOGICAL :: newmt,   &! TRUE if material constant computation necessary
            elast,   &! TRUE if material is strictly elastic
            plast     ! TRUE if plasticity is to be considered
 INTEGER (kind=4) :: osec,  & !associated material of previous element
                     numpt, & !number of points in curve
                     mtype, & !associated material type
                     is       !isotropic hardening model
 REAL (kind=8) :: fr(7),      & !Flow rule coefficients
                  km,gm,      & !bulk modulus and shear modulus
                  propi(7),   & !Material properties
                  chi(12),chib(4),  & !Hill-48 coeff
                  deatht        !end time for plasticity
 INTEGER (kind=4) ielem, & !element number
                  isec,  & !associated material
                  i,j,k,n  !different indexes

 REAL (kind=8) stres(4),   & !stresses (different measures)
               sigma(4),   & !stresses (different measures)
               stran(4),   & !C=U^2  also Log strains
               dstpl(4),   & !plastic strain increments
               efpst,      & !effective plastic strain
               strpl(4),   & !plastic strain (previous step)
               bmat(8,4),  & !B matrix
               x(2,4),     & !nodal coordinates of the element
               x0(4),      & !nodal original 'r' coordinates
               g(2,2),     & !in-plane deformation gradient
               r0,r,       & !original and current 'r' coordinate
               r1,r2,      & !eigevector components in local system
               lb(3),      & !eigenvalues
               aux,        & !auxiliar value
               c(4),       & !Elastic constitutive matrix for plane stress
               delta,      & !volumetric strain/3 and Jacobian
               deltm,      & !Average volumetric strain
               press,      & !first stress invariant (mean pressure)
               rl(2,4)       !nodal internal equivalent forces
 REAL (kind=8) :: rv(2,ngaus),lbv(3,ngaus),str(4,ngaus),gv(2,2,ngaus),rf(ngaus)

 INTEGER (kind=4), PARAMETER :: ngaud = 2
 REAL (kind=8) :: shape(nnode,ngaus),deriv(nnode,2),posgp(ngaud),weigp(ngaud)
 REAL (kind=8), POINTER :: val(:,:)
 LOGICAL :: pflag

 TYPE (ele17), POINTER :: e    !pointer to an element data
 TYPE (section), POINTER :: sec    !pointer to section data
 TYPE (mater), POINTER :: mat  !pointer to an material data

 INTERFACE
   INCLUDE 'stre14.h'
   INCLUDE 'corr17.h'
   INCLUDE 'rubberps.h'
   INCLUDE 'rubber2d.h'
 END INTERFACE

 IF( ntype == 3 )THEN
   CALL gaussq (ngaud ,posgp ,weigp )

   !       gauss points shape and derivatives of nodal functions

   k = 0
   DO i=1,ngaud
     DO j=1,ngaud
       k = k+1
       CALL shape3(deriv(1,1),shape(1,k),posgp(i),posgp(j),nnode)
     END DO
   END DO
 END IF

 !***  Computes nodal internal forces

 sigma(4) = 0d0                      !initializes for plane stress
 stres(4) = 0d0                      !initializes for plane stress
 bmat(:,4) = 0d0                     !initializes for plane stress or strain
 osec = -1                           !no previous material
 e => heade                          !point to first element in the list
 DO                                  !loop over all the elements
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   ielem = e%numel                   !element number
   isec = e%matno                   !associated material
   !     SET MATERIAL CONSTANTS
   newmt = osec /= isec             !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     sec => psecs(isec)%p
     osec = isec                    !keep present material
     mat => sec%mtbas                !point to associated material
     mtype = mat%mtype               !type of base material
     is   = mat%matdef(4)           !isotropic hardening model
     ! ONLY MTYPES 1, 5, 6) ARE ALLOWED
     IF( mtype == 1)THEN              !Isotropic Elastic - Anisotropic Plastic
       IF( mat%matdef(3) /= 1 )THEN
         elast = .FALSE.                !not elastic
         deatht = mat%propp(5)          !end of plasticity
         plast = deatht > ttime         !consider plasticity ?
         propi(1:7) = mat%propp(1:7)    !isotropic hardening parameters
         propi(5) = REAL(is,8)          !isotropic hardening model
         chi = mat%propp(16:27)         !Hill yield function coefficients
         IF ( is == 5 )THEN             !points defined yield curve
           val => mat%chead%val
           numpt = mat%chead%np
         ELSE
           NULLIFY (val)
           numpt = 0
         END IF
         ! for anisotropic (PLANAR) plasticity
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
         ! A matrix, derivative of yield function (associative only)
         CASE ( 2 )   !plane strain
           fr(1) =  chi(8) + chi(9) + 0.5d0  !0.5d0 added to make it diagonal
           fr(2) = -chi(8)          + 0.5d0  ! in case of von Mises
           fr(3) =  chi(7) + chi(8) + 0.5d0  ! it is not essential
           fr(4) = 2d0*chi(11)
           fr(5) = -chi(9)          + 0.5d0
           fr(6) = -chi(7)          + 0.5d0
           fr(7) =  chi(7) + chi(9) + 0.5d0
         CASE ( 3 )  !axilsymmetric
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
     ELSE IF( mtype == 6)THEN            !hyperelastic isotropic
       chi(1:12) = mat%prope(7:18)
       elast = .TRUE.
       plast = .FALSE.                   !consider plasticity ?
       !km = mat%prope(4)*3d0               ! bulk modulus * 3
       !gm = mat%prope(3)                   ! initial shear modulus
       !propi(1:2) = mat%propp(1:2)         ! plastic properties
       !plast = propi(1) > 0 ! .AND. propp(5) > ttime
     END IF
     newmt = .FALSE.
   END IF

   x = coora(:,e%lnods)                   !*** computes actual configuration
   IF( ntype == 3 ) x0 = coord(1,e%lnods) !compute original 'r' coordinates
   rl = 0d0                               !Initializes equivalent nodal forces
   deltm = 0d0                            !initializes meam volumetric strain
   !   Gauss point loop
   DO k=1,ngaus
     ! compute deformation gradient F and U^2=tp(F) F
     g = MATMUL(x,e%cartd(:,:,k))  !deformation gradient
     stran(1) = g(1,1)*g(1,1) + g(2,1)*g(2,1)  !U^2(1,1)
     stran(2) = g(1,2)*g(1,2) + g(2,2)*g(2,2)  !U^2(2,2)
     stran(3) = g(1,1)*g(1,2) + g(2,1)*g(2,2)  !U^2(1,2)
     SELECT CASE (ntype)
     CASE (1)  !   plane stress
       stran(4) = 1d0/(stran(1)*stran(2) - stran(3)**2)  !isochoric
     CASE (2)  !   plane strain
       stran(4) = 1d0
     CASE (3)  !   axilsymmetric
       r0 = DOT_PRODUCT(x0,shape(:,k))         !original r coordinate
       r  = DOT_PRODUCT(x(1,:),shape(:,k))     !present r coordinate
       stran(4) = (r/r0)**2                    !Lambda 3 squared
       rf(k) = stran(4)/r                      !Coefficient for B matrix
     END SELECT
     !computes Hencky (logarithmic) in-plane strains
     CALL lgst14(stran,r1,r2,lb,'RESV17',ierr)
     IF( ierr == 1 )THEN
       WRITE(55,"(4i5,4e12.4)",ERR=9999) e%lnods(1:4),stran(1:4)
       WRITE(55,"(2e15.4)",ERR=9999) x
       CALL runen3('RESV17: negative eigenvalues of U^2')
     END IF
     lb(3) = SQRT(stran(4))
     stran(4) = LOG(lb(3))

     IF( ntype == 1 )THEN            !plane stress
       IF( mtype == 1 )THEN          !one component material
         stres(1:2) = 0d0
         CALL stre14(stran,stres,c,propi,chib,fr,e%gausv(:,k),ierr,4, &
                     plast,elast,val,numpt,aux,mat%matdef(3),pflag)
       ELSE IF ( mtype == 5 )THEN ! Orthotropic elastic
         stres(1) = c(1)*stran(1)+c(2)*stran(2)
         stres(2) = c(2)*stran(1)+c(3)*stran(2)
         stres(3) = c(4)*stran(3)
       ELSE IF ( mtype == 6 )THEN ! Hyperelastic
         !IF(plast)THEN
           !!     elastic-plastic behavior (NEW) isochoric
           !efpst = e%gausv(5,k)                !effect plastic strain
           !IF( efpst > 0 )THEN
           !  !elastic (trial) strains
           !  strpl(1:4) = e%gausv(1:4,k)            !previous (twice) plastic strains
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
           !    e%gausv(5,k) = efpst      !effective plastic strain
           !    e%gausv(1:4,k) = stran(1:4) - dstpl(1:4)   !plastic  log strains
           !  END IF
           !END IF
           !
           !! compute (Total) stress tensor
           !stres(1) = s1*bet(1)*s1 + s2*bet(2)*s2 - bet(3)  !(1,1)
           !stres(2) = s2*bet(1)*s2 + s1*bet(2)*s1 - bet(3)  !(2,2)
           !stres(3) = s1*bet(1)*s2 - s2*bet(2)*s1           !(1,2)
           !stres(4) = 0d0                                   !(3,3)
        !ELSE
          CALL rubberps(chi,lb,mat%matdef(8),stres,r1=r1,r2=r2)
        !END IF

       END IF
     ELSE    !(ntype == 2 or 3) Plane Strain or Axisymmetric
       IF( mtype == 1 )THEN          !one component material
         ! compute volumetric and deviatoric parts
         delta = ( stran(1) + stran(2) + stran(4) )/3d0
         deltm = deltm + delta*e%dvol(k)     !use mean volumetric strain
         stran(1) = 2d0*(stran(1) - delta)  !twice the deviatoric shear strain
         stran(2) = 2d0*(stran(2) - delta)
         stran(4) = 2d0*(stran(4) - delta)
         !elastic (trial) strains
         IF( .NOT.elast )THEN
           strpl(1:4) = e%gausv(1:4,k)             !previous (twice) plastic strains
           stran(1:4) = stran(1:4) - strpl(1:4)    !trial Elastic log strains
         END IF
         stres= gm*stran                         !Trial elastic shear stresses
         IF( plast ) THEN              !Consider plasticity
           efpst = e%gausv(5,k)        !effect plastic strain
           CALL corr17(stres(1),stres(2),stres(3),stres(4),efpst,gm, &
                       propi,fr,ierr,dstpl,pflag,is,ielem,numpt,val)
           IF(ierr == 1) RETURN              !no convergence in plasticity
           IF( pflag )THEN                      !if plastic flow
             e%gausv(1:4,k) = e%gausv(1:4,k) + dstpl(1:4)  !total plastic shear strains
             e%gausv(5,k) = e%gausv(5,k) + efpst          !Total effect. plastic strain
           END IF
         END IF
       ELSE IF ( mtype == 5 )THEN ! Orthotropic elastic
         stres(1) = propi(1)*stran(1)+propi(2)*stran(2)+propi(3)*stran(3)
         stres(2) = propi(2)*stran(1)+propi(4)*stran(2)+propi(5)*stran(3)
         stres(3) = propi(6)*stran(3)
         stres(4) = propi(3)*stran(1)+propi(5)*stran(2)+propi(7)*stran(3)
       ELSE IF ( mtype == 6 )THEN ! Hyperelastic
         !IF( plast )THEN
           !efpst = e%gausv(5,k)                !effect plastic strain
           !IF( efpst > 0 )THEN
           !  !elastic (trial) strains
           !  strpl(1:4) = e%gausv(1:4,k)            !previous (twice) plastic strains
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
           !    e%gausv(5,k) = efpst      !effective plastic strain
           !    e%gausv(1:4,k) = stran(1:4) - dstpl(1:4)   !plastic  log strains
           !  END IF
           !END IF
           !
           !! compute (deviatoric) stress tensor
           !stres(1) = s1*bet(1)*s1 + s2*bet(2)*s2    !(1,1)
           !stres(2) = s2*bet(1)*s2 + s1*bet(2)*s1    !(2,2)
           !stres(4) = bet(3)                         !(3,3)
           !stres(3) = s1*bet(1)*s2 - s2*bet(2)*s1    !(1,2)
         !ELSE
          CALL rubber2d(chi,lb,mat%matdef(8),stres,r1,r2)
         !END IF
       END IF
     END IF
     ! store auxiliar values for second loop
     str(:,k)  =  stres
     rv(1,k)   =  r1
     rv(2,k)   =  r2
     lbv(:,k)  =  lb
     gv(:,:,k) =  g
   END DO

   IF( ntype /= 1 )THEN  !only for one-component material
     IF( mtype == 1 )THEN
       deltm = deltm/SUM(e%dvol)
       press = km*deltm
       str(1,:) = str(1,:) + press
       str(2,:) = str(2,:) + press
       str(4,:) = str(4,:) + press
     ELSE
       ! nothing yet
     END IF
   END IF
   DO k=1,ngaus         !second Gauss loop (for one component material)
     ! recovers Gauss point values
     stres = str(:,k)   !stresses
     e%stint(:,k) = stres   !average Kirchhoff stress
     IF( mtype == 6 )THEN
       sigma = stres        !2nd Piola Kirchhoff
     ELSE
       r1 = rv(1,k)       !first component of first eigenvector
       r2 = rv(2,k)       !second component of first eigenvector
       lb = lbv(:,k)      !eigenvalues (principal stretchings)
       g  = gv(:,:,k)     !metric tensor
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
    END IF
    !***   compute Int(B**t*sigma) on element level

     sigma = sigma*e%dvol(k)     !Original dvol used (TLF)

     !       actual configuration B Matrices
     CALL bmat17(e%cartd(:,:,k),bmat,g,nnode,rf(k),shape(:,k),ntype)
     CALL proma3(rl(1,1),bmat(1,1),sigma(1),8,4,.FALSE.)  !nodal contributions
   END DO

   DO i=1,4             !for each node in the element
     n = e%lnods(i)     !node
     resid(1:2,n) = resid(1:2,n)+rl(1:2,i) !sum internal force
   END DO

   e => e%next          !point to next node in the list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE resv17
