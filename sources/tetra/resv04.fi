 SUBROUTINE resv04( heade, coora, resid, ierr, ttime, small)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT T-L 3-D SOLID
 !
 !-------------------------------------------------------------------
 USE ele18_db, ONLY : rubber3d,lgst18,eige18
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(OUT) :: ierr        !flag for error detection
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             ttime            !present total time
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
 TYPE (ele04), POINTER  :: heade              !pointer to first element
 LOGICAL, INTENT(IN) :: small      !.TRUE. if Green strains will be used

 ! local variables
 LOGICAL, SAVE :: newmt,   & ! TRUE if material constant computation necessary
                  elast,   & ! TRUE if material is elastic
                  plast      ! TRUE if plasticity is to be considered
 INTEGER (kind=4), SAVE :: osec = -1,  & !associated material of previous element
                           numpt,      & !number of points in curve
                           mtype         !associated material type
 REAL (kind=8), SAVE ::  fr(9),      & !Flow rule coefficients
                         km,gm,      & !bulk modulus and shear modulus
                         eprop(12),  & !elastic Material properties
                         props(7),   & !plastic Material properties
                         ch(6)         !Hill-48 coefficients

 INTEGER (kind=4) ielem,             & !element number
                  is,                & !isotropic hardening model
                  isec,              & !associated material
                  eflag,             & !flag for eigen-decomposition
                  i,n          !different indexes

 INTEGER (kind=4) nl,nlayr,  & !number of layers
                  oldm,      & !old material label
                  mate,      & !material label
                  ilv,jlv,   & !indexes for layer varin
                  il,        & !layer index
                  nvar         !total number of internal variables
 INTEGER (kind=4), ALLOCATABLE ::  lvari(:)
 REAL (kind=8), ALLOCATABLE ::  thickl(:),zposl(:),rr(:,:)
 REAL (kind=8), POINTER :: varin(:)                      !internal variables
 REAL (kind=8) stral(6),  & !local strains
               strel(6),  & !local streses
               staux(6),  & !local streses auxiliar
               thick,     & !thickness (original)
               minstr,    & !minimum strain to integrate trougth the thickness
               dm(36),    & !Elastic integrated constitutive matrix
               aux          !auxiliar value

 REAL (kind=8) str(6),     & !stresses (different measures)
               sig(6),     & !stresses (different measures)
               stran(6),   & !C=U^2  also Log strains
               dstpl(6),   & !plastic strain increments
               efpst,      & !effective plastic strain
               strpl(6),   & !plastic strain (previous step)
               bmat(12,6), & !B matrix
               x(3,nnode), & !nodal coordinates of the element
               g(3,3),     & !deformation gradient
               r(3,3),lb(3),& !eigen-pair in local system
               delta,      & !volumetric strain/3, element volume
               press,      & !first stress invariant (mean pressure)
               rl(3,nnode)   !nodal internal equivalent forces
 !REAL (kind=8) :: rv(3,3,ngaus),lbv(3,ngaus),stres(6,ngaus),gv(3,3,ngaus)
 LOGICAL :: pflag,   &! TRUE if plastic flow in the step
            ttti,    &! TRUE if Trougth The Thickness Integration is necessary
            flag2     ! TRUE if plastic at any gauss point


 REAL (kind=8), POINTER :: val(:,:)
 REAL (kind=8)          :: elen    ! Damage length - For damage subroutine - XMG
 TYPE (ele04), POINTER :: e      !pointer to an element data
 TYPE (section), POINTER, SAVE :: sec  !pointer to an section data
 TYPE (mater), POINTER, SAVE :: mat    !pointer to an material data

 TYPE (pmater), POINTER :: matsl(:,:)    !array of pointers to materials data


 INTERFACE
   INCLUDE 'corr18.h'
 END INTERFACE


 !***  Computes nodal internal forces

 osec = -1                           !no previous material
 e => heade                          !point to first element in the list
 DO                                  !loop over all the elements
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   ielem = e%numel                   !element number
   isec = e%matno                    !associated material
   !     SET MATERIAL CONSTANTS
   newmt = osec /= isec              !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     osec = isec                     !keep present material
     sec => psecs(isec)%p             !point to section
     nlayr = sec%iprop(1)             !number of layers
     nvar  = sec%iprop(2)             !total number of internal variables

     !************************** STANDARD SOLID *****************************************
     IF( nlayr == 0 )THEN            !standard solid section (not layered)

       mat => sec%mtbas                !point to associated material
       mtype = mat%mtype               !type of base material
       elast = mat%matdef(3) == 1      !elastic material
       is    = mat%matdef(4)           !isotropic hardening model
       ! ------------------------------------------------------------------------------------
       ! ONLY MTYPES 1, 5, 6, 7 and 92 ARE ALLOWED
       SELECT CASE (mtype)
       ! ------------------------------------------------------------------------------------
       CASE (1)                        !Isotropic Elastic - Anisotropic Plastic  ( mtype == 1)
         km = mat%prope(4)*3d0              ! bulk modulus * 3
         gm = mat%prope(3)                  ! shear modulus
         IF( .NOT.elast ) THEN
           ! e1, nu1, uniaxial, efren, consn, r, exponent m, hill 79
           props(1:7) = mat%propp(1:7)        !
           ! for anisotropic (PLANAR) plasticity
           ch = mat%propp(16:21)              !Hill yield function coefficients
           ! A matrix, derivative of yield function (associative only)
           fr(1) =  ch(2)+ch(3) + 0.5d0       !0.5d0 added to make it diagonal
           fr(2) = -ch(3)       + 0.5d0       ! in case of von Mises
           fr(3) = -ch(2)       + 0.5d0       ! it is not essential
           fr(4) =  ch(1)+ch(3) + 0.5d0
           fr(5) = -ch(1)       + 0.5d0
           fr(6) =  ch(1)+ch(2) + 0.5d0
           fr(7) = 2d0*ch(6)
           fr(8) = 2d0*ch(4)
           fr(9) = 2d0*ch(5)
           plast = props(5) > ttime                 !consider plasticity ?
           IF ( is == 5 )THEN             !points defined yield curve
             val => mat%chead%val
             numpt = mat%chead%np
           ELSE
             NULLIFY (val)
             numpt = 0
           END IF
         ELSE
           plast = .FALSE.
         END IF
         ! ------------------------------------------------------------------------------------
       CASE (5)                        !Orthotropic Elastic - No Plasticity      ( mtype == 5)
         eprop(1) = mat%prope(45)      !(1,1)   11-11
         eprop(2) = mat%prope(46)      !(2,1)   11-22
         eprop(3) = mat%prope(48)      !(3,1)   11-33
         eprop(4) = mat%prope(52)      !(2,2)   22-22
         eprop(5) = mat%prope(54)      !(3,2)   22-33
         eprop(6) = mat%prope(66)      !(3,3)   33-33
         eprop(7) = mat%prope(59)      !(4,4)   12-12
         eprop(8) = mat%prope(73)      !(5,5)   13-13
         eprop(9) = mat%prope(80)      !(6,6)   23-23
         plast = .FALSE.
       ! ------------------------------------------------------------------------------------
       CASE (6)                        !Hyperelastic isotropic rubber            ( mtype == 6)
         eprop(1:12) = mat%prope(7:18)
         elast = .TRUE.
         plast = .FALSE.                   !consider plasticity ?
         !km = mat%prope(4)*3d0               ! bulk modulus * 3
         !IF( .NOT.elast ) THEN
         !  props(1:2) = mat%propp(1:2)         ! plastic properties
         !  plast = props(1) > 0 ! .AND. props(5) > ttime
         !ELSE
         !  plast = .FALSE.
         !END IF
       ! ------------------------------------------------------------------------------------
       END SELECT

       !************************** LAYERED SOLID *************************************************
     ELSE  !nlayr >= 1      layered solid (at integration point level)
       thick = sec%rprop(1)             !original thickness
       minstr= sec%rprop(4)             !strain threshold to use TTTI

       nl=1
       elast = sec%iprop(3) == 0               !elastic problem
       plast = .NOT.elast                      !consider plasticity ?
       dm = sec%rprop(6:41)                 !linear elastic integrated matrix
       IF( plast )THEN  !plastic lamina
         IF( ALLOCATED ( rr ) )DEALLOCATE( thickl,zposl,rr,matsl,lvari )
         ALLOCATE( thickl(nlayr),zposl(nlayr),rr(5,nlayr),matsl(3,nlayr),lvari(nlayr) )
         i = 42                                !pointer to layer thickness
         CALL vecasi(nlayr,sec%rprop(i),thickl(1))        !layer thickness
         i = i+2*nlayr                         !pointer to layer z position
         !CALL vecasi(nlayr,sec%rprop(i),zposl(1))         !layer z position
         i = i+nlayr                           !pointer to layer rotation matrix
         CALL vecasi(nlayr*5,sec%rprop(i),rr(1,1))        !layer rotation matrix

         lvari(1:nlayr)= sec%iprop(nlayr+5:2*nlayr+4)
         DO il=1,nlayr
           mate= sec%iprop(4+il)
           mat => pmats(mate)%p        !search layer material
           matsl(1,il)%p => mat
           SELECT CASE (mat%mtype)
           CASE ( 91,92 )
             !point to component materials data bases
             matsl(2,il)%p => pmats(mat%matdef(4))%p  !matm
             matsl(3,il)%p => pmats(mat%matdef(5))%p  !matf
           CASE DEFAULT
             NULLIFY( matsl(2,il)%p )
             NULLIFY( matsl(3,il)%p )
           END SELECT
         END DO
       END IF
       oldm = -1                              !recompute constants
     END IF
     newmt = .FALSE.                        !same material than previous ?

   END IF !newmt


   !*** computes actual configuration
   x = coora(:,e%lnods(1:4))
   g = MATMUL(x,e%cd) !x(1:3) Main element derivatives (deformation gradient)
   !       actual configuration B Matrices
   CALL bmat04(e%cd,bmat,g,nnode)

   !   check
   IF( ASSOCIATED(e%gausv) )THEN
     varin => e%gausv
   ELSE
     ALLOCATE( varin(nvar) )  !
     varin = 0d0
     flag2 = .FALSE.
   END IF
   ! Rigth Cauchy-Green Tensor
   stran(1) = g(1,1)*g(1,1) + g(2,1)*g(2,1) + g(3,1)*g(3,1)  !U^2(1,1)
   stran(2) = g(1,2)*g(1,2) + g(2,2)*g(2,2) + g(3,2)*g(3,2)  !U^2(2,2)
   stran(3) = g(1,3)*g(1,3) + g(2,3)*g(2,3) + g(3,3)*g(3,3)  !U^2(3,3)
   stran(4) = g(1,1)*g(1,2) + g(2,1)*g(2,2) + g(3,1)*g(3,2)  !U^2(1,2)
   stran(5) = g(1,1)*g(1,3) + g(2,1)*g(2,3) + g(3,1)*g(3,3)  !U^2(1,3)
   stran(6) = g(1,2)*g(1,3) + g(2,2)*g(2,3) + g(3,2)*g(3,3)  !U^2(2,3)

   eflag = 1                             !initializes to no spectral decomposition
   elen = e%dvol**(1d0/3d0)      ! ELEMENT FRACTURE LENGHT (initial)
   !************************** STANDARD SOLID *****************************************
   IF( nlayr == 0 )THEN            !standard solid section (not layered)
     ! ------------------------------------------------------------------------------------
     SELECT CASE (mtype)
     ! ------------------------------------------------------------------------------------
     CASE (1)                        !Isotropic Elastic - Anisotropic Plastic  ( mtype == 1)
       IF( small )THEN               !Green Lagrange strains (twice)
         delta = (SUM(stran(1:3))-3d0)/6d0       !volumetric strain /3
         stran(1:3) = stran(1:3)-delta*2d0-1d0   !shear strains (twice)
       ELSE   ! Hencky (logarithmic) strains (volumetric and deviatoric parts)
         !computes Hencky (logarithmic) strains (volumetric and deviatoric parts)
         ! stran = twice the deviatoric strains    Delta = Ln(J)/3
         CALL lgst18(stran,r,lb,'RESV04',ierr,eflag,delta)
         IF( ierr == 1 )THEN
           WRITE(55,"(8i5,6e12.4)",ERR=9999) e%lnods(1:8),stran(1:6)
           WRITE(55,"(3e15.4)",ERR=9999) x
           CALL runen3('RESV16: negative eigenvalues of U^2')
         END IF
       END IF
       !elastic (trial) strains
       IF( .NOT.elast )THEN
         strpl(1:6) = varin(1:6)              !previous (twice) plastic strains
         stran(1:6) = stran(1:6) - strpl(1:6)   !trial Elastic log strains
       END IF
       str = gm*stran                         !Trial elastic shear stresses
       press = km*delta
       IF( plast ) THEN              !Consider plasticity
         efpst = varin(7)        !effect plastic strain
         CALL corr18(str(1),str(2),str(3),str(4),str(5),str(6), &
                     efpst,gm,props,fr,ierr,dstpl,pflag,is,numpt,val)
         IF(ierr == 1) RETURN              !no convergence in plasticity
         IF( pflag )THEN                      !if plastic flow
           varin(1:6) = varin(1:6) + dstpl(1:6)  !total plastic shear strains
           varin(7) = varin(7) + efpst           !Total effect. plastic strain
           flag2 = .TRUE.
         END IF
       END IF
       str(1:3) = str(1:3) + press
     ! ------------------------------------------------------------------------------------
     CASE (5)                        !Orthotropic Elastic - No Plasticity      ( mtype == 5)
       stran(1:3) = (stran(1:3)-1d0)/2d0   !2nd Piola Kirchhoff stress
       str(1) = eprop(1)*stran(1)+eprop(2)*stran(2)+eprop(3)*stran(3)
       str(2) = eprop(2)*stran(1)+eprop(4)*stran(2)+eprop(5)*stran(3)
       str(3) = eprop(3)*stran(1)+eprop(5)*stran(2)+eprop(6)*stran(3)
       str(4) = eprop(7)*stran(4)
       str(5) = eprop(8)*stran(5)
       str(6) = eprop(9)*stran(6)
     ! ------------------------------------------------------------------------------------
     CASE (6)                        !Hyperelastic isotropic rubber            ( mtype == 6)
       CALL rubber3d(eprop,stran,mat%matdef(8),stre=str)
     ! ------------------------------------------------------------------------------------
     END SELECT

   !************************** LAYERED SOLID *************************************************
   ELSE  !nlayr >= 1      layered solid (at integration point level)

     ! check if laminate integration is necessary
     IF (elast) THEN
       ttti = .FALSE.
     ELSE
       aux = ABS(stran(1)-1d0)+ABS(stran(2)-1d0)+ABS(stran(3)-1d0)+ &
             ABS(stran(4))+ABS(stran(5))+ABS(stran(6))
       ttti = aux > minstr   !laminate integration?
     END IF
     pflag = .FALSE.                    !initializes flag

     IF( ttti )THEN
       IF( small )THEN               !Green Lagrange strains (twice)
         stran(1:3) = (stran(1:3)-1d0)*0.5   !shear strains (twice)
       ELSE   ! Hencky (logarithmic) strains (volumetric and deviatoric parts)
         !computes Hencky (logarithmic) strains (volumetric and deviatoric parts)
         ! stran = twice the deviatoric strains    Delta = Ln(J)/3
         CALL lgst18(stran,r,lb,'RESV16',ierr,eflag,delta)
         IF( ierr == 1 )THEN
           WRITE(55,"(6i5,6e12.4)",ERR=9999) e%lnods(1:6),stran(1:6)
           WRITE(55,"(3e15.4)",ERR=9999) x
           CALL runen3('RESV16: negative eigenvalues of U^2')
         END IF
         stran(1:3) = (stran(1:3)/2d0+delta)
       END IF
       str = 0d0                         !Initialization
       !Trought laminate integration loop (rule of mixtures)
       DO il=1,nlayr                        !for each layer
         !Rotate strains from element basis to layer basis (if necessary)
         IF( rr(1,il) < 0.9999999999 )THEN  ! Angle is not zero
           stral(1) = rr(1,il)*stran(1)+ rr(2,il)*stran(2)+ rr(3,il)*stran(4) !e11
           stral(2) = rr(2,il)*stran(1)+ rr(1,il)*stran(2)- rr(3,il)*stran(4) !e22
           stral(3) = stran(3)                                                !e33
           stral(4) =-rr(4,il)*stran(1)+ rr(4,il)*stran(2)+ rr(5,il)*stran(4) !e12
           stral(5) = stran(5)                                                !e13
           stral(6) = stran(6)                                                !e23
         ELSE
           stral(1:6) = stran(1:6)
         END IF

         IF (il==nlayr) THEN
           ilv=lvari(il)
           jlv=nvar
         ELSE
           ilv=lvari(il)
           jlv=lvari(il+1)-1
         END IF

         mate  = sec%iprop(2+2*il)
         newmt = mate /= oldm
         oldm = mate
         mat => matsl(1,il)%p        ! pointer to layer material

         SELECT CASE (mat%mtype)
         ! ------------------------------------------------------------------------------------
         CASE (5)                        !Orthotropic Elastic - No Plasticity      ( mtype == 5)
           strel(1) = eprop(1)*stral(1)+eprop(2)*stral(2)+eprop(3)*stral(3)
           strel(2) = eprop(2)*stral(1)+eprop(4)*stral(2)+eprop(5)*stral(3)
           strel(3) = eprop(3)*stral(1)+eprop(5)*stral(2)+eprop(6)*stral(3)
           strel(4) = eprop(7)*stral(4)
           strel(5) = eprop(8)*stral(5)
           strel(6) = eprop(9)*stral(6)
         ! ------------------------------------------------------------------------------------
         END SELECT

         !Rotate stresses from layer basis to element basis (if necessary)
         IF( rr(1,il) < 0.9999999999 )THEN  ! Angle is not zero
           staux(1) = strel(1)*rr(1,il)+ strel(2)*rr(2,il)- strel(4)*rr(4,il)
           staux(2) = strel(1)*rr(2,il)+ strel(2)*rr(1,il)+ strel(4)*rr(4,il)
           staux(3) = strel(3)
           staux(4) = strel(1)*rr(3,il)- strel(2)*rr(3,il)+ strel(4)*rr(5,il)
           staux(5:6) = strel(5:6)
         ELSE
           staux(1:6) = strel(1:6)
         END IF

         !***   compute laminate stress on integration point level
         str(1:6) = str(1:6) + staux(1:6)*thickl(il)/thick   ! laminate integrated stresses
       END DO !l (loop to next layer)

     ELSE !elastic response of integrated laminated solid
       stran(1:3) = (stran(1:3) - 1d0 )*0.5   !2nd Piola Kirchhoff stress
       str(1:6) = MATMUL(RESHAPE((dm),(/6,6/)),stran(1:6))
     END IF

   END IF

   IF(.NOT.ASSOCIATED(e%gausv))THEN
     IF( flag2 )THEN
       e%gausv => varin
       !IF( egst < rmin)THEN
       !  rmin = egst
       !  WRITE(55,"(' minimum value changed', e12.4, 'elem',i5)")rmin,ielem
       !  WRITE(55,"(8e13.4)")varin
       !END IF
     ELSE
       DEALLOCATE(varin)
     END IF
   END IF
   e%stint(:) = str ! keep average Kirchhoff stress for output


   IF( eflag == 0)THEN
     ! Computes Hencky stress on the natural Frame T(L) = R^T T R
     sig(1) = r(1,1)*str(1)*r(1,1)+r(2,1)*str(2)*r(2,1)+r(3,1)*str(3)*r(3,1)+ &
         2d0*(r(1,1)*str(4)*r(2,1)+r(1,1)*str(5)*r(3,1)+r(2,1)*str(6)*r(3,1))
     sig(2) = r(1,2)*str(1)*r(1,2)+r(2,2)*str(2)*r(2,2)+r(3,2)*str(3)*r(3,2)+ &
         2d0*(r(1,2)*str(4)*r(2,2)+r(1,2)*str(5)*r(3,2)+r(2,2)*str(6)*r(3,2))
     sig(3) = r(1,3)*str(1)*r(1,3)+r(2,3)*str(2)*r(2,3)+r(3,3)*str(3)*r(3,3)+ &
         2d0*(r(1,3)*str(4)*r(2,3)+r(1,3)*str(5)*r(3,3)+r(2,3)*str(6)*r(3,3))
     sig(4) = r(1,1)*str(1)*r(1,2)+r(2,1)*str(2)*r(2,2)+r(3,1)*str(3)*r(3,2)+ &
              r(1,1)*str(4)*r(2,2)+r(1,1)*str(5)*r(3,2)+r(2,1)*str(6)*r(3,2)+ &
              r(2,1)*str(4)*r(1,2)+r(3,1)*str(5)*r(1,2)+r(3,1)*str(6)*r(2,2)
     sig(5) = r(1,1)*str(1)*r(1,3)+r(2,1)*str(2)*r(2,3)+r(3,1)*str(3)*r(3,3)+ &
              r(1,1)*str(4)*r(2,3)+r(1,1)*str(5)*r(3,3)+r(2,1)*str(6)*r(3,3)+ &
              r(2,1)*str(4)*r(1,3)+r(3,1)*str(5)*r(1,3)+r(3,1)*str(6)*r(2,3)
     sig(6) = r(1,2)*str(1)*r(1,3)+r(2,2)*str(2)*r(2,3)+r(3,2)*str(3)*r(3,3)+ &
              r(1,2)*str(4)*r(2,3)+r(1,2)*str(5)*r(3,3)+r(2,2)*str(6)*r(3,3)+ &
              r(2,2)*str(4)*r(1,3)+r(3,2)*str(5)*r(1,3)+r(3,2)*str(6)*r(2,3)
     ! Computes 2nd P-K stress on the natural Frame S(L) = F(T(L)
     str(1) = sig(1)/lb(1)**2
     str(2) = sig(2)/lb(2)**2
     str(3) = sig(3)/lb(3)**2
     IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
       str(4) = sig(4)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
     ELSE                                !lb(1) = lb(2)
       str(4) = sig(4)/lb(1)/lb(2)
     END IF
     IF( ABS(lb(1)-lb(3)) > 1.d-6)THEN   !lb(1) /= lb(3)
       str(5) = sig(5)*2d0*LOG(lb(1)/lb(3))/(lb(1)**2-lb(3)**2)
     ELSE                                !lb(1) = lb(3)
       str(5) = sig(5)/lb(1)/lb(3)
     END IF
     IF( ABS(lb(2)-lb(3)) > 1.d-6)THEN   !lb(2) /= lb(3)
       str(6) = sig(6)*2d0*LOG(lb(2)/lb(3))/(lb(2)**2-lb(3)**2)
     ELSE                                !lb(2) = lb(3)
       str(6) = sig(6)/lb(2)/lb(3)
     END IF
     ! Computes 2nd P-K on the Lagrangian Frame S = R S(L) R^T
     sig(1) = r(1,1)*str(1)*r(1,1)+r(1,2)*str(2)*r(1,2)+r(1,3)*str(3)*r(1,3)+ &
         2d0*(r(1,1)*str(4)*r(1,2)+r(1,1)*str(5)*r(1,3)+r(1,2)*str(6)*r(1,3))
     sig(2) = r(2,1)*str(1)*r(2,1)+r(2,2)*str(2)*r(2,2)+r(2,3)*str(3)*r(2,3)+ &
         2d0*(r(2,1)*str(4)*r(2,2)+r(2,1)*str(5)*r(2,3)+r(2,2)*str(6)*r(2,3))
     sig(3) = r(3,1)*str(1)*r(3,1)+r(3,2)*str(2)*r(3,2)+r(3,3)*str(3)*r(3,3)+ &
         2d0*(r(3,1)*str(4)*r(3,2)+r(3,1)*str(5)*r(3,3)+r(3,2)*str(6)*r(3,3))
     sig(4) = r(1,1)*str(1)*r(2,1)+r(1,2)*str(2)*r(2,2)+r(1,3)*str(3)*r(2,3)+ &
              r(1,1)*str(4)*r(2,2)+r(1,1)*str(5)*r(2,3)+r(1,2)*str(6)*r(2,3)+ &
              r(1,2)*str(4)*r(2,1)+r(1,3)*str(5)*r(2,1)+r(1,3)*str(6)*r(2,2)
     sig(5) = r(1,1)*str(1)*r(3,1)+r(1,2)*str(2)*r(3,2)+r(1,3)*str(3)*r(3,3)+ &
              r(1,1)*str(4)*r(3,2)+r(1,1)*str(5)*r(3,3)+r(1,2)*str(6)*r(3,3)+ &
              r(1,2)*str(4)*r(3,1)+r(1,3)*str(5)*r(3,1)+r(1,3)*str(6)*r(3,2)
     sig(6) = r(2,1)*str(1)*r(3,1)+r(2,2)*str(2)*r(3,2)+r(2,3)*str(3)*r(3,3)+ &
              r(2,1)*str(4)*r(3,2)+r(2,1)*str(5)*r(3,3)+r(2,2)*str(6)*r(3,3)+ &
              r(2,2)*str(4)*r(3,1)+r(2,3)*str(5)*r(3,1)+r(2,3)*str(6)*r(3,2)
   ELSE
     sig = str
   END IF

   !***   compute Int(B**t*sigma) on element level

   sig = sig*e%dvol     !Original dvol used (TLF)

   CALL proma3(rl(1,1),bmat(1,1),sig(1),12,6,.TRUE.)  !nodal contributions

   DO i=1,nnode             !for each node in the element
     n = e%lnods(i)     !node
     resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
   END DO

   e => e%next          !point to next node in the list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE resv04
