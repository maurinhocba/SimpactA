 SUBROUTINE resv15( heade, iffix, coora, resid, logst, ierr, ttime, &
                    bottom, top, coorb, coort, ifact, nrf, rhead, stint, &
                    shear,shears,factors,ninv)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT TL BBST++
 !
 !-------------------------------------------------------------------
 !USE npoi_db, ONLY : tempe, itemp
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN) :: logst       !flags for small or large strain analysis
 LOGICAL, INTENT(IN) :: bottom,top  !flags for top and bottom surface computation
 INTEGER (kind=4), INTENT(IN) :: iffix(:), &  !nodal fixities
                                 nrf          !number of reinforced sides
 INTEGER (kind=4), INTENT(OUT) :: ierr,     & !flag for error detection
                               ifact(:)       !weigths for surface computation
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                                  ttime       !present total time
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:),& !internal nodal forces
                                 coorb(:,:),& !bottom surface coordinates
                                 coort(:,:)   !top surface coordinates
 TYPE (ele15), POINTER  :: heade              !pointer to first element
 TYPE (reinf), POINTER  :: rhead              !pointer to first element
 REAL(kind=8), POINTER :: stint(:,:)          !moments and shear
 INTEGER(kind=4), INTENT (IN) :: shear        !/=0 if shear are to be computed
 REAL(kind=8), POINTER :: shears(:,:),factors(:)
 INTEGER(kind=4), POINTER :: ninv(:)

 ! local variables

 TYPE (ele15), POINTER :: e    !pointer to an element data
 TYPE (reinf), POINTER  :: rf    !pointer to branching side

 LOGICAL :: isidf(3)  ! TRUE if side is clamped or symmetry line

 INTEGER (kind=4) sides(3)  ! 1 if side element exist

 INTEGER (kind=4)  :: iel         !element order

 REAL (kind=8) stra1(6),  & !first and second fundamental forms
               bmat(3,6,6),& !membrane and bending B matrix
               x(3,6),    & !nodal coordinates of the patch
               t(3,15),   & !Local triad + local in-plane vector at side elms
               rl(3,6)      !nodal internal equivalent forces

 REAL (kind=8) ::  mn(3,3),   & !matrix to compute shear at boundaries
                   n1,n2,ls,  & !normal to the side and length
                   n11,n12,n22  !factors


 LOGICAL :: large,   &! TRUE if log strains must be used
            rf2,rfn,     &! TRUE if normal distribution is to be done
            ttti      ! TRUE if Trougth The Thickness Integration is necessary

 LOGICAL :: pflag(7) ! TRUE if plastic flow in the step

 LOGICAL :: newmt,   &! TRUE if material constant computation necessary
            natst,   &! TRUE for large strain analysis
            elast,   &! TRUE if material is strictly elastic
            shell,   &! TRUE if bending included
            plast,   &! TRUE if plasticity is to be considered
            coupled   ! TRUE if membrane-bending coupling for laminates

 INTEGER (kind=4) ielem, & !element number
                  i,j,k,n,l,& !different indexes
                  error    !error flag

 INTEGER (kind=4) isec, isrf,  & !associated material
                  nl,nlayr, & !number of layers
                  mtype, & !associated material type
                  secty, & !section type
                  numpt, & !number of points in curve
                  osec, osrf, & !associated material of previous element
                  nvar     !number of internal variables per layer

 REAL (kind=8) stres(3),  & !stresses (different measures)
               sigma(3),  & !stresses (different measures)
               stran(3),  & !C=U^2  also Log strains
               stine(10), & !t.t.t integrated stresses (forces and moments)
               r1,r2,     & !eigevenctor components in local system
               lb(3),     & !eigenvalues
               thnew,     & !present thickness
               zk,        & !distance to mid surface
               aux,area1, & !auxiliar value
               !! t0,t1,j0,  & !thermical dilatation coeff
               a,b,rfp(3)   !reinforced B matrix
 REAL (kind=8) thick,     & !thickness (original)
               c(4),gh,   & !Elastic constitutive matrix for plane stress
               dm(21),    & !Elastic integrated constitutive matrix
               alpha,     & !thermical dilatation coeff
               propi(15), & !Plastic Material properties
               chi(12),chid(6),chib(4),   & !Hill 48 coefficients
               deatht,    & !end time for plasticity
               minstr,    & !minimum strain to integrate trougth the thickness
               min_tr,    & !Minimum Allowable Thickness ratio
               max_tr       !Maximum Allowable Thickness ratio
 ! Gauss points throught the thickness
 REAL (kind=8) :: thf(7),wei(7)

 !REAL (kind=8), ALLOCATABLE :: cm(:,:), prop(:,:), volfr(:), rr(:,:)
 REAL (kind=8), POINTER :: val(:,:)
 REAL (kind=8), POINTER :: varin(:,:)                      !internal variables

 CHARACTER(len=1 ) letra  !for console input

 TYPE (section), POINTER :: sec, srf  !pointer to a section data
 TYPE (mater), POINTER :: mat    !pointer to a section data


 INTERFACE
   INCLUDE 'gentbs.h'
   INCLUDE 'stre14.h'
   INCLUDE 'rubberps.h'
   !INCLUDE 'ud_shell_2.h'
 END INTERFACE


 IF( nrf > 0 )THEN    !if reinforcement exist
   rf => rhead          !point to first
   DO i=1,nrf              !loop
     rf%stint = 0d0          !null forces and moments
     rf => rf%next           !poin to next
   END DO
 END IF

 !***  Computes nodal internal forces, integrating trought the thickness

 osec = -1                           !no previous material
 osrf = -1                           !no previous material
 error = 0
 x = 0d0                             !initializes
 e => heade                          !point to first element in the list
 iel = 0
 DO                                  !loop over all the elements
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   ielem = e%numel                   !element number
   iel = iel + 1
   isec = e%matno                    !associated material
   newmt = osec /= isec              !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     sec => psecs(isec)%p             !point to section
     nlayr = sec%iprop(1)             !number of layers
     nl = nlayr
     nvar  = sec%iprop(2)             !number of internal variables per layer
     secty = sec%secty                !section type
     thick = sec%rprop(1)             !original thickness
     min_tr= sec%rprop(2)             !Minimum Thickness ratio
     max_tr= sec%rprop(3)             !Maximum Thickness ratio
     minstr= sec%rprop(4)             !strain threshold to use TTTI
     osec = isec                      !keep present section
     IF( secty == 12 )THEN            !standard solid section
       mat => sec%mtbas                 !point to associated material
       shell = nlayr > 1                !membrane or shell
       IF( .NOT.shell )rl(:,4:6) = 0d0
       mtype = mat%mtype                !type of base material
       natst = logst .OR. mtype == 6    !use log strains
       elast = mat%matdef(3) == 1       !elastic
       CALL gaussq(nlayr,thf(1),wei(1)) !integration points through the thickness
       thf(1:nlayr) = thf(1:nlayr)/2d0  !positions
       wei(1:nlayr) = wei(1:nlayr)/2d0  !weights
       alpha = mat%prope(6)             !Thermical Dilatation Coeff
       !IF( .NOT.itemp ) alpha = 0d0
       dm(1:8) = sec%rprop(5:12)             !integrated elasticity matrix

       IF( mtype == 1)THEN              !standard isotropic material
         IF( elast )THEN
           plast = .FALSE.
         ELSE   !for elasto-plastic mats
           ! e1, nu1, uniaxial, efren, consn, r, exponent m, hill 79
           propi(1:4) = mat%propp(1:4)       ! isotropic hardening parameters
           propi(5) = REAL( mat%matdef(4),8) ! isotropic hardening model
           chi    = mat%propp(16:27)         ! hill coefficients
           deatht = mat%propp(5)             !end of plasticity
           IF (mat%matdef(4) == 5 )THEN
             val => mat%chead%val
             numpt = mat%chead%np
           ELSE
             NULLIFY (val)
             numpt = 0
           END IF
           plast = deatht > ttime .AND. propi(1) > 0  !consider plasticity ?
           IF(plast)THEN
             SELECT CASE (mat%matdef(3))
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
           END IF
         END IF
         c(1:4) = mat%prope(7:10)          ! plane stress elasticity matrix
         gh = 5d0/6d0*mat%prope(3)*sec%rprop(1)
       ELSE IF( mtype == 5)THEN            ! orhthotropic material
         c(1:4) = mat%prope(16:19)         ! plane stress elasticity matrix
         elast = .TRUE.
         plast = .FALSE.                        !consider plasticity ?
         !IF( .NOT.elast )THEN              ! if elasto-plastic
         !  propi(1:15) = mat%propp(17:29)    !orthotropic hardening parameters
         !  deatht = mat%propp(5)             !end of plasticity
         !  plast = deatht > ttime .AND. propi(1) > 0  !consider plasticity ?
         !END IF
         newmt = .FALSE.                        !same material than previous ?
       ELSE IF( mtype == 6)THEN            !hyperelastic isotropic
         chi(1:12) = mat%prope(7:18)
         elast = .TRUE.
         plast = .FALSE.                   !consider plasticity ?
       !ELSE IF( mtype == 30)THEN             !user defined material
       !  plast = .TRUE.   !default
       END IF

     ELSE  !secty == 15  composite laminae
       elast = .TRUE. !elast = sec%iprop(3) == 0               !elastic problem
       plast = .FALSE.!plast = .NOT.elast                      !consider plasticity ?
       shell = .TRUE.
       dm = sec%rprop(6:26)                 !linear elastic integrated matrix
       coupled = ANY(dm(15:21) /= 0d0 )     !coupling terms
       !IF( plast )THEN  !plastic lamina
       !  IF( ALLOCATED ( cm ) )DEALLOCATE( volfr,cm,prop,rr )
       !  ALLOCATE( volfr(nlayr),cm(4,nlayr),prop(17,nlayr),rr(5,nlayr))
       !  i = 27                                !pointer to volume fractions
       !  CALL vecasi(nlayr,sec%rprop(i),volfr)        !volume fractions
       !  i = i+2*nlayr                         !pointer to C-matrices
       !  CALL vecasi(nlayr*4,sec%rprop(i),cm(1,1))    !elasticity matrix
       !  i = i+nlayr*4                         !pointer to R-matrices
       !  CALL vecasi(nlayr*5,sec%rprop(i),rr(1,1))    !rotation matrix
       !  i = i+nlayr*5                         !pointer to plastic properties
       !  CALL vecasi(nlayr*17,sec%rprop(i),prop(1,1)) !plastic properties
       !  oldm = -1                             !recompute constants
       !END IF
       newmt = .FALSE.                        !same material than previous ?
     END IF
   END IF

   !*** computes actual configuration

   isidf = .FALSE.                  !initializes
   x(:,1:3) = coora(:,e%lnods(1:3)) ! get coordinates of the main triangle
   DO i =1,3                        !for each node in the patch
     n = e%lnods(i+3)               ! node number
     IF(n > 0) THEN               ! adjacent element exist ?
       sides(i) = 1                 ! side element exist
       x(:,i+3) = coora(:,n)        ! get coordinates of the side node
     ELSE ! IF( n == 0 )THEN        ! boundary node
       sides(i) = 0                 ! side does NOT exists
       j = e%lnods(hh(3,i))         ! first node of the side
       k = e%lnods(hh(2,i))         ! second node of the side
       isidf(i) = iffix(j) == 1 .AND. iffix(k) == 1 !free or clamped ?
     END IF
   END DO
   ! compute first and second fundamental forms
   CALL stra15(e%a,e%b,e%c,e%cd,e%ci,e%a0,e%area1,t,e%lb,x,sides,isidf, &
               e%mems,.FALSE.,stra1,e%gamma)

   !check maximum and miminum thickness ratio
   ! uses
   !
   !   INTEGER (kind=4), INTENT(OUT) :: ierr    !flag for error detection
   !   REAL (kind=8) min_tr,    & !Minimum Allowable Thickness ratio
   !                 max_tr       !Maximum Allowable Thickness ratio
   !
   !   CHARACTER(len=1 ) letra  !for console input
   !
   !   TYPE (ele14), POINTER :: e    !pointer to an element data
   !   TYPE (section), POINTER :: sec    !pointer to a section data

   IF( e%lb < min_tr .OR. e%lb > max_tr)THEN      !check
     WRITE(*,"(' thickness ratio limit exceeded, TH_RATIO: ',f5.3, &
           &  /,' Do you want to continue [Y/N]')")e%lb
     DO                            !pause until answer
       READ(*,'(a1)')letra         !read the answer
       SELECT CASE (letra)         !process answer
       CASE ( 'y', 'Y')         !continue
         IF( e%lb < min_tr )THEN    !minimum exceeded
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

   e%stra1 = stra1         ! keep Mid Surface metric tensors
   ! consider thermal strains
   !IF( alpha > 0 )THEN     ! compute
   !  t0 = 1d0+alpha*SUM(tempe(e%lnods(1:3),1))
   !  e%lb = e%lb*t0        ! modify thickness ratio
   !END IF

   thnew = thick*e%lb                  !thickness at new configuration

!   IF( alpha > 0 )t1 = alpha*SUM(tempe(e%lnods(1:3),2))*2d0/thnew

   IF( bottom .OR. top) &      !bottom or top surfaces required
      CALL gentbs(isidf,e%lnods(1:3),bottom,top,coorb,coort,ifact,t(:,3),thnew,e%c(:,1:3,2),3)

   area1 = e%area1
   stine = 0d0                         !Integrated forces, moments, shears & VMises stresses

   ! check if trought the thickness integration is necessary

   IF( ASSOCIATED(e%gausv) )THEN
     ttti = .TRUE.
     varin => e%gausv
   ELSE
     aux = ABS(stra1(1)-1d0)+ABS(stra1(2)-1d0)+ABS(stra1(3))+ &
          (ABS(stra1(4))+ABS(stra1(5))+ABS(stra1(6)))*thick/2d0
     ttti = aux > minstr
     IF( ttti )THEN
       ALLOCATE( varin(nvar,nlayr) )
       varin = 0d0
     END IF
   END IF

   IF( secty == 12 )THEN              !for standard solid section
     !IF( ttti )THEN
     IF( .FALSE. )THEN
       !Trought the thickness (t.t.t.) integration loop
       pflag = .FALSE.
       DO l=1,nlayr                        !for each layer
         large = natst                     !use log strains
         IF( plast ) large = large .OR. ASSOCIATED(e%gausv)  !use log strains
         zk = thf(l)*thnew                     !Z coordinate
         stran = stra1(1:3)+stra1(4:6)*zk      !layer U^2
         !IF( alpha > 0 )THEN                 !consider thermal strains
         !  j0 = (t0 + zk*t1)**(2d0/3d0)      !modifies jacobian due to temperature
         !  stran = stran/j0                  !modifies C
         !END IF
         IF( large )THEN
           CALL lgst14(stran,r1,r2,lb,'RESVPL',error)    !Hencky (logarithmic) strains
           IF( error == 1 )THEN
             WRITE(55,"(i7,3i6,7e12.4)",ERR=9999) e%numel,e%lnods(1:3),stra1(1:6),zk
             WRITE(55,"(3e15.4)",ERR=9999) x
             ierr = 1
             error = 0
             !CALL runend('RESVPL: negative eigenvalues of U^2')
           END IF
         ELSE
           stran(1) = (stran(1) - 1d0 )/2d0
           stran(2) = (stran(2) - 1d0 )/2d0
         END IF
         IF( mtype == 1 )THEN     !one component material
           CALL stre14(stran,stres,c,propi,chib,chid,varin(1:4,l),ierr,3, &
                       plast,elast,val,numpt,aux,mat%matdef(3),pflag(l))
         ELSE IF( mtype == 5  )THEN
           !elastic (trial) stresses
           !aux = SQRT(2d0*area1)   ! characteristic length
           !IF( ASSOCIATED(e%gausv)) THEN
           !  IF(e%gausv(4,l) > 0d0 )stran(1:3) = stran(1:3) - e%gausv(1:3,l)
           !END IF
           stres(1)= c(1)*stran(1)+c(2)*stran(2)
           stres(2)= c(2)*stran(1)+c(3)*stran(2)
           stres(3)= c(4)*stran(3)

         ELSE IF ( mtype == 6 )THEN ! Hyperelastic
           lb(3) = 1d0/lb(1)/lb(2)
           CALL rubberps(chi,lb,mat%matdef(8),stres,r1=r1,r2=r2)

         END IF
         IF( l == 1     ) stine( 9) = aux
         IF( l == nlayr ) stine(10) = aux
         IF( natst .AND. mtype /= 6 )THEN
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
       stine = stine*thick                  !Original thick because it's TLF
       IF(.NOT.ASSOCIATED(e%gausv))THEN
         IF(ANY(pflag(1:nl)))THEN
           e%gausv => varin
         ELSE
           DEALLOCATE(varin)
         END IF
       END IF
     ELSE          !no t.t.t.i.
       stra1(1) = (stra1(1) - 1d0 )/2d0
       stra1(2) = (stra1(2) - 1d0 )/2d0
       stra1(4:5) = stra1(4:5)/2d0
       stine(1) = dm(1)*stra1(1) + dm(2)*stra1(2)
       stine(2) = dm(2)*stra1(1) + dm(3)*stra1(2)
       stine(3) = dm(4)*stra1(3)
       stine(4) = dm(5)*stra1(4) + dm(6)*stra1(5)
       stine(5) = dm(6)*stra1(4) + dm(7)*stra1(5)
       stine(6) = dm(8)*stra1(6)
     END IF

   ELSE   !secty = 15
     !IF( ttti )THEN
     !  stine = 0d0                         !Integrated forces and moments
     !  !Trought the thickness (t.t.t.) integration loop
     !  zk = -(thick+volfr(1))/2d0*e%lb
     !  DO l=1,nlayr                        !for each layer
     !    zk = zk+volfr(l)*e%lb                 !Z coordinate
     !    stran = stra1(1:3)+stra1(4:6)*zk      !layer U^2
     !    !IF( alpha > 0 )THEN                 !consider thermal strains
     !    !  j0 = (t0 + zk*t1)**(2d0/3d0)      !modifies jacobian due to temperature
     !    !  stran = stran/j0                  !modifies C
     !    !END IF
     !    IF( natst )THEN
     !      CALL lgst14(stran,r1,r2,lb,'RESVPL',ierr)    !Hencky (logarithmic) strains
     !      IF( ierr == 1 )THEN
     !        write(55,"(3i5,7e12.4)",ERR=9999) e%lnods(1:3),stra1(1:6),zk
     !        write(55,"(3e15.4)",ERR=9999) x
     !        CALL runen3('RESVPL: negative eigenvalues of U^2')
     !      END IF
     !    ELSE
     !      stran(1) = (stran(1) - 1d0 )/2d0
     !      stran(2) = (stran(2) - 1d0 )/2d0
     !    END IF
     !    mate  = sec%iprop(2+2*l)
     !    mtype = sec%iprop(3+2*l)
     !    IF( mtype == 1 )THEN     !one component material
     !      newmt = mate /= oldm
     !      oldm = mate
     !      CALL stre14(stran,stres,cm(:,l),prop(1:5,l),prop(6:17,l),e%gausv(:,l),ierr,3, &
     !                  newmt,plast,elast,val,numpt,aux,mat%matdef(3))
     !    ELSE IF( mtype == 5  )THEN
     !      !elastic (trial) stresses
     !      stres(1)= cm(1,l)*stran(1)+cm(2,l)*stran(2)
     !      stres(2)= cm(2,l)*stran(1)+cm(3,l)*stran(2)
     !      stres(3)= cm(4,l)*stran(3)
     !      IF( plast )THEN
     !        ! ROTATE STRAINS BEFORE
     !        aux = SQRT(2d0*area1)   ! characteristic length
     !        CALL corr05(stres(1),stres(2),stres(3),e%gausv(4,l),cm(:,l),prop(:,l),ierr, &
     !                    e%gausv(1:3,l),aux)
     !        ! ROTATE STRESSES AFTER
     !      END IF
     !
     !    END IF
     !
     !    IF( natst )THEN
     !     ! Computes Hencky stress on the natural Frame
     !      sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2+2d0*stres(3)*r1*r2
     !      sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1-2d0*stres(3)*r1*r2
     !      sigma(3) =(stres(2)-stres(1))*r1*r2+stres(3)*(r1*r1-r2*r2)
     !     ! Computes 2nd P-K stress on the natural Frame
     !      stres(1) = sigma(1)/lb(1)**2
     !      stres(2) = sigma(2)/lb(2)**2
     !      IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
     !        stres(3) = sigma(3)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
     !      ELSE                                !lb(1) = lb(2)
     !        stres(3) = sigma(3)/lb(1)/lb(2)
     !      END IF
     !     ! Computes 2nd P-K on the Lagrangian Frame
     !      sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2-2d0*stres(3)*r1*r2
     !      sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1+2d0*stres(3)*r1*r2
     !      sigma(3) =(stres(1)-stres(2))*r1*r2+stres(3)*(r1*r1-r2*r2)
     !    ELSE
     !      sigma = stres  !small strain
     !    END IF
     !
     !    !***   compute Int(B**t*sigma) on element level
     !
     !    stine(1:3) = stine(1:3)+sigma*volfr(l)   !t.t.t. integrated forces
     !    stine(4:6) = stine(4:6)+sigma*zk*volfr(l)  !integrated moments
     !  END DO
     !
     !
     !ELSE
       stra1(1) = (stra1(1) - 1d0 )/2d0
       stra1(2) = (stra1(2) - 1d0 )/2d0
       stra1(4:5) = stra1(4:5)/2d0
       stine(1) = dm(1)*stra1(1) + dm(2)*stra1(2) + dm(3)*stra1(3)
       stine(2) = dm(2)*stra1(1) + dm(4)*stra1(2) + dm(5)*stra1(3)
       stine(3) = dm(3)*stra1(1) + dm(5)*stra1(2) + dm(6)*stra1(3)
       stine(4) = dm(7)*stra1(4) + dm(8)*stra1(5) + dm(9)*stra1(6)
       stine(5) = dm(8)*stra1(4) + dm(10)*stra1(5) + dm(11)*stra1(6)
       stine(6) = dm(9)*stra1(4) + dm(11)*stra1(5) + dm(12)*stra1(6)
       IF( coupled )THEN
         stine(1) = stine(1) +dm(15)*stra1(4) + dm(14)*stra1(5) + dm(15)*stra1(6)
         stine(2) = stine(2) +dm(16)*stra1(4) + dm(17)*stra1(5) + dm(18)*stra1(6)
         stine(3) = stine(3) +dm(19)*stra1(4) + dm(20)*stra1(5) + dm(21)*stra1(6)
         stine(4) = stine(4) +dm(15)*stra1(1) + dm(16)*stra1(2) + dm(19)*stra1(3)
         stine(5) = stine(5) +dm(14)*stra1(1) + dm(17)*stra1(2) + dm(20)*stra1(3)
         stine(6) = stine(6) +dm(15)*stra1(1) + dm(18)*stra1(2) + dm(21)*stra1(3)
       END IF
     !END IF
   END IF

   IF( shear /= 0 )THEN

     DO i=1,3
       n = e%lside(i)                 ! element side
       IF(n > 0) THEN                 ! element exist ?
         stres(1:3) = stint(4:6,n)*e%ci(i) + stint(4:6,iel)*(1-e%ci(i))
       ELSE IF( n == 0 ) THEN !boundary side
         ls  = SQRT(e%b(i)**2+e%a(i)**2)     !side pseudo length
         n1 =  e%b(i)/ls                     !normalizes side to compute normal
         n2 = -e%a(i)/ls
         n11 = n1*n1                         !auxiliar factors
         n12 = n1*n2
         n22 = n2*n2
         IF (isidf(i) )THEN !clamped/symmetry
           mn = RESHAPE ((/   n11*n11+n22*n22,       2d0*n11*n22, n12*(n11-n22),    &
                                  2d0*n11*n22,   n11*n11+n22*n22, n12*(n22-n11),    &
                            2d0*n12*(n11-n22), 2d0*n12*(n22-n11),   4d0*n11*n22 /), &
                        (/3,3/))
         ELSE               !simple supported
           mn = RESHAPE ((/ 2d0*n11*n22+n22*n22,            -n11*n22,        -n12*n11,    &
                                       -n11*n22, 2d0*n11*n22+n11*n11,        -n12*n22,    &
                                   -2d0*n12*n11,        -2d0*n12*n22, n11*n11+n22*n22 /), &
                        (/3,3/))
         END IF
         stres(1:3) = MATMUL (mn,stint(4:6,iel))
       ELSE !branching side
         stres(1:3) = stint(4:6,iel)
       END IF
       stine(7) = stine(7) + stres(1) * e%b(i) - stres(3) * e%a(i)
       stine(8) = stine(8) + stres(3) * e%b(i) - stres(2) * e%a(i)
     END DO
     stine(7:8) = 2d0*stine(7:8)
   END IF

   stint(1:10,iel) = stine                      !store N, M, Q & vMS  for post-processing

  !       average configuration B Matrices
   stine(1:3) = stine(1:3)*area1/3d0  ! /3 to use with average matrix
   CALL bmem15(e%a,e%b,e%cd,e%a0,e%gamma,bmat(:,:,1:3),t,e%mems)
   CALL proma3(rl(1,1),bmat(1,1,1),stine(1),18,3,.TRUE.)  !membrane contributions
   IF(shell) THEN
     ! compute bending B matrix
     stine(4:6) = stine(4:6)*area1      !
     CALL bfle15(e%a,e%b,e%c,e%ci,sides,isidf,bmat(:,:,4:6),t(:,1:9))
     CALL proma3(rl(1,1),bmat(1,1,4),stine(4),18,3,.FALSE.) !bending contributions
   END IF

   IF( nrf > 0 )THEN  !if any side is reinforced
     DO i=1,3           !for each side
       IF( ASSOCIATED(e%si(i)%p))THEN  !if side is reinforced
         l = e%si(i)%p%lnods(3)        !left adjacent element + side
         n = e%si(i)%p%lnods(4)        !right adjacent element + side
         j = 1                                    !left triangle
         IF(ABS(e%si(i)%p%lnods(4))/10 == iel) j = 2
         isrf = e%si(i)%p%secno     !associated material
         newmt = osrf /= isrf       !same material than previous ?
         IF(newmt)THEN                     !if new material, get properties
           srf => psecs(isrf)%p             !point to section
           rfp(1:3) = srf%rprop((/ 5,9,8 /))  !axial,bending and torsion stiffnesses
         END IF
         !compute axial strain
         rf2 = l >= 0 .AND. n >= 0  !both adjacent elements
         IF( .NOT.rf2 )THEN
           rfn =  (l/10) == iel .OR.  (n/10) == iel
         ELSE
           rfn = .FALSE.
         END IF
         IF( rf2 .OR. rfn )THEN
           lb = x(:,hh(2,i)) - x(:,hh(3,i))           !side vector
           aux = DOT_PRODUCT(lb,lb)                   !actual length squared
           stran(1) = (aux/e%si(i)%p%l0**2 - 1d0)/2d0 !G-L axial strain
           stres(1) = rfp(1)*stran(1)         !axial force
           e%si(i)%p%stint(1) = stres(1)      !keep axial force
           IF( .NOT.rfn )stres(1) = stres(1)*e%si(i)%p%fc(j)   !use factor
           ! normal stress
           stres(1) = stres(1)/e%si(i)%p%l0                !N/L0
           rl(:,hh(3,i)) = rl(:,hh(3,i)) - lb*stres(1)     ! add to first side node
           rl(:,hh(2,i)) = rl(:,hh(2,i)) + lb*stres(1)     ! add to second side node
         ELSE
           stres(1) = 0d0
         END IF
         ! get normalized direction
         a = e%si(i)%p%dc(1)
         b = e%si(i)%p%dc(2)
         ! compute local strains (curvatures)
         n1 = (e%stra1(4)*a + e%stra1(6)*b)/2d0
         n2 = (e%stra1(5)*b + e%stra1(6)*a)/2d0
         stran(2) = a*n1+b*n2              !bending curvature
         stran(3) = b*n1-a*n2              !torsion curvature
         ! compute stresses
         stres(2) = rfp(2)*stran(2)         !bending moment
         stres(3) = rfp(3)*stran(3)         !torsion moment
         e%si(i)%p%curv(1:2,j) = stran(2:3)     !store present curvatures
         stres(2:3) = stres(2:3)*e%si(i)%p%fc(j)          !use factor
         e%si(i)%p%stint(2:3) = e%si(i)%p%stint(2:3) + stres(2:3)   !sum on reinforcement
         ! bending and torsion moment
         stres(2:3) = stres(2:3)*e%si(i)%p%l0            !M T /L0
         n1 = a*a                                        !n1*n1
         n2 = b*b                                        !n2*n2
         aux= a*b                                        !n1*n2
         ! compute auxiliar values to evaluate internal forces
         sigma(1) = n1* stres(2) +  aux*       stres(3)  !
         sigma(2) = n2* stres(2) -  aux*       stres(3)
         sigma(3) = aux*stres(2) + (n2-n1)/2d0*stres(3)
         CALL proma3(rl(1,1),bmat(1,1,4),sigma(1),18,3,.FALSE.) !reinforced contributions
       END IF
     END DO

   END IF

   DO i=1,6             !for each node in the patch
     n = e%lnods(i)     !node
     IF(n > 0) resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
   END DO

   e => e%next          !point to next node in the list
 END DO

 IF( shear < 0 ) THEN
   ! compute smoothed shears
   shears = 0d0
   e => heade
   iel = 1
   DO
     DO i=1,3
       j = ninv(e%lnods(i))
       shears(:,j) = shears(:,j) + stint(7:8,iel)/e%area1*factors(j)
     END DO
     e => e%next
     IF( .NOT.ASSOCIATED(e) )EXIT
     iel = iel+1
   END DO
   stint(7:8,:) = 0d0
   e => heade
   iel = 1
   DO
     DO i=1,3
       j = ninv(e%lnods(i))
       stint(7:8,iel) = stint(7:8,iel) + shears(:,j)/3d0
     END DO
     e => e%next
     IF( .NOT.ASSOCIATED(e) )EXIT
     iel = iel+1
   END DO
 END IF
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE resv15
