  SUBROUTINE resv16( nele0,head0, coora, resid, ierr, ttime, small, ngaus, quad, shell, sbbar)
  !-------------------------------------------------------------------
  !
  !          EVALUATION OF INTERNAL FORCES ELEMENT T-L 3-D SOLID
  !
  !-------------------------------------------------------------------
  USE gvar_db, ONLY : updiv,newiv,elastic
  USE ele18_db, ONLY : lgst18,eige18,rubber3dn
  !$ USE omp_lib
  IMPLICIT NONE

  ! dummy arguments
  INTEGER (kind=4), INTENT(OUT) :: ierr        !flag for error detection
  INTEGER (kind=4), INTENT(IN) :: ngaus,     & !number of integration points
                                  nele0        !number of elements in the set
  REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                              ttime            !present total time
  REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
  TYPE (ele16), POINTER  :: head0              !pointer to first element
  LOGICAL, INTENT(IN) :: small,   & !.TRUE. if Green strains will be used
                         quad,    & !.TRUE. if quadratic approach
                         shell,   & !.TRUE. for assumed transverse shear
                         sbbar      !.TRUE. for avergar volumetric strain

  ! local variables
  ! variables for assumed strain transverse shear
  REAL (kind=8) :: tas(3,3),      & ! der x/x_3 at Assumed Strain points
                   gas(3,3),      & ! der x/tan at Assumed Strain points
                   vas(3,2),      & ! strains at Assumed Strain Points
                   bs(3,6,3,2)      ! B matrix at Assumed strain Points
  TYPE (ele16), POINTER :: e => NULL()        !pointer to an element data
  TYPE (ele16), POINTER :: heade => NULL()    !pointer to an element data
  LOGICAL :: newmt,   & ! TRUE if material constant computation necessary
             elast,   & ! TRUE if material is elastic
             bbar,    & ! TRUE to decouple shear and volumetric terms
             plast      ! TRUE if plasticity is to be considered
  INTEGER (kind=4) :: osec,       & !associated material of previous element
                      numpt,      & !number of points in curve
                      mtype         !associated material type
  REAL (kind=8) ::  fr(9),yf(9),& !Flow rule and yield function coefficients
                    km,gm,      & !bulk modulus and shear modulus
                    eprop(12),  & !elastic Material properties
                    props(7),   & !plastic Material properties
                    ch(6)         !Hill-48 coefficients

  INTEGER (kind=4) ielem,nelem,       & !element number
                   ndofe,             & !number of DOFs of the elemnt
                   nn,                & !number of nodes
                   is,kh,             & !isotropic and kinematic hardening models
                   isec,              & !associated material
                   flag(ngaus),       & !flag for eigen-decomposition
                   i,ii,jj,j,k,n        !different indexes

  INTEGER (kind=4) nl,nlayr,  & !number of layers
                   oldm,      & !old material label
                   mate,      & !material label
                   ilv,jlv,   & !indexes for layer varin
                   il,        & !layer index
                   nvar         !total number of internal variables
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
                g(3,3),     & !deformation gradient
                r(3,3),lb(3),& !eigen-pair in local system
                delta,      & !volumetric strain/3, element volume
                deltm,      & !Average volumetric strain
                press         !first stress invariant (mean pressure)
  REAL (kind=8) :: rv(3,3,ngaus),lbv(3,ngaus),stres(6,ngaus),gv(3,3,ngaus)
  REAL (kind=8), PARAMETER :: sh(2,2) = (/ 0.7886751345948129, 0.2113248654051871,  &
                                           0.2113248654051871, 0.7886751345948129   /)
  LOGICAL :: pflag,   &! TRUE if plastic flow in the step
             ttti,    &! TRUE if Trougth The Thickness Integration is necessary
             flag2     ! TRUE if plastic at any gauss point

  REAL (kind=8)          :: elen    ! Damage length - For damage subroutine - XMG

  TYPE (section), POINTER :: sec => NULL()  !pointer to a section data
  TYPE (mater), POINTER   :: mat => NULL()  !pointer to a material data

  LOGICAL :: sides(3)
  ! allocatable arrays : MUST be replaced by any STATIC variables to avoid run-time memory increase
  REAL (kind=8), ALLOCATABLE :: bmat(:,:,:), &  !B matrix
                                x(:,:),    &  !element coordinates
                                rl(:,:),   &  !nodal internal equivalent forces
                                tq(:,:,:,:), &  !in-plane gradient at mid-side points
                                cq(:,:),   &  !average in-plane metric tensor at faces
                                bq(:,:,:,:)   !in-plane BB matrix at faces
  REAL (kind=8), POINTER :: val(:,:) => NULL()
  TYPE (pmater), POINTER :: matsl(:) => NULL()    !array of pointers to materials data
  INTEGER (kind=4), ALLOCATABLE :: lvari(:)
  REAL (kind=8), ALLOCATABLE :: thickl(:),zposl(:),rr(:,:)
  REAL (kind=8), ALLOCATABLE :: varin(:,:)                      !internal variables

  INTERFACE
    INCLUDE 'corr18.h'
    INCLUDE 'corr18k.h'
  END INTERFACE
  !Auxiliar variables for parallel code
  !!$ INTEGER(kind=4), PARAMETER :: maxel=34816,minel=32768  !8 threads
  !!$ INTEGER(kind=4), PARAMETER :: maxel=18432,minel=16384  !4 threads
  !$ INTEGER(kind=4), PARAMETER :: maxel=10240,minel= 8192    !2 threads
  !$ INTEGER(kind=4):: nloop,iloop,irest,m    !
  INTEGER(kind=4):: piter,    & !
                    lim_inf,  & ! first element for each processor
                    lim_sup,  & ! last element for each processor
                    num_chk     ! number of elements chunks

  !$ INTEGER(kind=4),PARAMETER:: MAXLOCK=2048               !number of keys
  !$ INTEGER(kind=4),SAVE :: lock1(0:MAXLOCK-1) !keys
  !!$ INTEGER(kind=OMP_lock_kind),SAVE :: lock1(0:MAXLOCK-1) !keys
  !$ LOGICAL,SAVE             :: inilock=.FALSE.           !initialization flag


  heade => head0
  nelem = nele0
  IF( quad ) THEN   !for quadratic approach
    nn = 12                          !number of nodes in the patch
  ELSE
    nn = nnode                       !element nodes
  END IF
  ndofe = 3*nn                       !number of DOFs

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

  !initialize variables
  num_chk =  1    !only one processor or thread
  lim_inf =  1    !first element
  lim_sup = nelem !last element

  !$OMP  PARALLEL                                                                           &
  !$OMP& SHARED (psecs,pmats, updiv, ngaus, coora, ttime, small, quad, shell, sbbar,        &
  !$OMP&         newiv, ierr, resid, heade, nelem, ndofe, nn, num_chk,lock1,elastic )       &
  !$OMP& PRIVATE(tas, gas, vas, bs, e, newmt, elast, bbar, plast, osec, numpt, mtype, fr, yf, kh, &
  !$OMP&         km,gm, eprop, props, ch, ielem, is, isec, flag, i,ii,jj,j,k,n,m,nl,nlayr,  &
  !$OMP&         oldm, mate, ilv,jlv, il, nvar, stral, strel, staux, thick, minstr, dm,     &
  !$OMP&         aux, str, sig, stran, dstpl, efpst, strpl, g, r,lb, delta, deltm, press,   &
  !$OMP&         rv,lbv,stres,gv, pflag, ttti, flag2, elen, sec, mat, sides,     &
  !$OMP&         bmat, x, rl, tq, cq, bq, val, matsl, lvari, thickl,zposl,rr, varin,        &
  !$OMP&         piter, lim_inf, lim_sup)

  !$ IF( omp_get_num_threads() > 1 )THEN
    !$ num_chk = (nelem+omp_get_num_threads()-1)/omp_get_num_threads()
    !$ lim_inf = MIN(nelem,1+omp_get_thread_num()*num_chk)     !Find the first element in the list
    !$ lim_sup = MIN(nelem,num_chk*(1+omp_get_thread_num()))   !Find the last element in the list
  !$ END IF

  e => heade                          !point to head of the list
  DO piter=1,lim_inf-1                !Loop over the list to find the first element in the list
    e => e%next                       !Point the first element of the list
  END DO

  !***  Computes nodal internal forces
  IF( quad ) THEN   !for quadratic approach
    ALLOCATE(tq(3,2,3,2),         &  !local in-plane deformation gradient at face sides
             cq(3,2),             &  !average metric tensor at each face
             bq(3,6,3,2))            !in-plane B-matrix at face centers
  END IF

  ALLOCATE( bmat(3,nn,6), x(3,nn), rl(3,nn) )   !Global B-matrix, coordinates and equivalent nodal forces
  bmat = 0d0                          !initializes B-matrix
  osec    = -1    !no previous section

  ELEMS: DO piter=lim_inf,lim_sup    !loop over the elements ()
    IF (ierr == 1) CYCLE ELEMS       ! if there is an error CYCLE TO END (RETURN or EXIT are not allowed in parallel)

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
        kh    = mat%matdef(5)           !kinematic hardening model
        ! ------------------------------------------------------------------------------------
        ! ONLY MTYPES 1, 5, 6, 8  ARE ALLOWED
        SELECT CASE (mtype)
        ! ------------------------------------------------------------------------------------
        CASE (1)                        !Isotropic Elastic - Anisotropic Plastic  ( mtype == 1)
          km = mat%prope(4)*3d0              ! bulk modulus * 3
          gm = mat%prope(3)                  ! shear modulus
          bbar = sbbar
          IF( .NOT.elast ) THEN
            ! e1, nu1, uniaxial, efren, consn, r, exponent m, hill 79
            props(1:7) = mat%propp(1:7)        !
            ! for anisotropic (PLANAR) plasticity
            ch = mat%propp(16:21)              !Hill yield function coefficients
            ! A matrix, derivative of yield function (associative only)
            yf(1) =  ch(2)+ch(3) + 0.5d0       !0.5d0 added to make it diagonal
            yf(2) = -ch(3)       + 0.5d0       ! in case of von Mises
            yf(3) = -ch(2)       + 0.5d0       ! it is not essential
            yf(4) =  ch(1)+ch(3) + 0.5d0
            yf(5) = -ch(1)       + 0.5d0
            yf(6) =  ch(1)+ch(2) + 0.5d0
            yf(7) = 2d0*ch(6)
            yf(8) = 2d0*ch(4)
            yf(9) = 2d0*ch(5)
            ! for non-associative plasticity
            ch = mat%propp(22:27)              !Hill yield function coefficients
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
          bbar  = .FALSE.                   !coupled shear and volumetric part
          plast = .FALSE.
        ! ------------------------------------------------------------------------------------
        CASE (6)                        !Hyperelastic isotropic rubber            ( mtype == 6)
          eprop(1:12) = mat%prope(7:18)
          bbar  = .FALSE.                   !coupled shear and volumetric part
          elast = .TRUE.
          plast = .FALSE.                   !consider plasticity ?
        ! ------------------------------------------------------------------------------------
        CASE (8)                        !Isotropic visco-elastic                  ( mtype == 8)
          km = mat%prope(4)*3d0              ! bulk modulus * 3
          gm = mat%prope(3)                  ! shear modulus
          aux= mat%prope(11)/2d0             ! viscosity parameter
          bbar = .TRUE.
          plast = .FALSE.
          elast = .TRUE.
          !viscous = .TRUE.
        ! ------------------------------------------------------------------------------------
        END SELECT
        IF( plast ) plast = .NOT.elastic  !for elastic behavior with previous internal variables

        !************************** LAYERED SOLID *************************************************
      ELSE  !nlayr >= 1      layered solid (at integration point level)
        thick = sec%rprop(1)             !original thickness
        minstr= sec%rprop(4)             !strain threshold to use TTTI

        nl=1
        elast = sec%iprop(3) == 0               !elastic problem
        plast = .NOT.elast .OR. .NOT.elastic    !consider plasticity ?
        dm = sec%rprop(6:41)                 !linear elastic integrated matrix
        IF( plast )THEN  !plastic lamina
          IF( ALLOCATED ( rr ) )DEALLOCATE( thickl,zposl,rr,matsl,lvari )
          ALLOCATE( thickl(nlayr),zposl(nlayr),rr(5,nlayr),matsl(nlayr),lvari(nlayr) )
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
            matsl(il)%p => mat
          END DO
        END IF
        oldm = -1                              !recompute constants
      END IF
      newmt = .FALSE.                        !same material than previous ?
      IF( ALLOCATED(varin) ) DEALLOCATE (varin)
      IF( nvar > 0 ) ALLOCATE( varin(nvar,ngaus) )  !

    END IF !newmt

    x(:,1:nnode) = coora(:,e%lnods(1:nnode))       !*** computes actual configuration
    IF( quad )THEN         !for quadratic approach
      DO i=1,3               !for each adyacent element
        sides(i) = ( e%lnods(i+6) /= 0 )     !side exist?
        IF( sides(i) )THEN   !if side exist get extra nodes coordinates
          x(:,i+6:i+9:3) = coora(:,e%lnods(i+6:i+9:3))
        ELSE                 !else set those coordinates to zero
          x(:,i+6:i+9:3) = 0d0
        END IF
      END DO
      cq = 0d0  ! initializes average in-plane metric tensor at each face
      DO ii=1,2 ! for each face
        i = 3*ii-2      !first node of the face
        j = i+5
        DO jj=1,3  !for each mid side point
          tq(:,1:2,jj,ii) = MATMUL(x(:,i:i+2),e%cdq(1:3,1:2,jj,ii))
          IF( sides(jj) )THEN  !add extra nodes if side element exist
            tq(:,1,jj,ii) = tq(:,1,jj,ii) + x(:,j+jj)*e%cdq(4,1,jj,ii)
            tq(:,2,jj,ii) = tq(:,2,jj,ii) + x(:,j+jj)*e%cdq(4,2,jj,ii)
          END IF
          cq(1,ii) = cq(1,ii) + DOT_PRODUCT(tq(:,1,jj,ii),tq(:,1,jj,ii)) !C11
          cq(2,ii) = cq(2,ii) + DOT_PRODUCT(tq(:,2,jj,ii),tq(:,2,jj,ii)) !C22
          cq(3,ii) = cq(3,ii) + DOT_PRODUCT(tq(:,1,jj,ii),tq(:,2,jj,ii)) !C12
        END DO
      END DO
      cq = cq/3d0               !average metric tensor at faces
      CALL bmat16q(e%cdq,bq,tq) !in-plane b-matrix at faces
    END IF
    IF( shell )THEN !compute transvere shear strain components
      DO k=1,ngaus
        tas(:,:) = MATMUL(x(:,1:6),e%nfdas(:,:,k))  !der x/x_3
        gas(:,1) = (-x(:,2)+x(:,3))*sh(1,k)+(-x(:,5)+x(:,6))*sh(2,k) !point 1 der x/tan SQR(2)
        gas(:,2) = (-x(:,3)+x(:,1))*sh(1,k)+(-x(:,6)+x(:,4))*sh(2,k) !point 2 der -x/eta
        gas(:,3) = (-x(:,1)+x(:,2))*sh(1,k)+(-x(:,4)+x(:,5))*sh(2,k) !point 3 der x/xita
        DO i=1,3  !for each assumed strain point compute x_tan . x_X_3
          vas(i,k) = DOT_PRODUCT(tas(:,i),gas(:,i))
        END DO
        CALL bsma16(gas,tas,e%nfdas(:,:,k),bs(:,:,:,k),sh(:,k)) !compute B_s matrix at sampling points
      END DO
    END IF
    rl = 0d0                               !Initializes equivalent nodal forces
    IF(bbar)deltm = 0d0                    !initializes meam volumetric strain
    !   check
    flag2 = .FALSE.    !initializes to elastic step
    IF( ASSOCIATED(e%gausv) )THEN  !if internal variables already exist
      varin = e%gausv                !pass to auxiliar array
    ELSE IF( nvar > 0 ) THEN
      varin = 0d0                    !set internal variables to zero
    END IF
    !   Gauss point loop
    DO k=1,ngaus
      ! compute deformation gradient F and U^2=tp(F) F
      g = MATMUL(x(:,1:nnode),e%cartd(:,:,k))  !deformation gradient
      !  in-plane components
      IF( quad )THEN !interpolate from lower and upper surfaces
        stran(1) = sh(1,k)*cq(1,1) + sh(2,k)*cq(1,2)              !U^2(1,1)
        stran(2) = sh(1,k)*cq(2,1) + sh(2,k)*cq(2,2)              !U^2(2,2)
        stran(4) = sh(1,k)*cq(3,1) + sh(2,k)*cq(3,2)              !U^2(1,2)
      ELSE
        stran(1) = g(1,1)*g(1,1) + g(2,1)*g(2,1) + g(3,1)*g(3,1)  !U^2(1,1)
        stran(2) = g(1,2)*g(1,2) + g(2,2)*g(2,2) + g(3,2)*g(3,2)  !U^2(2,2)
        stran(4) = g(1,1)*g(1,2) + g(2,1)*g(2,2) + g(3,1)*g(3,2)  !U^2(1,2)
      END IF
      ! out-of-plane components
      stran(3) = g(1,3)*g(1,3) + g(2,3)*g(2,3) + g(3,3)*g(3,3)  !U^2(3,3)
      IF( shell )THEN  !shell version
        staux(3) = SUM(vas(:,k))/3d0                               !auxiliar value
        staux(1) = -staux(3) + vas(3,k)                            !natural xita-x_3
        staux(2) = +staux(3) - vas(2,k)                            !natural eta-x_3
        stran(5) = e%jacin(1,1,k)*staux(1)+e%jacin(2,1,k)*staux(2) !x_1-x_3
        stran(6) = e%jacin(1,2,k)*staux(1)+e%jacin(2,2,k)*staux(2) !x_2-x_3
      ELSE  !standard version
        stran(5) = g(1,1)*g(1,3) + g(2,1)*g(2,3) + g(3,1)*g(3,3)  !U^2(1,3)
        stran(6) = g(1,2)*g(1,3) + g(2,2)*g(2,3) + g(3,2)*g(3,3)  !U^2(2,3)
      END IF

      flag(k) = 2                   ! initializes to no spectral decomposition

      elen = e%dvol(k)**(1d0/3d0)   ! ELEMENT FRACTURE LENGHT (initial)
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
            ! stran = twice the deviatoric strains    Delta = Ln(J)/3
            CALL lgst18(stran,r,lb,'RESVPL',ierr,flag(k),delta)
            IF( ierr == 1 )CYCLE ELEMS
          END IF
          !elastic (trial) strains
          IF( .NOT.elast )THEN
            strpl(1:6) = varin(1:6,k)            !previous (twice) plastic strains
            stran(1:6) = stran(1:6) - strpl(1:6)   !trial Elastic log strains
          END IF
          str = gm*stran                         !Trial elastic shear stresses

          IF( plast ) THEN              !Consider plasticity
            efpst = varin(7,k)        !effect plastic strain
            IF( kh == 1 )THEN
              CALL corr18(str(1),str(2),str(3),str(4),str(5),str(6), &
                          efpst,gm,props,yf,fr,ierr,dstpl,pflag,is,numpt,val)
            ELSE
              CALL corr18k(str,dstpl,efpst,varin(8:13,k),gm,props,yf,fr,ierr,pflag,is,kh,numpt,val)
            END IF
            IF(ierr == 1) CYCLE ELEMS            !no convergence in plasticity
            IF( pflag )THEN                      !if plastic flow
              IF( updiv )THEN
                varin(1:6,k) = varin(1:6,k) + dstpl(1:6)  !total plastic shear strains
                varin(7,k) = varin(7,k) + efpst           !Total effect. plastic strain
                flag2 = .TRUE.
              END IF
              newiv = .TRUE.
            END IF
          END IF
          IF(bbar) THEN
            deltm = deltm + delta*e%dvol(k)     !use mean volumetric strain
          ELSE
            str(1:3) = str(1:3) + delta*km      !add press to deviatoric stress tensor
          END IF
        ! ------------------------------------------------------------------------------------
        CASE (5)                        !Orthotropic Elastic - No Plasticity      ( mtype == 5)
          IF( small )THEN               !Green Lagrange strains (twice)
            stran(1:3) = (stran(1:3)-1d0)/2d0   !2nd Piola Kirchhoff stress
          ELSE   ! Hencky (logarithmic) strains (volumetric and deviatoric parts)
            ! stran = twice the deviatoric strains    Delta = Ln(J)/3
            CALL lgst18(stran,r,lb,'RESVPL',ierr,flag(k),delta)
            stran(1:3) = stran(1:3)/2 + delta
          END IF
          str(1) = eprop(1)*stran(1)+eprop(2)*stran(2)+eprop(3)*stran(3)
          str(2) = eprop(2)*stran(1)+eprop(4)*stran(2)+eprop(5)*stran(3)
          str(3) = eprop(3)*stran(1)+eprop(5)*stran(2)+eprop(6)*stran(3)
          str(4) = eprop(7)*stran(4)
          str(5) = eprop(8)*stran(5)
          str(6) = eprop(9)*stran(6)
        ! ------------------------------------------------------------------------------------
        CASE (6)                        !Hyperelastic isotropic rubber            ( mtype == 6)
          CALL eige18(stran,r,lb,flag(k))          !squared eigenvalues
          CALL rubber3dn(eprop,lb,mat%matdef(8),str)   !compute principal stresses and pressure
          IF( bbar ) deltm = deltm + str(4)*e%dvol(k)
        ! ------------------------------------------------------------------------------------
        CASE (8)                        !Isotropic DAMAGE MODEL - XMG            ( mtype == 7)
           delta = (SUM(stran(1:3))-3d0)/6d0       !volumetric strain /3
           stran(1:3) = stran(1:3)-delta*2d0-1d0   !shear strains (twice)
           deltm = deltm + delta*e%dvol(k)         !use mean volumetric strain
           str = gm*stran + dstpl                  !visco-elastic shear stresses

        ! ------------------------------------------------------------------------------------
        END SELECT

      !************************** LAYERED SOLID *************************************************
      ELSE  !nlayr >= 1      layered solid (at integration point level)

        ! check if laminate integration is necessary
        IF (elast) THEN
          ttti = .FALSE.
        ELSE
          aux  = ABS(stran(1)-1d0)+ABS(stran(2)-1d0)+ABS(stran(3)-1d0)+ &
                 ABS(stran(4))+ABS(stran(5))+ABS(stran(6))
          ttti = aux  > minstr   !laminate integration?
        END IF
        pflag = .FALSE.                    !initializes flag

        IF( ttti )THEN
          IF( small )THEN               !Green Lagrange strains (twice)
            stran(1:3) = (stran(1:3)-1d0)*0.5   !shear strains (twice)
          ELSE   ! Hencky (logarithmic) strains (volumetric and deviatoric parts)
            !computes Hencky (logarithmic) strains (volumetric and deviatoric parts)
            ! stran = twice the deviatoric strains    Delta = Ln(J)/3
            CALL lgst18(stran,r,lb,'RESVPL',ierr,flag(k),delta)
            IF( ierr == 1 )   CYCLE ELEMS
            stran(1:3) = (stran(1:3)/2d0+delta)   !recompute axial strains
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

            mat => matsl(il)%p        ! pointer to layer material

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

        ELSE !elastic response of integrated laminated solid (small strains only)
          stran(1:3) = (stran(1:3) - 1d0 )*0.5   !2nd Piola Kirchhoff stress
          str(1:6) = MATMUL(RESHAPE((dm),(/6,6/)),stran(1:6))
        END IF

      END IF
      ! store auxiliar values for second loop
      stres(:,k)  =  str  !shear stresses
      gv(:,:,k) =  g      !deformation gradient
      rv(:,:,k) =  r      !eigenvectors and
      lbv(:,k)  =  lb     !eigenvalues
    END DO  !k (loop to next integration point)

    IF( updiv )THEN        !if internal variables have to be updated
      IF( flag2 )THEN         !if there are variables to modify
        IF(.NOT.ASSOCIATED(e%gausv))   ALLOCATE( e%gausv(nvar,ngaus))
        e%gausv = varin
      END IF
    END IF

    IF( bbar )THEN
      IF( mtype == 1 .OR. mtype == 8)THEN
        deltm = deltm/SUM(e%dvol)   !mean volumetric strain / 3
        press = km*deltm            !mean pressure
        stres(1:3,:) = stres(1:3,:) + press  !sum on diagonal stresses for all Gauss points
      ELSE IF( mtype == 6 )THEN
        press = deltm/SUM(e%dvol)   !mean volumetric press / 3
        stres(4,:) = press
      END IF
    END IF

    DO k=1,ngaus         !second Gauss loop
      ! recovers Gauss point values
      str = stres(:,k)   !stresses
      lb = lbv(:,k)      !eigenvalues (principal stretchings)
      r  = rv(:,:,k)     !eigenvectors
      IF(mtype == 6)THEN  !rubbers
        !Second Piola Kirchhoff stress tensor
        DO i=1,3
         str(i) = (str(i)+str(4)/lb(i))/lb(i)     ! S = [ Phid_i + (Li Phiv_i)/Li ]/Li
        END DO
        ! store the stress vector (2nd Piola-Kirchhoff stress tensor)
        !Compute 3-D stresses with the spectral decomposition  R^T * Beta * R
        sig(1) = r(1,1)*str(1)*r(1,1) + r(1,2)*str(2)*r(1,2) + r(1,3)*str(3)*r(1,3) !(1,1)
        sig(2) = r(2,1)*str(1)*r(2,1) + r(2,2)*str(2)*r(2,2) + r(2,3)*str(3)*r(2,3) !(2,2)
        sig(3) = r(3,1)*str(1)*r(3,1) + r(3,2)*str(2)*r(3,2) + r(3,3)*str(3)*r(3,3) !(3,3)
        sig(4) = r(1,1)*str(1)*r(2,1) + r(1,2)*str(2)*r(2,2) + r(1,3)*str(3)*r(2,3) !(1,2)
        sig(5) = r(1,1)*str(1)*r(3,1) + r(1,2)*str(2)*r(3,2) + r(1,3)*str(3)*r(3,3) !(1,3)
        sig(6) = r(2,1)*str(1)*r(3,1) + r(2,2)*str(2)*r(3,2) + r(2,3)*str(3)*r(3,3) !(2,3)
        e%stint(:,k) = sig ! keep 2nd Piola Kirchhoff
      ELSE
        e%stint(:,k) = str ! keep average Kirchhoff stress for output
        IF( flag(k) == 0)THEN
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
        ELSE IF(flag(k) == 1)THEN  !all eigenvalues are almost equal
          sig(:) = str(:)/lbv(1,k) !lambda squared
        ELSE
          sig = str
        END IF
      END IF

      !***   compute Int(B**t*sigma) on element level

      sig = sig*e%dvol(k)     !Original dvol used (TLF)

      !       actual configuration B Matrices
      IF( shell ) THEN
        CALL bmat16(e%cartd(:,:,k),bmat,gv(:,:,k),quad,bq,sh(:,k),shell,bs(:,:,:,k),e%jacin(:,:,k))
      ELSE
        CALL bmat16(e%cartd(:,:,k),bmat,gv(:,:,k),quad,bq,sh(:,k),shell)
      END IF
      CALL proma3(rl(1,1),bmat(1,1,1),sig(1),ndofe,6,.FALSE.)  !nodal contributions
    END DO

    DO i=1,nn             !for each node in the element
      n = e%lnods(i)     !node
      IF( n > 0 ) THEN
        !$ m = MOD(n,MAXLOCK)
        !$ CALL omp_set_lock(lock1(m))    !Start blocking access in shared variables
        resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
        !$ CALL omp_unset_lock(lock1(m))  !End blocking access in shared variables
      END IF
    END DO

    e => e%next          !point to next element in the list
  END DO ELEMS
  DEALLOCATE (bmat,x,rl)
  IF( quad ) DEALLOCATE(tq,cq,bq)
  IF(ASSOCIATED(val)) NULLIFY(val)
  IF(ASSOCIATED(matsl)) DEALLOCATE(matsl)
  IF(ALLOCATED(lvari)) DEALLOCATE(lvari)
  IF(ALLOCATED(thickl)) DEALLOCATE(thickl,zposl,rr)
  IF(ALLOCATED(varin)) DEALLOCATE(varin)
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

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE resv16
