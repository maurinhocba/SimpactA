 SUBROUTINE resv04( nele0,head0, coora, resid, ierr, ttime, small, ngaus, beta)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT T-L 3-D SOLID
 !
 !-------------------------------------------------------------------
 USE gvar_db, ONLY : updiv,newiv,elastic
 USE ele18_db, ONLY : lgst18,eige18,rubber3d
 !$ USE omp_lib
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(OUT) :: ierr        !flag for error detection
 INTEGER (kind=4), INTENT(IN) :: ngaus,     & !number of integration points
                                 nele0        !number of elements in the set
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             beta(:),       & !stabilization factors
                             ttime            !present total time
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
 TYPE (ele04), POINTER  :: head0              !pointer to first element
 LOGICAL, INTENT(IN) :: small      !.TRUE. if Green strains will be used

 ! local variables
 INTEGER (kind=4), PARAMETER :: ndime=3
 ! variables for assumed strain transverse shear
 REAL (kind=8) :: tas(3,4),      & ! der x/x_3 at Assumed Strain points
                  gas(3,4),      & ! der x/tan at Assumed Strain points
                  vas(4,2),va(4),& ! strains at Assumed Strain Points
                  bs(3,8,2,2),bh(24,2,2)   ! Bs matrix at face centers, Stab matrix
 TYPE (ele04), POINTER :: e => NULL() !pointer to an element data
 TYPE (ele04), POINTER :: heade => NULL() !pointer to an element data
 LOGICAL :: newmt,   & ! TRUE if material constant computation necessary
            elast,   & ! TRUE if material is elastic
            plast      ! TRUE if plasticity is to be considered
 INTEGER (kind=4) :: osec,     & !associated material of previous element
                     numpt,    & !number of points in curve
                     ndofe,    & !number of DOFs
                     mtype         !associated material type
 REAL (kind=8) ::  fr(9),yf(9),& !Flow rule and yield function coefficients
                   km,gm,      & !bulk modulus and shear modulus
                   eprop(12),  & !elastic Material properties
                   props(7),   & !plastic Material properties
                   thick,      & !thickness (original)
                   minstr,     & !minimum strain to integrate trougth the thickness
                   dm(36),     & !Elastic integrated constitutive matrix
                   ch(6)         !Hill-48 coefficients

 INTEGER (kind=4) ielem,nelem,       & !element number
                  is,kh,             & !isotropic and kinematic hardening models
                  isec,              & !associated material
                  flag,              & !flag for eigen-decomposition
                  i,ii,jj,j,g,n,gg     !different indexes

 INTEGER (kind=4) nl,nlayr,  & !number of layers
                  oldm,      & !old material label
                  mate,      & !material label
                  ilv,jlv,   & !indexes for layer varin
                  il,        & !layer index
                  nvar         !total number of internal variables

 REAL (kind=8) stral(6),   & !local strains
               strel(6),   & !local streses
               staux(6),   & !local streses auxiliar
               str(6),     & !stresses (different measures)
               sig(6),     & !stresses (different measures)
               stran(6),   & !C=U^2  also Log strains
               aux,voli,   & !auxiliar value
               dstpl(6),   & !plastic strain increments
               efpst,      & !effective plastic strain
               r(3,3),lb(3),& !eigen-pair in local system
               delta         !volumetric strain/3, element volume
 REAL (kind=8) :: f3(3),c33,c33z,b3(3,nnode),ka,fac, &  !EAS
                  seq(11),ralpha,aprim,d33,d33e,stm,sts
 !REAL (kind=8) :: u(ndime,nnode),          &  !incremental displacements
 REAL (kind=8) :: x(ndime,nnode),          &  !element coordinates
                  rl(ndime,nnode),         &  !nodal internal equivalent forces
                  cq(3,2)                     !average metric tensor at each face
 REAL (kind=8) :: tq(3,2,4,2), &  !(ndime,2,3,2) local in-plane deformation gradient at face sides
                  bq(3,4,3,2), &  !(ndime,4,3,2) in-plane B-matrix at face centers
                  bm(12,4,2),eh(4,2),ehb(2),ff(3,2),fff(3)  !in-plane B12-matrix at mid-point sides for stabilization

 REAL (kind=8) :: thf(ngaus), wei(ngaus), sh(3,ngaus), shape(2,ngaus)   !thickness gauss integration

 LOGICAL :: pflag,   &! TRUE if plastic flow in the step
            ttti,    &! TRUE if Trougth The Thickness Integration is necessary
            flag2     ! TRUE if plastic at any gauss point

 REAL (kind=8), POINTER :: val(:,:) => NULL() ! curve for IHPOINT
 REAL (kind=8)          :: elen     ! Damage length - For damage subroutine - XMG

 TYPE (section), POINTER :: sec => NULL() !pointer to a section data
 TYPE (mater), POINTER   :: mat => NULL() !pointer to a material data

 TYPE (pmater), POINTER :: matsl(:) => NULL()   !array of pointers to materials data
 ! allocatable arrays : MUST be replaced by any STATIC variables to avoid run-time memory increase
 INTEGER (kind=4), ALLOCATABLE ::  lvari(:)
 REAL (kind=8), ALLOCATABLE ::  thickl(:),zposl(:),rr(:,:)
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


 heade => head0          !auxiliar pointer
 nelem = nele0           !auxiliar counter
 ndofe = ndime*nnode     ! number of DOFs
 !integration points through the thickness
 CALL gaussq(ngaus,thf(1),wei(1))
 wei = wei*4d0
 ! compute shape(2,ngaus) to interpolate
 ! compute sh(3,ngaus) to evaluate jacobian
 DO gg=1,ngaus
   sh(1,gg) = -0.5d0*(1d0-thf(gg))*thf(gg)  !quadratic 1-D Lagrangian polynomials (for jacobian)
   sh(2,gg) =  1d0-thf(gg)**2
   sh(3,gg) =  0.5d0*(1d0+thf(gg))*thf(gg)
   shape(1,gg) =  0.5d0*(1d0-thf(gg))      !linear 1-D Lagrangian polynomials (for interpolation)
   shape(2,gg) =  0.5d0*(1d0+thf(gg))
 END DO

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

  !$OMP  PARALLEL                                                                             &
  !$OMP& SHARED (psecs,pmats, updiv, ngaus, coora, ttime, small, beta, newiv, ierr, resid,    &
  !$OMP&          heade, nelem, ndofe, num_chk,lock1,elastic,thf,wei,sh,shape,gg)             &
  !$OMP& PRIVATE(tas, gas, vas, bs, e, newmt, elast, plast, osec, numpt, mtype, fr, yf, kh,   &
  !$OMP&         va, bh, stm, sts, bm, eh, ff, fff, ehb,                                      &
  !$OMP&         km, gm, eprop, props, thick, minstr, dm, ch, ielem, is, isec, flag, i,       &
  !$OMP&         ii, jj, j, g, n, m, nl, nlayr, oldm, mate, ilv, jlv, il, nvar, stral, strel, &
  !$OMP&         staux, str, sig, stran, aux, voli, dstpl, efpst, r, lb, delta, f3,           &
  !$OMP&         c33, c33z, b3, ka, fac, seq, ralpha, aprim, d33, d33e, x, rl, cq, tq,        &
  !$OMP&         bq, pflag, ttti, flag2, val, elen, sec, mat, matsl, lvari,                   &
  !$OMP&         thickl,zposl,rr,varin,piter, lim_inf, lim_sup)

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
 NULLIFY(val); NULLIFY(matsl) !; NULLIFY(lvari); NULLIFY(thickl)
 osec = -1                           !no previous material

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
       ! ONLY MTYPES 1, 5, 6 ARE ALLOWED
       SELECT CASE (mtype)
       ! ------------------------------------------------------------------------------------
       CASE (1)                        !Isotropic Elastic - Anisotropic Plastic  ( mtype == 1)
         km = mat%prope(4)*3d0              ! bulk modulus * 3
         gm = mat%prope(3)                  ! shear modulus
         stm = 2d0*gm*beta(1)
         sts = gm*beta(2)/3d0
         d33e = (km + gm*4d0)/3d0           ! Elastic D33
         IF( .NOT.elast ) THEN
           ! 1:4 iso-hard pars, 5: death-time, 6:7 kin-hard pars
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
         stm = 4d0*eprop(7)*beta(1)
         sts = 2d0/3d0*eprop(8)*beta(2)
         d33 = eprop(6)
         plast = .FALSE.
       ! ------------------------------------------------------------------------------------
       CASE (6)                        !Hyperelastic isotropic rubber            ( mtype == 6)
         eprop(1:12) = mat%prope(7:18)
         elast = .TRUE.
         plast = .FALSE.                   !consider plasticity ?
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
     IF( plast) ALLOCATE( varin(nvar,ngaus) )  !

   END IF !newmt

   !*** computes actual configuration
   x(:,1:nnode) = coora(:,e%lnods(1:nnode))               !coordinates
   !u(:,1:nnode) = x(:,1:nnode) - coori(:,e%lnods(1:nnode))  !incremental displacements
   cq = 0d0  ! initializes average in-plane metric tensor at each face
   ! compute variables associated with in-plane behavior
   DO ii=1,2 ! for each face
     i = 4*ii-3         !first node of the face (1 or 5)
     DO jj=1,4  !for each mid side point
       tq(:,1:2,jj,ii) = MATMUL(x(:,i:i+3),e%cdq(1:4,1:2,jj,ii))
       ! add to average metric tensor
       cq(1,ii) = cq(1,ii) + DOT_PRODUCT(tq(:,1,jj,ii),tq(:,1,jj,ii))*e%area(jj,ii) !C11
       cq(2,ii) = cq(2,ii) + DOT_PRODUCT(tq(:,2,jj,ii),tq(:,2,jj,ii))*e%area(jj,ii) !C22
       eh(jj,ii) =           DOT_PRODUCT(tq(:,1,jj,ii),tq(:,2,jj,ii)) !C12
       cq(3,ii) = cq(3,ii) + eh(jj,ii)*e%area(jj,ii)                                !C12
     END DO
     eh(:,ii) =  eh(:,ii) - cq(3,ii) !C12-H
     !compute transvere shear strain components
     tas(:,:) = MATMUL(x,e%nfdas(:,:,ii))  !der x/x_3  at the 4 sides
     IF( ii == 1 )THEN
       gas(:,1) = (x(:,4) -x(:,1))/2d0 !point A der x/eta
       gas(:,2) = (x(:,2) -x(:,1))/2d0 !point B der x/xita
       gas(:,3) = (x(:,3) -x(:,2))/2d0 !point C der x/eta
       gas(:,4) = (x(:,3) -x(:,4))/2d0 !point D der x/xita
       CALL bsma04(gas,tas,e%nfdas(:,:,1),bs(:,:,:,1),e%jacin(:,:,1),1,bh(:,:,1))     !compute B_s matrix at sampling points
     ELSE
       gas(:,1) = (x(:,8) -x(:,5))/2d0 !point A der x/eta
       gas(:,2) = (x(:,6) -x(:,5))/2d0 !point B der x/xita
       gas(:,3) = (x(:,7) -x(:,6))/2d0 !point C der x/eta
       gas(:,4) = (x(:,7) -x(:,8))/2d0 !point D der x/xita
       CALL bsma04(gas,tas,e%nfdas(:,:,2),bs(:,:,:,2),e%jacin(:,:,2),2,bh(:,:,2))     !compute B_s matrix at sampling points
     END IF
     !ehb(ii) = 0d0
     DO g=1,4  !for each assumed strain point compute x_tan . x_X_3
       va(g) = DOT_PRODUCT(tas(:,g),gas(:,g))
       !ehb(ii) = ehb(ii) + ABS(eh(g,ii))*e%area(g,ii)
     END DO
     vas(3,ii) =(va(2)+va(4))/2d0                                    !C_xi3
     vas(4,ii) =(va(1)+va(3))/2d0                                    !C_eta3
     vas(1,ii) = e%jacin(1,1,ii)*vas(3,ii)+e%jacin(2,1,ii)*vas(4,ii) !x_1-x_3
     vas(2,ii) = e%jacin(1,2,ii)*vas(3,ii)+e%jacin(2,2,ii)*vas(4,ii) !x_2-x_3
     vas(3,ii) =-va(2)+va(4)                                         !C_xi3_H
     vas(4,ii) =-va(1)+va(3)                                         !C_eta3_H
   END DO

   CALL bmat04q(e%cdq,bq,tq,bm,e%area) !in-plane b-matrix at faces and stabilization
   ! EAS
   f3 = MATMUL(x(:,1:nnode),e%cartd)     !f_3 at element center
   c33 = DOT_PRODUCT(f3,f3)              !C(3,3) at element center
   DO i=1,nnode
     b3(:,i) = f3*e%cartd(i)           !b(:,3) matrix at element center
   END DO
   !   check
   flag2 = .FALSE.    !initializes to elastic step
   IF( ASSOCIATED(e%gausv) )THEN  !if internal variables already exist
     varin = e%gausv                !pass to auxiliar array
   ELSE IF( plast ) THEN
     varin = 0d0                    !set internal variables to zero
   END IF
   !   Gauss point loop
   ! initializes variables to be computed over GP loop
   ka = 0d0                 !factor Ka
   ralpha = 0d0             !EAS residual
   seq = 0d0                !initializes integrated equivalent stresses
   ff = 0d0                 !initializes stabilization factors due to plasticitiy
   DO g=1,ngaus      !for each Gauss point
     fac = EXP(2d0*e%alpha*thf(g))                             !
     voli = DOT_PRODUCT(sh(:,g),e%dvol)*wei(g)                 !gauss point volume
     !interpolate metric tensor from lower and upper surfaces
     stran((/1,2,4/)) = shape(1,g)*cq(:,1) + shape(2,g)*cq(:,2) !  in-plane components
     stran(5:6) = shape(1,g)*vas(1:2,1) + shape(2,g)*vas(1:2,2) !transverse shear components
     c33z = c33*fac                                             ! out-of-plane components
     stran(3) = c33z                                            !U^2(3,3)
     IF( small )THEN
       flag = 2                   ! no spectral decomposition
     ELSE
       flag = 1                   ! initializes to no spectral decomposition
     END IF
     elen = voli**0.333d0           ! ELEMENT FRACTURE LENGHT (initial)
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
           CALL lgst18(stran,r,lb,'RESVPL',ierr,flag,delta)
           IF( ierr == 1 )CYCLE ELEMS
         END IF
         ! in-plane stabilization factors
         fff = 1d0
         IF( plast ) THEN
           IF( ABS(stran(4)) > 1d-12 ) fff(1) = MAX(1d0 - ABS(varin(4,g)/stran(4)),beta(3))
           IF( ABS(stran(5)) > 1d-12 ) fff(2) = MAX(1d0 - ABS(varin(5,g)/stran(5)),beta(3))
           IF( ABS(stran(6)) > 1d-12 ) fff(3) = MAX(1d0 - ABS(varin(6,g)/stran(6)),beta(3))
         END IF
         ff(:,1) = ff(:,1) + fff*shape(1,g)*wei(g)
         ff(:,2) = ff(:,2) + fff*shape(2,g)*wei(g)
          !elastic (trial) strains
         IF( .NOT.elast ) stran(1:6) = stran(1:6) - varin(1:6,g) !trial Elastic log strains
         str = gm*stran                         !Trial elastic shear stresses
         d33 = d33e
         IF( plast ) THEN              !Consider plasticity
           efpst = varin(7,g)        !effect plastic strain
           IF( kh == 1 )THEN
             CALL corr18(str(1),str(2),str(3),str(4),str(5),str(6), &
                         efpst,gm,props,yf,fr,ierr,dstpl,pflag,is,numpt,val,aprim)
           ELSE
             CALL corr18k(str,dstpl,efpst,varin(8:13,g),gm,props,yf,fr,ierr,pflag,is,kh,numpt,val,aprim)
           END IF
           IF(ierr == 1) CYCLE ELEMS            !no convergence in plasticity
           IF( pflag )THEN                      !if plastic flow
             varin(1:6,g) = varin(1:6,g) + dstpl(1:6)  !total plastic shear strains
             varin(7,g) = varin(7,g) + efpst           !Total effect. plastic strain
             flag2 = .TRUE.
             !d33 = km/3d0+aprim                 ! commented in last version
           END IF
           d33 = d33/c33z**2
           newiv = .TRUE.
         END IF
         str(1:3) = str(1:3) + delta*km      !add press to deviatoric stress tensor
       ! ------------------------------------------------------------------------------------

       CASE (5)                        !Orthotropic Elastic - No Plasticity      ( mtype == 5)
         IF( small )THEN               !Green Lagrange strains (twice)
           stran(1:3) = (stran(1:3)-1d0)/2d0   !2nd Piola Kirchhoff stress
         ELSE   ! Hencky (logarithmic) strains (volumetric and deviatoric parts)
           ! stran = twice the deviatoric strains    Delta = Ln(J)/3
           CALL lgst18(stran,r,lb,'RESVPL',ierr,flag,delta)
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
         CALL rubber3d(eprop,stran,mat%matdef(8),stre=str,d33=d33)

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
           CALL lgst18(stran,r,lb,'RESVPL',ierr,flag,delta)
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

     END IF  !nlayr == 0

     IF(mtype == 6)THEN  !rubbers
       sig = str
       e%stint(:,g) = sig ! keep 2nd Piola Kirchhoff
     ELSE
       e%stint(:,g) = str ! keep Kirchhoff stress for output
       IF( flag == 0)THEN
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
       ELSE IF(flag == 1)THEN  !all eigenvalues are almost equal
         sig(:) = str(:)/lb(1) !lambda squared
       ELSE
         sig = str
       END IF
     END IF

     !***   compute Int(B**t*sigma) on element level

     sig = sig*voli          !Original dvol used (TLF)

     ! Compute equivalent stresses and constitutive matrix
     seq( 1: 3) = seq( 1: 3) + sig((/1,2,4/))*shape(1,g)
     seq( 4: 6) = seq( 4: 6) + sig((/1,2,4/))*shape(2,g)
     seq( 7: 8) = seq( 7: 8) + sig(5:6)      *shape(1,g)
     seq( 9:10) = seq( 9:10) + sig(5:6)      *shape(2,g)
     seq(11)    = seq(11)    + sig(3)        *fac
     ! compute  EAS residual and stiffness
     ralpha = ralpha + sig(3)*thf(g)*c33z
     ka = ka + (d33*c33z*thf(g)*voli + sig(3)*2d0*thf(g))*c33z *thf(g)

   END DO

    IF( updiv )THEN        !if internal variables have to be updated
      IF( flag2 )THEN         !if there are variables to modify
        IF(.NOT.ASSOCIATED(e%gausv))   ALLOCATE( e%gausv(nvar,ngaus))
        e%gausv = varin
      END IF
    END IF

   DO i=1,4 !for each face node
     j  = i+4  !for each other face node
     rl(:,i) = seq( 1)*bq(:,i,1,1) + seq( 2)*bq(:,i,2,1) + seq( 3)*bq(:,i,3,1) + &
               seq( 7)*bs(:,i,1,1) + seq( 8)*bs(:,i,2,1) +                       &
               seq( 9)*bs(:,i,1,2) + seq(10)*bs(:,i,2,2) + seq(11)*b3(:,i)
     rl(:,j) = seq( 4)*bq(:,i,1,2) + seq( 5)*bq(:,i,2,2) + seq( 6)*bq(:,i,3,2) + &
               seq( 7)*bs(:,j,1,1) + seq( 8)*bs(:,j,2,1) +                       &
               seq( 9)*bs(:,j,1,2) + seq(10)*bs(:,j,2,2) + seq(11)*b3(:,j)
   END DO
   e%alpha = e%alpha - ralpha/ka  !first update of EAS parameter
               !Stabilization
   IF( mtype == 6 )THEN
      stm = d33e*beta(1)  !or d33 only ?
      sts = d33e*beta(2)
   END IF
   ! compute in-plane stabilization
   ka = e%dvol(2)*stm
   DO i=1,2
     !fac = ehb(ii)/MAX(ABS(cq(3,ii)),1e-4)
     !IF( fac < beta(3) )THEN
     !  fac = 4d0
     !ELSE
     !  fac = 4d0*ehb(ii)/beta(3)
     !END IF
     DO j=1,4
       !eh(j,i) = eh(j,i)*ka*e%area(j,i)*fac
       eh(j,i) = eh(j,i)*ka*e%area(j,i)*ff(1,i)
     END DO
   END DO
   CALL stbm04r(rl,bm,eh)
   ! compute shear stabilization
   ka = e%dvol(2)*sts
   DO i=1,2
     DO j=3,4
       vas(j,i) = vas(j,i)*ka*ff(j-1,i)                             !C_xi3_H
     END DO
   END DO
   CALL stbs04r(rl,bh,vas(3:4,:),e%jacin)
   DO i=1,nnode             !for each node in the element
     n = e%lnods(i)     !node
     !$ m = MOD(n,MAXLOCK)
     !$ CALL omp_set_lock(lock1(m))    !Start blocking access in shared variables
     resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
     !$ CALL omp_unset_lock(lock1(m))  !End blocking access in shared variables
   END DO
   IF( flag2 )THEN         !if there are variables to modify
     IF(.NOT.ASSOCIATED(e%gausv))   ALLOCATE( e%gausv(nvar,ngaus))
     e%gausv = varin
   END IF
   e => e%next          !point to next element in the list
 END DO ELEMS
 IF(ALLOCATED(varin)) DEALLOCATE(varin)
  IF(ASSOCIATED(val)) NULLIFY(val)
  IF(ASSOCIATED(matsl)) DEALLOCATE(matsl)
  IF(ALLOCATED(lvari)) DEALLOCATE(lvari)
  IF(ALLOCATED(thickl)) DEALLOCATE(thickl,zposl,rr)
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
 END SUBROUTINE resv04


 SUBROUTINE stbm04r(s,b,eh)
 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT(IN) :: eh(:,:),b(:,:,:)
 REAL(kind=8), INTENT(IN OUT) :: s(:,:)
 ! local variables
 INTEGER(kind=4) :: j,k,l,m,i0,i
 REAL(kind=8) :: f

 i0 = 1   !first position in r for first face
 DO m=1,2  !for each face
   DO k=1,4        !for each mid-side point
     l = 0         !initializes position in B mat
     f = eh(k,m)
     DO i=i0,i0+3
       DO j=1,3     !for each DOF in the face
         l = l+1    !row
         s(j,i) = s(j,i) + f*b(l,k,m)
       END DO
     END DO
   END DO
   i0 = 5        !first position in r for second face
 END DO

 RETURN

 END SUBROUTINE stbm04r


 SUBROUTINE stbs04r(s,bs,eh,ji)
 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT(IN) :: eh(:,:),bs(:,:,:),ji(:,:,:)
 REAL(kind=8), INTENT(IN OUT) :: s(:,:)
 ! local variables
 INTEGER(kind=4) :: i,j,l,m
 REAL(kind=8) :: a,b

 DO m=1,2  !for each face
   a = (ji(1,1,m)*ji(1,1,m) + ji(2,1,m)*ji(2,1,m))*eh(1,m)
   b = (ji(1,2,m)*ji(1,2,m) + ji(2,2,m)*ji(2,2,m))*eh(2,m)
   l = 0         !initializes position in B mat
   DO i=1,8       !for each node
     DO j=1,3       !for each DOF
       l = l+1    !row
       s(j,i) = s(j,i) + a*bs(l,1,m) + b*bs(l,2,m)
     END DO
   END DO
 END DO

 RETURN

 END SUBROUTINE stbs04r
