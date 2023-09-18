 SUBROUTINE resv25( nele0,head0, iffix, coora, resid, logst, ierr, ttime,&
                    bottom, top, coorb, coort, ifact, nbs, bhead,  &
                    stm, stb, stint)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT TL NBSQ
 !
 !-------------------------------------------------------------------
 !$ USE omp_lib
 USE gvar_db, ONLY : updiv,newiv,elastic
 USE kinc_db, ONLY : ndepd
 !USE npoi_db, ONLY : tempe, itemp
 !USE npo_db, ONLY : coorc
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN) :: logst       !flags for small or large strain analysis
 LOGICAL, INTENT(IN) :: bottom,top  !flags for top and bottom surface computation
 INTEGER (kind=4), INTENT(IN) ::  nele0,    & !number of elements in the set
                                  iffix(:), &  !nodal fixities
                                  nbs         !number of branching sides
 INTEGER (kind=4), INTENT(OUT) :: ierr,     & !flag for error detection
                               ifact(:)       !weigths for surface computation
 REAL(kind=8),INTENT (IN) :: coora(:,:),    & !present coordinates
                             ttime,         & !present total time
                             stm,           & !membrane stabilization factor
                             stb              !bending stabilization factor
 REAL(kind=8),INTENT (IN OUT) :: resid(:,:),& !internal nodal forces
                                 coorb(:,:),& !bottom surface coordinates
                                 coort(:,:)   !top surface coordinates
 TYPE (ele25), POINTER  :: head0              !pointer to first element
 TYPE (sideb), POINTER  :: bhead              !pointer to first element
 REAL(kind=8), POINTER :: stint(:,:)       !moments and shear

 ! local variables

 TYPE (ele25), POINTER :: e    !pointer to an element data
 TYPE (ele25), POINTER  :: heade  !pointer to first element
 TYPE (sideb), POINTER  :: si    !pointer to branching side

 LOGICAL :: is(4),solsi(4),mti,   &! TRUE if side is clamped or symmetry line
            branch,  &! TRUE if a branching side exists at the element
            pstr
 INTEGER (kind=4) sides(4), &! 1 if side element exist
                  m,nbn,nbo,lnb(18) !nodal connectivities for branching element

 INTEGER (kind=4)  :: nelem,iel,iel0,i3, &    !element order
                      lna(4)      !solid side connectivities

 REAL (kind=8) fff(3,4),     & !C=U^2  hourglass strains
               stra1(9),     & !first and second fundamental forms
               strag(4),     &
               bmatx(3,12,7),& !B matrix (average)
               bb(12,3,4),   & !membrane B matrix at gauss points
               bbk(3,12,4),  &
               bbg(3,12,4),  &
               x(3,12),      & !nodal coordinates of the patch
               tb(3,8),      & !normal of the element 8adyacent element)
               tns(3,2,8),   &
               t(3,2,0:8),   & !Local triad + local in-plane vector of side elms
               t3(3),        & !element normal
               h(3,3),       & !integrated vectors
               hb(3,2),      & !integrated side vectors
               ratio,        & !stabilization modulus
               stbf,stra2(3),& !stabilization factor for bending
               rl(3,12),     & !nodal internal equivalent forces
               y(3,4),yl(2,4),yp(3,4),sh(4) !quad coordinates
 REAL (kind=8), PARAMETER :: fac(4) = (/ -1d0, 1d0, 1d0, -1d0 /)
 REAL (kind=8), ALLOCATABLE :: bbr(:,:,:),bbrk(:,:,:),rb(:,:)                     !branching sides
 REAL (kind=8) ::   dt(4),f2

 INCLUDE '../lbst/resv_p0.fpp'           !common definition

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

 heade => head0
 nelem = nele0
 iel0 = 0

 !***  compute data associated to branching sides  only if branching sides exists
 IF( nbs > 0 )  CALL bran25(nbs,bhead,coora)         !compute angles and B matrices

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
 !$OMP&         updiv, logst, bottom, top, coort, coorb, iffix, coora, ttime,              &
 !$OMP&         heade, ierr, newiv, resid, stint, ifact,stm,stb,                           &
 !$OMP&         nbs, letra, iel0                 )                                         &
 !$OMP& PRIVATE(lim_inf, lim_sup, e, si, is,solsi,branch, mti, pstr,sides, nbn, nbo, lnb,  &
 !$OMP&               iel, i3, lna, fff, stra1, strag, bmatx, x, tb, tns, t, t3, rl, sh,   &
 !$OMP&         bb, bbk, bbg, rb, h, hb, ratio, stbf, stra2, y, yl, yp, bbr, bbrk, dt, f2, &
 !$OMP&         large, ttti, pflag,                                                        &
 !$OMP&         newmt,  natst, elast, shell, plast, shear, visco, coupled,                 &
 !$OMP&         ielem, i, j, k, m, n, l, ilv, jlv, error,                                  &
 !$OMP&         isec, nl, nlayr, mtype, secty, oldm, mate, numpt, osec, nvar,       &
 !$OMP&         stres, sigma, stran, stral, strel, stine, r1,r2, lb, lc, thnew, zk,        &
 !$OMP&         aux, area1, s1, s2, efpst, strpl, dstpl, eql,                              &
 !$OMP&         thick, poiss, db, c, gh, dm, alpha, propi, chi, chib, chid,        &
 !$OMP&         deatht, minstr, min_tr, max_tr,                                            &
 !$OMP&         thf, wei, shf,                                                             &
 !$OMP&         thickl, zposl, rr, lvari, val, varin,                     &
 !$OMP&         sec, mat, matsl)
 !$ IF( omp_get_num_threads() > 1 )THEN
   !$ num_chk = (nelem+omp_get_num_threads()-1)/omp_get_num_threads()
   !$ lim_inf = MIN(nelem,1+omp_get_thread_num()*num_chk) + iel0    !Find the first element in the list
   !$ lim_sup = MIN(nelem,num_chk*(1+omp_get_thread_num())) + iel0  !Find the last element in the list
 !$ END IF
 !***  compute data associated to branching sides

 IF( nbs > 0 ) THEN                     !only if branching sides exists
   nbo = 12                              !initializes
   ALLOCATE( bbr(3,nbo,6),bbrk(3,nbo,4),rb(3,nbo) )    !get memory for auxiliar matrices
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
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   ielem = e%numel                   !element number
   isec = e%matno                   !associated material
   newmt = osec /= isec             !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     INCLUDE '../lbst/resv_p1.fpp'           !get material constant
     newmt = .FALSE.                        !same material than previous ?
     stbf = stb*dm(5)               !bending hourglass control
   END IF

   !*** computes actual configuration
   is = .FALSE.
   solsi = .FALSE.
   branch = ANY(e%lnods < 0 )    !a branching side exists
   IF( branch )THEN
     nbn = 4                     !main quadrilateral nodes
     lnb(1:4) = e%lnods(1:4)     !fourth three nodes
   END IF
   DO i =1,4                       !for all adyacent nodes
     n = e%lnods(i)                !node number
     x(1:3,i) = coora(1:3,n)       !get actual coordinates
     j = fn(i)                     !first node of associated side element (local)
     l = e%lnods(j)                !first node of associated side element (global)
     IF(l > 0) THEN               !if node exist
       sides(i) = 1               !side element exists
       x(:,j) = coora(:,l)     !get first node actual coordinates
       j = j+1                     !second (local) node
       m = e%lnods(j)              !second (global) node
       x(1:3,j) = coora(1:3,m)     !get second node actual coordinates
       IF(branch)THEN
         nbn = nbn + 1
         lnb(nbn) = l
         nbn = nbn + 1
         lnb(nbn) = m
       END IF
     ELSE IF( l == 0 )THEN     !boundary side
       sides(i) = 0                !side element does not exists
       n = e%lnods(hh(3,i))        !first node of the side
       m = e%lnods(hh(2,i))        !second node of the side
       x(1:3,nextn(i)) = coora(1:3,m) !get second node original coordinates
       is(i) = ( iffix(n) /= 0 .AND. iffix(m) /= 0 )  !clamped side ?
       IF( iffix(n) < 0 .AND. iffix(m) < 0 )THEN !node on a side
         l = -iffix(n)            !compute position of the i+1 node
         lna(1:2) = nndpd(2:3,l)
         l = -iffix(m)            !compute position of the i+2 node
         lna(3:4) = (/ nndpd(3,l),nndpd(2,l) /)
         y =  coora(:,lna(:))  ! nodes of the quadrilateral
         t(:,1,1) = -y(:,1)+y(:,2)+y(:,3)-y(:,4) !phi,xita at the center
         t(:,2,1) = x(:,nextn(i))-x(:,i)  !side vector
         CALL vecpro(t(1,1,1),t(1,2,1),x(1,j)) !outward (to the solid) normal vector
         CALL vecuni(3,x(1,j),zk)          !unit normal to the quad
         e%normd(1:3,2,i+4) = -x(:,j)         !store unit normal to the solid side
         solsi(i) = .TRUE.
       ELSE
         x(1:3,j:j+1) = 0d0           !null values (unnecessary)
       END IF
     ELSE            !ELSE IF( l < 0 )THEN branching side
       sides(i) = -1   !side element is a branching side
       !j = -l          !position in side conns
       si => e%si(i)%p !point SI to the corresponding branching side
       ! assign present distorsion at side
       j = -l          !position
       IF( si%lnods(2*j-1) > 0 )THEN       !if order
         e%gamma(i) = +si%gamma(j)      !(+)
       ELSE                            !reversed order
         e%gamma(i) = -si%gamma(j)     !(-)
       END IF
       ! generate connectivities for enlarged B matrix
       k = nbn + 1                     !next position
       nbn = nbn + 2*si%nn-2             !update number of nodes
       lnb(k:nbn) = (/ ABS(si%lnods(1:2*j-2)), ABS(si%lnods(2*j+1:si%nn*2)) /) !+/-
     END IF
   END DO
   ! compute first and second fundamental forms

   CALL stra25(e%cartd,e%normd,e%ns,t,t3,tb,tns,e%a0,e%ci,x,sides,is,h, &
               hb,dt,e%area,e%gamma,e%lb,stra1,strag,fff)

   INCLUDE '../lbst/resv_p2.fpp'            !check thickness ratio


   ! Modify curvatures, substracting the original ones
   IF( shell ) stra1(4:6) = stra1(4:6) - e%stra0(1:3) + stra1(7:9)
   e%stra1 = stra1(1:6)    ! keep Mid Surface metric tensors
   stra2 = stra1(1:3)
   thnew = thick*e%lb                  !thickness at new configuration

   IF( bottom .OR. top) &      !bottom or top surfaces required
      CALL gentbs(is,e%lnods(1:4),bottom,top,coorb,coort,ifact,t3,thnew,e%normd(1:3,2,5:8),4)

   area1 = e%area(0)/2d0
   INCLUDE '../lbst/resv_p3.fpp'                !trougth the thickness integration

   stine(7:8) = 0d0                    !initializes shear
   DO i=1,4
     n = e%lside(i)                 ! element side
     IF(n > 0) THEN                ! element exist ?
       stine(7) = stine(7) + stint(4,n) * e%ns(1,i) + stint(6,n) * e%ns(2,i)
       stine(8) = stine(8) + stint(6,n) * e%ns(1,i) + stint(5,n) * e%ns(2,i)
     ELSE
       stine(7) = stine(7) + stint(4,iel) * e%ns(1,i) + stint(6,iel) * e%ns(2,i)
       stine(8) = stine(8) + stint(6,iel) * e%ns(1,i) + stint(5,iel) * e%ns(2,i)
     END IF
   END DO
   stine(7:8) = stine(7:8)/e%area(0)/2d0
   stint(1:10,iel) = stine                      !store N, M, Q & vMS  for post-processing

   stine(1:6) = stine(1:6)*e%area(0)            !Original area used (TLF)

   CALL bmem25(e%cartd(:,:,1:4),bmatx(:,1:4,1:3),bb,t(:,:,1:4),e%area(1:4),4)

   IF(shell) &
     CALL bfle25(e%normd,e%ci,e%ns,sides,is,solsi,bmatx(:,:,4:7),bbk,bbg,dt,t3,tb,e%area(0))
   

  IF ( branch ) THEN
    IF( nbn /= nbo )THEN
      DEALLOCATE( bbr,bbrk,rb )
      ALLOCATE( bbr(3,nbn,6),bbrk(3,nbn,4),rb(3,nbn) )
      nbo = nbn
    END IF
    CALL bbra25(bmatx,bbg,e%lb,sides,nbn,e%si,bbr,bbrk,e%lnods(5:12),dt,e%ns)
    CALL proma3(rb(1,1),bbr(1,1,1),stine(1),3*nbn,6,.TRUE.) !membrane + bending contributions

    ! stabilization contribution
    IF( stm > 0d0 )THEN          !membrane
      CALL lgst14(stra2,r1,r2,lb,'RESVPL',error)    !Hencky membrane (logarithmic) strains
      IF( secty == 12 .AND. ASSOCIATED(e%gausv) )THEN  !do not know how to handle section 13
        pstr = .TRUE.
        IF( MOD(nlayr,2) == 0 )THEN  !
          i = nlayr / 2
          strpl(3) = (e%gausv(3,i)+e%gausv(3,i+1))/2d0  !shear component only
        ELSE
          i = (nlayr+1) / 2
          strpl(3) = e%gausv(3,i)                       !shear component only
        END IF
      ELSE
        pstr = .FALSE.
      END IF
      DO k=1,4
        f2 = (fff(3,k) - e%stra1(3))                       !hourglass shear strain
        IF( pstr )THEN
          strpl(3) = stra2(3) - strpl(3)                   !elastic shear strain
          ratio = MAX(MIN((strpl(3)/e%stra1(3)),1d0),0.01) !elastic-plastic factor
          f2 = f2*ratio                                    !corrected stress
        END IF
        aux = e%area(k)*e%area(0)*stm*dm(1)         !area * C11 * Beta
        f2 = f2*aux                                 !hourglass stress
        CALL proma3(rb(1,1),bb(1,3,k),f2,12,1,.FALSE.)  !nodal contributions
      END DO
    END IF
    stint(9,iel) = f2                             !ratio for post-process
    IF( stbf > 0d0 )THEN          !bending
      aux = e%area(0)*stbf
      DO k=1,4
        sigma(1) = strag(k)*aux   !bending hourglass stress
        CALL proma3(rb(1,1),bbrk(:,:,k),sigma(1),3*nbn,1,.FALSE.)  !nodal contributions
      END DO
    END IF

    DO i=1,nbn            !for each node in the patch
      n = lnb(i)          !node
      !$ CALL omp_set_lock(lock1(MOD(n,MAXLOCK)))    !Start blocking acess in shared variables
      resid(1:3,n) = resid(1:3,n)+rb(1:3,i) !sum internal force
      !$ CALL omp_unset_lock(lock1(MOD(n,MAXLOCK)))  !End blocking acess in shared variables
    END DO
  ELSE
    CALL proma3(rl(1,1),bmatx(1:3,1:4,1:3),stine(1),12,3,.TRUE.)  !membrane contributions
    IF(shell)THEN
      rl(:,5:12) = 0d0
      CALL proma3(rl(1,1),bmatx(:,:,4:6),stine(4),36,3,.FALSE.) !bending contributions
    END IF

    ! stabilization contribution
    IF( stm > 0d0 )THEN          !membrane
      CALL lgst14(stra2,r1,r2,lb,'RESVPL',error)    !Hencky membrane (logarithmic) strains
      IF( secty == 12 .AND. ASSOCIATED(e%gausv) )THEN  !do not know how to handle section 13
        pstr = .TRUE.
        IF( MOD(nlayr,2) == 0 )THEN  !
          i = nlayr / 2
          strpl(3) = (e%gausv(3,i)+e%gausv(3,i+1))/2d0  !shear component only
        ELSE
          i = (nlayr+1) / 2
          strpl(3) = e%gausv(3,i)                       !shear component only
        END IF
      ELSE
        pstr = .FALSE.
      END IF
      DO k=1,4
        f2 = (fff(3,k) - e%stra1(3))                       !hourglass shear strain
        IF( pstr )THEN
          strpl(3) = stra2(3) - strpl(3)                   !elastic shear strain
          ratio = MAX(MIN((strpl(3)/e%stra1(3)),1d0),0.01) !elastic-plastic factor
          f2 = f2*ratio                                    !corrected stress
        END IF
        aux = e%area(k)*e%area(0)*stm*dm(4)         !area * C11 * Beta
        f2 = f2*aux                                 !hourglass stress
        CALL proma3(rl(1,1),bb(1,3,k),f2,12,1,.FALSE.)  !nodal contributions
      END DO
    END IF

    IF( stbf > 0d0 )THEN          !bending
      aux = e%area(0)*stbf

      DO k=1,4
        sigma(1) = strag(k)*aux   !bending hourglass stress
        CALL proma3(rl(1,1),bbk(:,:,k),sigma(1),36,1,.FALSE.)  !nodal contributions
      END DO
    END IF
    mti = .FALSE.
    DO i=4,1,-1 !invert order to use LNA and Y
      IF( solsi(i) )THEN
        j = nextn(i)       !i is the first node  j is the second node of the side
        t(:,1,1) = x(:,j)-x(:,i)       !side J-K
        CALL vecuni(3,t(1,1,1),zk)     !unit side vector
        ! bending moment times original side length / 2 / present thickness
        aux =-(stine(4)*e%ns(1,i)**2+stine(5)*e%ns(2,i)**2-2d0*stine(6)*e%ns(1,i)*e%ns(2,i)) &
              /e%cartd(4,2,i+4)/e%area(0)/thnew/2d0   !times the side length, divided by the thickness
        stres = x(:,fn(i))*aux            !forces normal to the face
        ! torsor moment  = forces applied on nodes J and K
        aux = ((stine(4)-stine(5))*e%ns(1,i)*e%ns(2,i) + stine(6)*(e%ns(1,i)**2-e%ns(2,2)**2))/e%area(0)
        lc = aux*t3                    !point loads leading to torsor moment
        !
        rl(:,i) = rl(:,i) - lc          !modify node J forces ????
        rl(:,j) = rl(:,j) + lc          !modify node K forces
        ! "S" forces = torsor moment * original length / 2 / present thickness
        aux = aux*e%cartd(4,2,i+4)/thnew/2d0
        stres = stres + aux*t(:,1,1)      !add to node forces
        CALL vecpro(t(1,1,1),x(1,fn(i)),t(1,2,1)) !both are unit vectors ==> t vector at the side
        t(:,2,1) = t(:,2,1)*thnew/2d0             !distance along normal
        yp(:,1) = x(:,i) - t(:,2,1)  !under i == 1
        yp(:,2) = x(:,i) + t(:,2,1)  !above i == 2
        yp(:,3) = x(:,j) + t(:,2,1)  !above j == 3
        yp(:,4) = x(:,j) - t(:,2,1)  !under j == 4
        yl = RESHAPE((/ -1d0,-1d0, 1d0,-1d0, 1d0, 1d0,-1d0, 1d0 /),(/2,4/))
        IF( mti )THEN !recompute LNA and Y
          n = -iffix(e%lnods(i))  !compute position of the i+1 node
          lna(1:2) = nndpd(2:3,n)
          n = -iffix(e%lnods(j))  !compute position of the i+2 node
          lna(3:4) = (/ nndpd(3,n),nndpd(2,n) /)
          y =  coora(:,lna(:))  ! nodes of the quadrilateral
        END IF
        DO n=1,4
          ! recompute proyections
          CALL projq4(yp(1,n),y(1,1),yl(1,n),yl(2,n),sh(1),.TRUE.)
          lb = fac(n)*stres     !equivalent force
          DO j=1,4
            i3 = lna(j)          !node
            !$ CALL omp_set_lock(lock1(MOD(i3,MAXLOCK)))    !Start blocking acess in shared variables
            resid(:,i3) = resid(:,i3) + sh(j)*lb
            !$ CALL omp_unset_lock(lock1(MOD(i3,MAXLOCK)))  !End blocking acess in shared variables
          END DO
        END DO
        mti = .TRUE.
      END IF
    END DO
    DO i=1,12             !for each node in the patch
      n = e%lnods(i)     !node
      !$ CALL omp_set_lock(lock1(MOD(n,MAXLOCK)))    !Start blocking acess in shared variables
      IF(n > 0) resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
      !$ CALL omp_unset_lock(lock1(MOD(n,MAXLOCK)))  !End blocking acess in shared variables
    END DO
  END IF
 e => e%next          !point to next node in the list
 END DO ELEMS
 IF( nbs > 0 ) DEALLOCATE( bbr,bbrk,rb )
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
 9999 CALL runen2('')
 END SUBROUTINE resv25
