 SUBROUTINE real20(nelem,heade,taile,iwrit,eulrf) !,eside)
 !=======================================================================
 ! Reallocates element set database after remeshing
 !=======================================================================
   USE npo_db, ONLY: oldlb,coora
   USE ctrl_db, ONLY: npoin
   IMPLICIT NONE
   !--- Dummy variables
   INTEGER(kind=4),INTENT(IN):: iwrit
   INTEGER(kind=4),INTENT(OUT):: nelem
   LOGICAL, INTENT(IN):: eulrf
   !INTEGER(kind=4)    :: eside(:,:)
   TYPE(ele20),POINTER:: heade, taile
   !!--- Functions
   !INTEGER(kind=4) chnode
   !--- Constants
   REAL(kind=8),PARAMETER:: cte=1d0/3d0,  d32=1.5d0
   !--- Local variables
   LOGICAL        :: plast,visco
   INTEGER(kind=4):: i, isec, mtype, ielem, nvare, lind(3), &
                     ipt, nvar, nstint, fyield, maxne, oldnn(1)
   REAL(kind=8)   :: lb(3), r1, r2, km, gm, jac, angle, press, delta
   REAL(kind=8), ALLOCATABLE :: strsg(:)
   TYPE(ele20),POINTER:: elm
   TYPE(ele20),POINTER:: posic,anter
   TYPE(section),POINTER:: sec     !auxiliar pointer

   !Initialization
   nstint = SIZE(heade%stint)
   maxne  = 0
   isec  = heade%matno              !material number is the same for all the set
   angle = heade%angle              !local angle is the same for all the set
   sec => psecs(isec)%p             ! pointer to material
   mtype = sec%mtbas%mtype          ! check material type for element
   fyield = sec%mtbas%matdef(3)     ! yield function
   plast = fyield > 1               ! plasticity
   visco = sec%mtbas%matdef(6) == 3 !Viscoplasticity
   km    = sec%mtbas%prope(4)       ! bulk modulus
   gm    = 2d0*sec%mtbas%prope(3)   ! twice the shear modulus

   IF( plast )THEN
     IF( mtype == 1 .AND. fyield == 2 )THEN
       nvare = 6
       IF( visco ) nvare = nvare + 1
     ELSE
       nvare = 5
     END IF
   ELSE
     nvare = 0
   END IF


   !------------------- DELETING ELEMENTS
   IF(r_elm_zone)THEN    !----> Delete distorted element in remesh zone
     NULLIFY(posic,anter) !initialize
     elm => heade ! point to head element
     DO
       IF(.NOT.ASSOCIATED(elm))EXIT
       posic => elm    ! element position
       elm => elm%next
       IF(posic%repla)THEN !delete element if repla = .TRUE.
         CALL del_ele20e(heade, taile, anter, posic)
       ELSE ! if not restore the original node label
         maxne = MAX(posic%numel,maxne) ! max element label number
         lind(1:3) = posic%lnods(1:3)
         DO i=1,3 ! return the original labels in conectivities
           posic%lnods(i  ) = oldlb(lind(i))
           posic%lnods(i+3) = 0
         END DO
         anter => posic ! point to 'next' previous element
       END IF
     END DO
     NULLIFY(posic,anter) !make null the auxiliar pointers
     maxne = ((maxne/1000)+1)*1000  ! same that FF

   ELSE                   !----> Delete ALL set of elements
     CALL delete_ele20e(heade,taile)
   END IF

   !Copy data to new set list of elements
   IF (iwrit == 1) WRITE(lures,"(//5x,' Element   Property      Node Numbers')",ERR=9999)
   nelem = ne_new               !number of new elements
   nvar = SIZE(sg_new,DIM=1)    !number of transfered variables
   ALLOCATE(strsg(nvar))

   !---------------------- ADDING NEW ELEMENTS
   DO ielem=1,nelem   !Loop over new set of elements
     CALL new_ele20e(elm)   !Create a new element

     elm%numel = maxne + ielem      !element label
     elm%matno = isec               !associated material
     elm%angle = angle              !

     !check de material for set variables to transfer
     IF ( nvare > 0) THEN
       ALLOCATE(elm%gausv(nvare))
       elm%gausv = 0d0
     END IF
     ! deal with connectivities
     IF(r_elm_zone)THEN      ! zone remeshing
       DO i=1,3  !for each node
         ipt   = l_new(i,ielem)            !new node
         oldnn = MAXVAL(nodlb(2,:),MASK=nodlb(1,:) == (ipt+maxnn))
         IF(oldnn(1) > 0)THEN        ! is a frontier node (exists in the original mesh)
           elm%lnods(i) = oldnn(1)     !associate with original label
         ELSE                        ! is a new node
           elm%lnods(i) = ipt + maxnn  !transfer and modify conectivities numeration
         END IF
       END DO
     ELSE                    ! full remeshing
       DO i=1,3
         ipt = l_new(i,ielem)
         elm%lnods(i) = ipt + maxnn     !transfer and modify conectivities numeration
       END DO
     END IF
     !Loading information into new Gauss Points
     strsg(1:nvar) = sg_new(1:nvar,ielem) !Gauss transfered values

     IF ( mtype == 1 .AND. fyield == 2 ) THEN  !interpolation to Gauss-Point
       IF(eulrf)THEN
         !strsg(1:4) deviatoric Finger tensor
         delta = strsg(5)     !volumetric part Finger tensor
         elm%gausv(1:nstint) = strsg(1:nstint) + (/ delta,delta,0d0,delta /) !Finger tensor
         elm%gausv(6) = strsg(6)     !eq. plastic strain
         IF( visco ) elm%gausv(5) = strsg(7)     !plastic strain rate
         jac   = elm%gausv(1) + elm%gausv(2) + elm%gausv(4)
         press = km*(3d0-jac)/2d0    !volum. mean stress
         elm%stint(1:nstint) = (gm/2d0)*strsg(1:nstint) + & !Kirchhoff Stress
                               (/ press,press,0d0,press /)  !
       ELSE
         !strsg(1:4) is deviatoric Kirchhoff stress
         press = strsg(5)     !mean press
         jac   = SQRT(d32*(strsg(1)**2+strsg(2)**2+strsg(4)**2)+2d0*strsg(3)**2) !vMises
         strsg(1:4) = strsg(1:4)*strsg(6)/jac  !correct to vMises stress
         elm%stint(1:nstint) = strsg(1:nstint) + (/ press,press,0d0,press /) !Kirchhoff Stress

         CALL edst20 (strsg(1),r1,r2,lb(1)) !,ierr) !eigen-decomposition of Kirchhoff stress
                                                    !tensor
         jac = (EXP(press/km))**cte   !jacobian^(1/3)
         DO i=1,3
           lb(i) = jac*EXP(lb(i)/gm)       !ppal stretches
         END DO

         elm%gausv(1) = lb(1)*r1*r1 + lb(2)*r2*r2   !
         elm%gausv(2) = (lb(1)-lb(2))*r1*r2         !
         elm%gausv(3) = (lb(1)-lb(2))*r1*r2         ! inverse plastic gradient
         elm%gausv(4) = lb(1)*r2*r2 + lb(2)*r1*r1   !
         elm%gausv(5) = lb(3)                       !
         elm%gausv(6) = strsg(7)                 ! effective plastic strain
       END IF
       IF( elm%gausv(6) < 0d0 )elm%gausv(6) = 0d0 !not negative

     ELSE IF( plast ) THEN !other materials and formulation
       !strsg(1:4) is deviatoric Kirchhoff stress          !interpolation to Gauss-Point
       press = strsg(5)                   !mean press
       jac   = SQRT(d32*(strsg(1)**2+strsg(2)**2+strsg(4)**2)+2d0*strsg(3)**2) !vMises
       strsg(1:4) = strsg(1:4)*strsg(6)/jac  !correct to vMises stress
       elm%stint(1:nstint) = strsg(1:nstint) + (/ press,press,0d0,press /) !Kirchhoff Stress

       ! Ep = -Ee
       jac = press / km * cte                    !volumetric strain / 3
       elm%gausv(1) = -strsg(1)/gm - jac
       elm%gausv(2) = -strsg(2)/gm - jac
       elm%gausv(3) = -strsg(3)*2d0 / gm
       elm%gausv(4) = -strsg(4)/gm - jac
       elm%gausv(5) =  strsg(7)                   ! effective plastic strain
       IF( elm%gausv(5) < 0d0 )elm%gausv(5) = 0d0 !not negative
     END IF

     CALL add_ele20e(elm,heade,taile)  !Add element to list
   END DO

   elm => heade !point to first element in the list
   ielem = 0
   DO
     IF(.NOT.ASSOCIATED(elm))EXIT
     ielem = ielem + 1
     !!IF (iwrit == 1) WRITE(lures,"(6x,i5,i9,6x,10i5)",ERR=9999) ielem,elm%matno,elm%lnods(1:3) !echo element
     !!!Original nodes are changed for internal association
     !!DO i=1,3
     !!  elm%lnods(i) = chnode(elm%lnods(i))
     !!  IF(elm%lnods(i) == 0) WRITE(*,*) 'label elpos nodpos: ', elm%numel, ielem, i
     !!END DO
     elm => elm%next
   END DO
   nelem = ielem ! updates the element number

 RETURN
  9999 CALL runen2('')
 END SUBROUTINE real20


 SUBROUTINE edst20 (stran,r1,r2,lb) !,ierr)
  ! compute eigen-decomposition for Kirchhoff stress tensor
  IMPLICIT NONE
  ! dummy arguments
  REAL(kind=8), INTENT (OUT) :: lb(3)       !eigenvalues of Tau
  REAL(kind=8), INTENT (IN)  :: stran(4)    !(IN) deviatoric stress
  REAL(kind=8), INTENT (OUT) :: r1,r2       !components of first eigenvector
  !INTEGER (kind=4), INTENT(OUT) :: ierr

  ! local variables
  REAL (kind=8) :: c,r,d

  !compute eigenvalues
  c = (stran(1)+stran(2))/2d0  !center of circle
  d = (stran(1)-stran(2))/2d0  !half the difference
  r = SQRT(d**2+stran(3)**2)   !circle radius
  !Compute principal deviatoric stresses stresses
  lb(1) = c+r                                 !first (maximum) eigenvalue
  lb(2) = c-r                                 !second (minimum) eigenvalue
  lb(3) = stran(4)
  !compute eigenvectors
  IF( r > 1d-15)THEN                       !check
    c = stran(1) - lb(1)                   !first diagonal element
    r = -stran(3)/c                        ! v = (r,1)
    c = SQRT(r**2+1d0)                     ! eigenvector length
    r1 = r/c                               !first component of eigevector
    r2 = 1d0/c                             !first component of eigevector
  ELSE
    r1 = 1d0                               !local direction is the vector
    r2 = 0d0
  END IF
  RETURN
 END SUBROUTINE edst20
