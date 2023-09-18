 SUBROUTINE gaus16(heade,coord,istop,gauss,angdf,ngaus,locax,quad,shell)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for solid element 16 PRISM
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(OUT) :: istop
 INTEGER (kind=4), INTENT(IN) :: ngaus,locax
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf(:)      !default angles Gl_X - Ort_Axes
 TYPE (ele16), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if gauss constants computed
 LOGICAL, INTENT (IN) :: quad         !TRUE for quadratic version
 LOGICAL, INTENT (IN) :: shell        !TRUE for assumed transverse shear

 INTEGER (kind=4), PARAMETER :: ngaud = 2   !two gauss point in zeta direction
 INTEGER (kind=4) i,ielem
 REAL (kind=8) :: shape(nnode,ngaus), &
                  deriv(nnode,3,ngaus),weigp(ngaud),posgp(ngaud),gpw(ngaus)
 REAL (kind=8) :: ang(3),factor,rm(3,3),lb(3,3)
 REAL (kind=8), ALLOCATABLE :: x(:,:)
 REAL (kind=8), PARAMETER :: a=0.333333333333333d0       !1/3

 TYPE (ele16), POINTER :: e           !pointer to an element
 LOGICAL :: use_gang, use_ang         !global non-zero angles,  a local system exist
 LOGICAL :: sides(3)                  !side element existence for quad approach

 IF( gauss ) RETURN
 !     ---------------------------------------------------------

 CALL gaussq (ngaud ,posgp ,weigp ) ! get standard  gauss points in local coordinates and weigths
 !  get memory for nodal coordinates
 IF( quad ) THEN              ! for quad approach
   ALLOCATE(x(3,12))      !extended patch
 ELSE                         !
   ALLOCATE(x(3,nnode))   !standard element
 END IF
 DO i=1,ngaud
   CALL shape4(nnode,shape(1,i),deriv(1,1,i),a,a,posgp(i))
 END DO
 gpw = 0.5d0          !1/2

 !     ---------------------------------------------------------

 factor = ASIN(1d0)/90d0         !pi/180
 ang = angdf*factor              !default angle in radians
 use_gang = ANY(angdf /= 0d0)    !only if not all angles are null
 CALL inrotm(ang,rm(1,1))        !compute global relative rotation matrix

 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit

   x(:,1:nnode) = coord(:,e%lnods(1:nnode))     !get original coordinates (element nodes)
   IF( quad )THEN   !for quadratic approach, get extra nodes
     DO i=1,3        !get original coordinates (extra patch nodes)
       sides(i) = ( e%lnods(i+6) /= 0 )         !side element exist
       IF( sides(i) )THEN
         x(:,i+6:i+9:3) = coord(:,e%lnods(i+6:i+9:3))
       ELSE
         x(:,i+6:i+9:3) = 0d0                   !may be unnecessary
       END IF
     END DO
   END IF

   IF( ANY(e%angle /= 0d0))THEN          !IF local angle is non zero
     ang = e%angle*factor                  ! in radians
     CALL inrotm(ang,lb(1,1))                !local relative rotation matriz
     lb = MATMUL(rm,lb)                         ! rotation matrix
     use_ang = .TRUE.
   ELSE
     lb = rm                                 !use global matrix
     use_ang = use_gang
   END IF
   IF( locax /= 0 )CALL lcsy16(lb,x,locax,use_ang)   !compute local shell system

   ielem = e%numel
   ! compute element local systems, nodal functions and derivatives
   CALL jacob16(quad,shell,deriv,e%cartd,e%dvol,x,ielem,istop, &
                e%cdq,lb,sides,e%nfdas,e%jacin)
   e%dvol= e%dvol*gpw(1)    !both points have the same weight
   e%stint = 0d0            !Initializes stresses

   e => e%next          !point to next element
 END DO
 gauss = .TRUE.

 RETURN
 END SUBROUTINE gaus16
