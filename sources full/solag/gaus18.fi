 SUBROUTINE gaus18(heade,coord,istop,gauss,angdf,ngaus,shell,gpc,check,locax)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 18
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(OUT) :: istop
 INTEGER (kind=4), INTENT(IN) :: ngaus,locax
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf(:)      !default angles Gl_X - Ort_Axes
 REAL (kind=8), INTENT(IN OUT) :: gpc(:,:)  !Gauss-point coordinates
 TYPE (ele18), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists
 LOGICAL, INTENT (IN OUT) :: check    !TRUE if connectivities order must be checked
 LOGICAL, INTENT (IN) :: shell        !TRUE to consider asummed shear strains

 INTEGER (kind=4) k,ielem
 REAL (kind=8) :: x(3,nnode),cd(nnode,3),shape(nnode,ngaus),dvol, &
                  deriv(nnode,3,ngaus)

 REAL (kind=8) :: ang(3),factor,rm(3,3),lb(3,3),gpw(ngaus),t(3)
 TYPE (ele18), POINTER :: e
 LOGICAL :: use_gang, use_ang

 !     ---------------------------------------------------------

 IF( gauss ) RETURN

 !       gauss points in local coordinates and weigths
 IF( ngaus == 4)THEN
   gpc = RESHAPE((/  sq3i, -sq3i, -sq3i,  -sq3i,  sq3i, -sq3i, &
                    -sq3i, -sq3i,  sq3i,   sq3i,  sq3i,  sq3i /),(/ 3,ngaus /))
   gpw(1:ngaus) = 2d0
 ELSE
   gpc = RESHAPE((/ -sq3i, -sq3i, -sq3i,   sq3i, -sq3i, -sq3i, &
                    -sq3i,  sq3i, -sq3i,   sq3i,  sq3i, -sq3i, &
                    -sq3i, -sq3i,  sq3i,   sq3i, -sq3i,  sq3i, &
                    -sq3i,  sq3i,  sq3i,   sq3i,  sq3i,  sq3i /),(/ 3,ngaus /))
   gpw(1:ngaus) = 1d0

 END IF
 !       gauss points shape and derivatives of nodal functions
 DO k=1,ngaus
   CALL shape5(deriv(1,1,k),shape(1,1),gpc(1,k),gpc(2,k),gpc(3,k),nnode)
 END DO


 !     ---------------------------------------------------------

 factor = ASIN(1d0)/90d0         !pi/180
 ang = angdf*factor
 use_gang = ANY(angdf /= 0d0)
 CALL inrotm(ang,rm(1,1))

 e => heade       !point to first element

 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   ielem = e%numel

   x = coord(:,e%lnods)     !get original coordinates

   IF( ANY(e%angle /= 0d0))THEN   !if local angles defined for the element
     ang = e%angle*factor
     CALL inrotm(ang,lb(1,1))
     lb = MATMUL(rm,lb)
     use_ang = .TRUE.
   ELSE
     lb = rm
     use_ang = use_gang
   END IF
   IF( shell )THEN  !for shell version
     IF( check ) THEN
       SELECT CASE (locax)
       CASE(1)
         t = (/ 1d0, 0d0, 0d0 /)
       CASE(2)
         t = (/ 0d0, 1d0, 0d0 /)
       CASE(3)
         t = (/ 0d0, 0d0, 1d0 /)
       CASE DEFAULT
         t = lb(:,3)
       END SELECT
       CALL check_orientation18(e%lnods,x,t) !lb(:,3))
     END IF
     CALL lcas18(deriv,lb,e%cartd,e%dvol,x,ngaus,use_ang,locax,e%nfdas,e%jacin)
     e%dvol= e%dvol*gpw(1)                !all GP have the same weigth
   ELSE             !standard version
     IF(locax == 0)THEN
       ! compute element local systems, nodal functions and derivatives
       DO k=1,ngaus
         CALL jacob5(cd(1,1),deriv(1,1,k),dvol,x(1,1),nnode,ielem,istop)
         e%cartd(:,:,k) = MATMUL(cd,lb)
         e%dvol(k)= dvol*gpw(k)
       END DO
     ELSE
       CALL lcas18(deriv,lb,e%cartd,e%dvol,x,ngaus,use_ang,locax)
       e%dvol= e%dvol*gpw(1)                !all GP have the same weigth
     END IF
   END IF
   e%stint = 0d0   !Initializes stresses

   e => e%next          !point to next element
 END DO
 gauss = .TRUE.

 RETURN
 END SUBROUTINE gaus18
 SUBROUTINE check_orientation18(lnods,x,t)
 ! check element connectivities
 IMPLICIT NONE
 !dummy arguments
 INTEGER(kind=4), INTENT(IN OUT) :: lnods(8) !connectivities
 REAL(kind=8), INTENT(IN OUT) :: x(3,8)      !coordinates
 REAL(kind=8), INTENT(IN) :: t(3)            !shell normal
 !local variables
 REAL(kind=8) :: t1(3),t2(3),t3(3),ca,z(3),y(3,8)
 INTEGER(kind=4) :: ln(8)
 INTERFACE
   INCLUDE 'vecuni.h'
   INCLUDE 'vecpro.h'
 END INTERFACE

 z = t                ! normal shell direction
 CALL vecuni(3,z,ca)  !unit shell normal
 t1 =  -x(:,1)+x(:,2)+x(:,3)-x(:,4)   !x_xita
 t2 =  -x(:,1)-x(:,2)+x(:,3)+x(:,4)   !x_eta
 CALL vecpro(t1,t2,t3)                !x_xita x x_eta
 CALL vecuni(3,t3,ca)                 !element normal
 ca = ABS( DOT_PRODUCT(z,t3) )        !COS between normals
 IF( ca > 0.866 )RETURN !O.K   angle is less than 30 degrees
 ln = lnods                           !assign connectivities to auxiliar array
 y = x                                !assign coordinates to auxiliar array
 t1 = x(:,2) - x(:,1)                 !vector along first arista
 CALL vecuni(3,t1,ca)                 !unit vector along first arista
 ca = ABS( DOT_PRODUCT(z,t1) )        !COS between shell normal and first arista
 IF( ca > 0.866 )THEN                   !if angle is less than 30 degrees
   lnods = ln((/ 1,4,8,5, 2,3,7,6 /))   !swap connectivities
   x     = y(:,(/ 1,4,8,5, 2,3,7,6 /))  !swap coordinates
   WRITE(58,"(8i6)")lnods
   RETURN !O.K
 END IF
 t1 = x(:,4) - x(:,1)                 !vector along last arista
 CALL vecuni(3,t1,ca)                 !unit vector along last arista
 ca = ABS( DOT_PRODUCT(z,t1) )        !COS between shell normal and last arista
 IF( ca > 0.866 )THEN                   !if angle is less than 30 degrees
   lnods = ln((/ 1,5,6,2, 4,8,7,3 /))   !swap connectivities
   x     = y(:,(/ 1,5,6,2, 4,8,7,3 /))  !swap coordinates
   WRITE(58,"(8i6)")lnods
   RETURN !O.K
 END IF
 WRITE(58,"(8i6,'WARNING: ORIENTATION COULD NOT BE CORRECTED')")lnods

 RETURN
 END SUBROUTINE check_orientation18
