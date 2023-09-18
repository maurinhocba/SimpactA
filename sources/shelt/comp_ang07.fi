 SUBROUTINE comp_ang07(x,t,lb,nangl,ang,locax,jas)
 !***********************************************************************
 !
 !****this routine compute the angle (cos,sin) between local cartesyan system
 !    and nodal system
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 INTEGER (kind=4), INTENT(IN) :: locax

 REAL (kind=8), INTENT(IN) :: x(:,:),lb(:,:,:),t(:,:),ang

 REAL (kind=8), INTENT(OUT) :: nangl(:,:)
 REAL (kind=8), POINTER, INTENT(OUT) :: jas(:,:)

 !                   local variables

 INTEGER (kind=4) :: n
 REAL (kind=8) :: t1(3),t2(3),t3(3),j0,lt,cb,gcova(3,2),l1(3),l2(3)

 INTERFACE
    INCLUDE 'vecpro.h'
    INCLUDE 'vecuni.h'
 END INTERFACE

 gcova(:,1) = x(:,2) - x(:,1) !  convective base
 gcova(:,2) = x(:,3) - x(:,1) !  convective base
 CALL vecpro(gcova(:,1),gcova(:,2),t3)   !  t3 normal to the plane
 CALL vecuni(3,t3,j0)                    !  normalizes t3   &   j0 = |t3| (jacobian)
 SELECT CASE(locax)           !according to user preference
 CASE (1)                !use X as the default normal
   lt = (t3(2)*t3(2)+t3(3)*t3(3))  !component in th Y-Z plane
   IF( lt  < 1.0d-5) THEN          !If t3 is almost orthogonal to  Y-Z plane
     t2 = (/ -t3(3), 0d0, t3(1) /) !choose t2 orthogonal to global Y direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local y=t1 in the global YZ plane
     t1 = (/ 0d0, -t3(3), t3(2)  /)
     t2 = (/ lt, -t3(2)*t3(1), -t3(3)*t3(1) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (2)                !use Y as the default normal
   lt = (t3(3)*t3(3)+t3(1)*t3(1)) !component in th Z-X plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  Z-Y plane
     t2 = (/ t3(2), -t3(1), 0d0 /) !choose t2 orthogonal to global Z direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local z=t1 in the global ZX plane
     t1 = (/ t3(3), 0d0, -t3(1)  /)
     t2 = (/  -t3(1)*t3(2), lt, -t3(3)*t3(2) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (3)                !use Z as the default normal
   lt = (t3(1)*t3(1)+t3(2)*t3(2)) !component in th X-Y plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  X-Y plane
     t2 = (/ 0d0, t3(3), -t3(2) /) !choose t2 orthogonal to global X direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local x=t1 in the global XY plane
     t1 = (/ -t3(2), t3(1) , 0d0 /)
     t2 = (/ -t3(1)*t3(3), -t3(2)*t3(3), lt /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 END SELECT

 l1 =  COS(ang)*t1 + SIN(ang)*t2  !t_x
 l2 =  COS(ang)*t2 - SIN(ang)*t1  !t_y
 DO n=1,3        !for each mid-side point
   !       compute the local cartesian system
   !local systems at nodes
   cb = DOT_PRODUCT(t3,t(:,n))                !COS beta
   nangl(1,n) = DOT_PRODUCT(lb(:,1,n),l1)/cb  !COS alpha
   nangl(2,n) = DOT_PRODUCT(lb(:,2,n),l1)/cb  !SIN alpha
   !local systems at mid-side points
   !local derivative along sides at shear mid-points
   SELECT CASE (n)
   CASE (1)               !(-xita+eta) (node 4)
     jas(1,1) =  -DOT_PRODUCT(gcova(:,1),l1) +DOT_PRODUCT(gcova(:,2),l1)
     jas(2,1) =  -DOT_PRODUCT(gcova(:,1),l2) +DOT_PRODUCT(gcova(:,2),l2)
   CASE (2)               !-eta  (node 5)
     jas(1,2) =  -DOT_PRODUCT(gcova(:,2),l1)
     jas(2,2) =  -DOT_PRODUCT(gcova(:,2),l2)
   CASE (3)               !xita (node 6)
     jas(1,3) =  DOT_PRODUCT(gcova(:,1),l1)
     jas(2,3) =  DOT_PRODUCT(gcova(:,1),l2)
   END SELECT
 END DO

 RETURN
 END SUBROUTINE comp_ang07
