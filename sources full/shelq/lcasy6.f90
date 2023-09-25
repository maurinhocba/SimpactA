 SUBROUTINE lcasy6(deriv,x,jacin,j0,t3,ang,locax)
 !***********************************************************************
 !
 !****this routine sets up the local cartesyan system for element
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 REAL (kind=8), INTENT(IN) :: deriv(:,:),x(:,:),ang
 REAL (kind=8), INTENT(OUT) :: t3(:),j0,jacin(:,:)
 INTEGER(kind=4), INTENT(IN) :: locax

 !                   local variables

 REAL (kind=8) t1(3),t2(3),ta(3),gcova(3,2),lt

  INTERFACE
    INCLUDE 'vecpro.h'
    INCLUDE 'vecuni.h'
  END INTERFACE
 !     evaluates the local cartesian system

 !         convective base
 gcova = MATMUL(x,deriv)
 !         t3 normal to the plane
 CALL vecpro(gcova(:,1),gcova(:,2),t3)
 !         normalizes t3   &   j0 = |t3| (jacobian)
 CALL vecuni(3,t3,j0)
 SELECT CASE(locax)
 CASE (1)
   lt = (t3(2)*t3(2)+t3(3)*t3(3)) !component in th Y-Z plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  Y-Z plane
     t2 = (/ -t3(3), 0d0, t3(1) /) !choose t2 orthogonal to global Y direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local y=t1 in the global YZ plane
     t1 = (/ 0d0, -t3(3), t3(2)  /)
     t2 = (/ lt, -t3(2)*t3(1), -t3(3)*t3(1) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (2)
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
 CASE (3)
   lt = (t3(1)*t3(1)+t3(2)*t3(2)) !component in th X-Y plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  X-Y plane
     t2 = (/ 0d0, t3(3), -t3(2) /) !choose t2 orthogonal to global X direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local x=t1 in the global xy plane
     t1 = (/ -t3(2), t3(1) , 0d0 /)
     t2 = (/ -t3(1)*t3(3), -t3(2)*t3(3), lt /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 END SELECT
 ta = t1
 t1 =  COS(ang)*t1 + SIN(ang)*t2
 t2 = -SIN(ang)*ta + COS(ang)*t2
 !                  evaluates inverse of jacobian matrix (transpose)
 jacin(1,1) =  DOT_PRODUCT(gcova(:,2),t2)/j0     ! phi^xita . t_x
 jacin(2,1) = -DOT_PRODUCT(gcova(:,2),t1)/j0     ! phi^xita . t_y
 jacin(1,2) = -DOT_PRODUCT(gcova(:,1),t2)/j0     ! phi^eta  . t_x
 jacin(2,2) =  DOT_PRODUCT(gcova(:,1),t1)/j0     ! phi^eta  . t_y

 RETURN
 END SUBROUTINE lcasy6
