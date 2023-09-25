 SUBROUTINE lcsy04(lb,x,locax,ang)
 !***********************************************************************
 !
 !****this routine sets up the local cartesyan system for element PRISM
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   dummy arguments

 REAL (kind=8), INTENT(IN) :: x(:,:),ang    !nodal coordinates & local angle
 REAL (kind=8), INTENT(OUT) :: lb(3,3)      !local system
 INTEGER(kind=4), INTENT(IN) :: locax       !local normal axiz


 !                   local variables
 INTEGER(kind=4) :: i
 REAL (kind=8) :: t(3,3),ta(3),lt,cosa,sina,y(3,4)

  INTERFACE
    INCLUDE 'vecpro.h'
    INCLUDE 'vecuni.h'
  END INTERFACE

 !     evaluates the local cartesian system

 DO i=1,4  !mid surface coordinates
   y(:,i) = (x(:,i) + x(:,i+4))/2d0
 END DO
 t(:,1) = - y(:,1) + y(:,2) + y(:,3) - y(:,4)           !xita dir
 t(:,2) = - y(:,1) - y(:,2) + y(:,3) + y(:,4)           !eta dir

 !*** evaluate the cross product => plane normal

 CALL vecpro(t(:,1),t(:,2),t(:,3))                    !normal at center
 CALL vecuni(3,t(:,3),lt)                             !set unit normal vector

 SELECT CASE(locax)
 CASE (1)     !IF X is the prefered normal vector
   lt = (t(2,3)*t(2,3)+t(3,3)*t(3,3)) !component in th Y-Z plane
   IF( lt  < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  Y-Z plane
     t(:,2) = (/ -t(3,3), 0d0, t(1,3) /) !choose t2 orthogonal to global Y direction
     CALL vecuni(3,t(:,2),lt)
     CALL vecpro(t(:,2),t(:,3),t(:,1))
   ELSE       !         SELECT local y=t(:,1) in the global YZ plane
     t(:,1) = (/ 0d0, -t(3,3), t(2,3)  /)
     t(:,2) = (/ lt, -t(2,3)*t(1,3), -t(3,3)*t(1,3) /)
     CALL vecuni(3,t(:,1),lt)   !     normalizes t1 & t2
     CALL vecuni(3,t(:,2),lt)
   END IF
 CASE (2)     !IF Y is the prefered normal vector
   lt = (t(3,3)*t(3,3)+t(1,3)*t(1,3)) !component in th Z-X plane
   IF( lt  < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  Z-X plane
     t(:,2) = (/ 0d0, t(3,3), -t(2,3) /) !choose t(:,2) orthogonal to global X direction
     CALL vecuni(3,t(:,2),lt)
     CALL vecpro(t(:,2),t(:,3),t(:,1))
   ELSE       !         SELECT local z=t(:,1) in the global ZX plane
     t(:,1) = (/ t(3,3), 0d0, -t(1,3)  /)
     t(:,2) = (/  -t(1,3)*t(2,3), lt, -t(3,3)*t(2,3) /)
     CALL vecuni(3,t(:,1),lt)   !     normalizes t1 & t2
     CALL vecuni(3,t(:,2),lt)
   END IF
 CASE (3)     !IF Z is the prefered normal vector
   lt = (t(1,3)*t(1,3)+t(2,3)*t(2,3)) !component in th X-Y plane
   IF( lt  < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  X-Y plane
     t(:,2) = (/ 0d0, t(3,3), -t(2,3) /) !choose t(:,2) orthogonal to global X direction
     CALL vecuni(3,t(:,2),lt)
     CALL vecpro(t(:,2),t(:,3),t(:,1))
   ELSE       !         SELECT local x=t(:,1) in the global xy plane
     t(:,1) = (/ -t(2,3), t(1,3) , 0d0 /)
     t(:,2) = (/ -t(1,3)*t(3,3), -t(2,3)*t(3,3), lt /)
     CALL vecuni(3,t(:,1),lt)   !     normalizes t(:,1) & t(:,2)
     CALL vecuni(3,t(:,2),lt)
   END IF
 END SELECT
 IF( ang /= 0d0 )THEN         !If lb is non trivial
   ! compute angle between local system t1-t2  and L1
   cosa = COS(ang)
   sina = SIN(ang)
   ! rotate local system t1-t2 to best fit L1-L2
   ta = t(:,1)
   t(:,1) =  cosa*t(:,1) + sina*t(:,2)
   t(:,2) = -sina*ta     + cosa*t(:,2)
 END IF
 lb = t       !computed local system at element center
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lcsy04
