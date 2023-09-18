 SUBROUTINE lcas18(deriv,lb,cd,dvol,x,ngaus,use_ang,locax,nfdas,jacin)
 !***********************************************************************
 !
 !****this routine sets up the local cartesyan system for element SOLAG
 !    version with assumed transverse shear for shells
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 REAL (kind=8), INTENT(IN) :: deriv(:,:,:),x(:,:),lb(3,3)
 REAL (kind=8), INTENT(OUT) :: cd(:,:,:),dvol(:)
 REAL (kind=8), OPTIONAL, INTENT(OUT) :: jacin(:,:,:),nfdas(:,:)
 INTEGER(kind=4), INTENT(IN) :: ngaus,locax
 LOGICAL, INTENT(IN) :: use_ang


 !                   local variables
 INTEGER(kind=4) :: g,i
 REAL (kind=8) :: t(3,3),ta(3),gcova(3,3),lt,jaci(3,3),j0,der(8,3),shape(8),cosa,sina
 REAL (kind=8), PARAMETER :: gpa(2,4) = (/ -1d0, 0d0,  0d0, -1d0,  1d0, 0d0,  0d0, 1d0 /)

  INTERFACE
    INCLUDE 'vecpro.h'
    INCLUDE 'vecuni.h'
    ! INCLUDE 'angeul.h' !to use in debug
  END INTERFACE
 !     evaluates the local cartesian system

 DO g=1,ngaus
 !         convective base
   gcova = MATMUL(x,deriv(:,:,g))  !in-plane covariant base
   !         t(:,3) normal to the plane
   CALL vecpro(gcova(:,1),gcova(:,2),t(:,3))
   !         normalizes t(:,3)   &   j0 = |t(:,3)| (jacobian)
   CALL vecuni(3,t(:,3),j0)

   SELECT CASE(locax)
   CASE (1)
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
   CASE (2)
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
   CASE (3)
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
   IF( use_ang )THEN
     cosa = DOT_PRODUCT(t(:,1),lb(:,1))
     sina = DOT_PRODUCT(t(:,2),lb(:,1))
     lt = SQRT(cosa**2+sina**2)
     cosa = cosa/lt
     sina = sina/lt
     ta = t(:,1)
     t(:,1) =  cosa*t(:,1) + sina*t(:,2)
     t(:,2) = -sina*ta     + cosa*t(:,2)
   END IF
   !CALL angeul(t,ta)  !to use in Debug
   dvol(g) =  gcova(1,1)*gcova(2,2)*gcova(3,3) +gcova(1,3)*gcova(2,1)*gcova(3,2)  &
             +gcova(3,1)*gcova(1,2)*gcova(2,3) -gcova(3,1)*gcova(2,2)*gcova(1,3)  &
             -gcova(3,3)*gcova(1,2)*gcova(2,1) -gcova(1,1)*gcova(2,3)*gcova(3,2)
   IF(dvol(g) <= 0) THEN
     !WRITE(*,600,ERR=9999) ielem
     !WRITE(*,900,ERR=9999)
     !WRITE(*,910,ERR=9999) (inode,elcod(1:ndime,inode),inode=1,nnode)
     !WRITE(lures,600,ERR=9999) ielem
     !WRITE(lures,900,ERR=9999)
     WRITE(lures,"(5x,i5,3e15.8)",ERR=9999) (i,x(:,i),i=1,8)
     !iSTOP = 1
     STOP 'null or negative jacobian'
     RETURN
   END IF
   gcova = MATMUL(TRANSPOSE(t),gcova)      !rotate to the local system
   ! compute inverse jacobian
   jaci(1,1) =  (gcova(2,2)*gcova(3,3)-gcova(2,3)*gcova(3,2))/dvol(g)
   jaci(2,1) = -(gcova(2,1)*gcova(3,3)-gcova(3,1)*gcova(2,3))/dvol(g)
   jaci(3,1) =  (gcova(2,1)*gcova(3,2)-gcova(2,2)*gcova(3,1))/dvol(g)
   jaci(1,2) = -(gcova(1,2)*gcova(3,3)-gcova(1,3)*gcova(3,2))/dvol(g)
   jaci(2,2) =  (gcova(1,1)*gcova(3,3)-gcova(3,1)*gcova(1,3))/dvol(g)
   jaci(3,2) = -(gcova(1,1)*gcova(3,2)-gcova(1,2)*gcova(3,1))/dvol(g)
   jaci(1,3) =  (gcova(1,2)*gcova(2,3)-gcova(1,3)*gcova(2,2))/dvol(g)
   jaci(2,3) = -(gcova(1,1)*gcova(2,3)-gcova(1,3)*gcova(2,1))/dvol(g)
   jaci(3,3) =  (gcova(1,1)*gcova(2,2)-gcova(2,1)*gcova(1,2))/dvol(g)

   !     calculate cartesian derivatives
   cd(:,:,g) = MATMUL(deriv(:,:,g),jaci)
   !der = MATMUL(deriv(:,:,g),jaci)
   !cd(:,1,g) = der(:,1)*t(1,1) + der(:,2)*t(2,1) + der(:,3)*t(3,1)
   !cd(:,2,g) = der(:,1)*t(1,2) + der(:,2)*t(2,2) + der(:,3)*t(3,2)
   !cd(:,3,g) = der(:,1)*t(1,3) + der(:,2)*t(2,3) + der(:,3)*t(3,3)
   !
   IF( PRESENT(jacin) ) jacin(:,:,g) = jaci(1:2,1:2)

 END DO

 IF( .NOT.PRESENT(jacin) .OR. .NOT.PRESENT(nfdas)) RETURN

 DO g=1,4   !shear assumed strain points
   CALL shape5(der,shape,gpa(1,g),gpa(2,g),0d0,8) !shape function derivatives

   gcova = MATMUL(x,der)  ! covariant base
   !         t(:,3) normal to the plane
   CALL vecpro(gcova(:,1),gcova(:,2),t(:,3))
   !         normalizes t(:,3)   &   j0 = |t(:,3)| (jacobian)
   CALL vecuni(3,t(:,3),j0)
   SELECT CASE(locax)
   CASE (1)
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
   CASE (2)
     lt = (t(3,3)*t(3,3)+t(1,3)*t(1,3)) !component in th Z-X plane
     IF( lt  < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  Z-Y plane
       t(:,2) = (/ t(2,3), -t(1,3), 0d0 /) !choose t2 orthogonal to global Z direction
       CALL vecuni(3,t(:,2),lt)
       CALL vecpro(t(:,2),t(:,3),t(:,1))
     ELSE       !         SELECT local z=t(:,1) in the global ZX plane
       t(:,1) = (/ t(3,3), 0d0, -t(1,3)  /)
       t(:,2) = (/  -t(1,3)*t(2,3), lt, -t(3,3)*t(2,3) /)
       CALL vecuni(3,t(:,1),lt)   !     normalizes t1 & t2
       CALL vecuni(3,t(:,2),lt)
     END IF
   CASE (3)
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
   IF( use_ang )THEN
     cosa = DOT_PRODUCT(t(:,1),lb(:,1))
     sina = DOT_PRODUCT(t(:,2),lb(:,1))
     lt = SQRT(cosa**2+sina**2)
     cosa = cosa/lt
     sina = sina/lt
     ta = t(:,1)
     t(:,1) =  cosa*t(:,1) + sina*t(:,2)
     t(:,2) = -sina*ta     + cosa*t(:,2)
   END IF
   !WRITE(55,"(3e21.10)")t
   gcova = MATMUL(TRANSPOSE(t),gcova)      !rotate to the local system

   j0 =  gcova(1,1)*gcova(2,2)*gcova(3,3) +gcova(1,3)*gcova(2,1)*gcova(3,2)  &
        +gcova(3,1)*gcova(1,2)*gcova(2,3) -gcova(3,1)*gcova(2,2)*gcova(1,3)  &
        -gcova(3,3)*gcova(1,2)*gcova(2,1) -gcova(1,1)*gcova(2,3)*gcova(3,2)

   IF( j0   <= 0d0 ) THEN
     !WRITE(*,600,ERR=9999) ielem
     !WRITE(*,900,ERR=9999)
     !WRITE(*,910,ERR=9999) (inode,elcod(1:ndime,inode),inode=1,nnode)
     !WRITE(lures,600,ERR=9999) ielem
     !WRITE(lures,900,ERR=9999)
     WRITE(lures,"(5x,i5,3e15.8)",ERR=9999) (i,x(:,i),i=1,8)
     !iSTOP = 1
     STOP 'null or negative jacobian'
     RETURN
   END IF
   ! compute inverse jacobian
   !jaci(1,1) =  (gcova(2,2)*gcova(3,3)-gcova(2,3)*gcova(3,2))/j0
   !jaci(2,1) = -(gcova(2,1)*gcova(3,3)-gcova(3,1)*gcova(2,3))/j0
   !jaci(3,1) =  (gcova(2,1)*gcova(3,2)-gcova(2,2)*gcova(3,1))/j0
   !jaci(1,2) = -(gcova(1,2)*gcova(3,3)-gcova(1,3)*gcova(3,2))/j0
   !jaci(2,2) =  (gcova(1,1)*gcova(3,3)-gcova(3,1)*gcova(1,3))/j0
   !jaci(3,2) = -(gcova(1,1)*gcova(3,2)-gcova(1,2)*gcova(3,1))/j0
   jaci(1,3) =  (gcova(1,2)*gcova(2,3)-gcova(1,3)*gcova(2,2))/j0
   jaci(2,3) = -(gcova(1,1)*gcova(2,3)-gcova(1,3)*gcova(2,1))/j0
   jaci(3,3) =  (gcova(1,1)*gcova(2,2)-gcova(2,1)*gcova(1,2))/j0

   DO i=1,8
     nfdas(i,g) = DOT_PRODUCT(der(i,:),jaci(:,3))
   END DO

 END DO
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lcas18
