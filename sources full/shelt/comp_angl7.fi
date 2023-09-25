 SUBROUTINE comp_angl7(nnode,    x,    t,   lb,nangl,ang,locax,jas)
 !***********************************************************************
 !
 !****this routine compute the angle (cos,sin) between local cartesyan system
 !    and nodal system
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 INTEGER (kind=4), INTENT(IN) :: nnode,locax

 REAL (kind=8), INTENT(IN) :: x(:,:),lb(:,:,:),t(:,:),ang

 REAL (kind=8), INTENT(OUT) :: nangl(:,:)
 REAL (kind=8), POINTER, INTENT(OUT) :: jas(:,:)

 !                   local variables

 INTEGER (kind=4) :: n,j,g
 REAL (kind=8) :: t1(3),t2(3),t3(3),j0,lt,cb,gcova(3,2)
 ! nodal shape functions derivatives 6: number of nodes (shape functions)
 !                                   2: number of in-plane coordinates
 !                                   6: number of sampling points
 !                                   2: 1-nodes
 !                                      2-side points
 REAL (kind=8), PARAMETER :: der4(6,2,6,2) = (/                                       &
                -3D0, -1D0, 0D0,  4D0,  0D0,  0D0, -3D0,  0D0, -1D0,  0D0,  0D0,  4D0,&
                 1D0,  3D0, 0D0, -4D0,  0D0,  0D0,  1D0,  0D0, -1D0, -4D0,  4D0,  0D0,&
                 1D0, -1D0, 0D0,  0D0,  4D0, -4D0,  1D0,  0D0,  3D0,  0D0,  0D0, -4D0,&
                -1D0,  1D0, 0D0,  0D0,  0D0,  0D0, -1D0,  0D0, -1D0, -2D0,  2D0,  2D0,&
                 1D0,  1D0, 0D0, -2D0,  2D0, -2D0,  1D0,  0D0,  1D0, -2D0,  2D0, -2D0,&
                -1D0, -1D0, 0D0,  2D0,  2D0, -2D0, -1D0,  0D0,  1D0,  0D0,  0D0,  0D0,&
 -2.1547005384,  -0.1547005384,  0.0000000000,   2.3094010768,   0.0000000000,   0.0000000000, &
 -2.1547005384,   0.0000000000, -1.0000000000,  -0.8452994616,   0.8452994616,   3.1547005384, &
  0.1547005384,   2.1547005384,  0.0000000000,  -2.3094010768,   0.0000000000,   0.0000000000, &
  0.1547005384,   0.0000000000, -1.0000000000,  -3.1547005384,   3.1547005384,   0.8452994616, &
  1.0000000000,   2.1547005384,  0.0000000000,  -3.1547005384,   0.8452994616,  -0.8452994616, &
  1.0000000000,   0.0000000000, -0.1547005384,  -3.1547005384,   3.1547005384,  -0.8452994616, &
  1.0000000000,  -0.1547005384,  0.0000000000,  -0.8452994616,   3.1547005384,  -3.1547005384, &
  1.0000000000,   0.0000000000,  2.1547005384,  -0.8452994616,   0.8452994616,  -3.1547005384, &
  0.1547005384,  -1.0000000000,  0.0000000000,   0.8452994616,   3.1547005384,  -3.1547005384, &
  0.1547005384,   0.0000000000,  2.1547005384,   0.0000000000,   0.0000000000,  -2.3094010768, &
 -2.1547005384,  -1.0000000000,  0.0000000000,   3.1547005384,   0.8452994616,  -0.8452994616, &
 -2.1547005384,   0.0000000000, -0.1547005384,   0.0000000000,   0.0000000000,   2.3094010768 /)

 INTERFACE
    INCLUDE 'vecpro.h'
    INCLUDE 'vecuni.h'
 END INTERFACE

 g = 0
 DO n=1,nnode    !for each sampling point
   !       compute the local cartesian system
   DO j=1,2         !1: node , 2:mid-side
     gcova = MATMUL(x,der4(:,:,n,j))         !  convective base
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
     IF( j == 1 )THEN     !local systems at nodes
       t1 =  COS(ang)*t1 + SIN(ang)*t2   !         !t_x
       cb = DOT_PRODUCT(t3,t(:,n))                 !COS beta
       nangl(1,n) = DOT_PRODUCT(lb(:,1,n),t1)/cb   !COS alpha
       nangl(2,n) = DOT_PRODUCT(lb(:,2,n),t1)/cb   !SIN alpha
     ELSE                 !local systems at mid-side points
       t3 = t1                          !auxiliar
       t1 =  COS(ang)*t1 + SIN(ang)*t2  !t_x
       t2 =  COS(ang)*t2 - SIN(ang)*t3  !t_y
       g = g+1
       !local derivative along sides at shear mid-points
       SELECT CASE (n)
       CASE (1:2)               !xita (A-B points)
         jas(1,g) =  DOT_PRODUCT(gcova(:,1),t1)
         jas(2,g) =  DOT_PRODUCT(gcova(:,1),t2)
       CASE (3:4)               !(-xita+eta) (C-D points)
         jas(1,g) =  -DOT_PRODUCT(gcova(:,1),t1) +DOT_PRODUCT(gcova(:,2),t1)
         jas(2,g) =  -DOT_PRODUCT(gcova(:,1),t2) +DOT_PRODUCT(gcova(:,2),t2)
       CASE (5:6)               !-eta  (A-C points)
         jas(1,g) =  -DOT_PRODUCT(gcova(:,2),t1)
         jas(2,g) =  -DOT_PRODUCT(gcova(:,2),t2)
       END SELECT
     END IF
   END DO
 END DO

 RETURN
 END SUBROUTINE comp_angl7
