 SUBROUTINE stra24(a,b,c,cd,ci,a0,a1,t,lb,x,sides,is,mems,first,stran,gamma,k0)

 !     Compute first and second fundamental forms for element BBST++

 IMPLICIT NONE

 REAL (kind=8), INTENT(IN) :: a(3),b(3),  & !side proyections on local x1
                              x(3,6),     & !nodal cooordinates
                              cd(4,2,3),  & !cartesian nodal derivatives at sides
                              ci(3),      & !stiffness factors
                              a1,         & !original area
                              a0(3)         !original angles between elements
 REAL (kind=8), INTENT(IN OUT) :: c(3,3,2)      !side proyections in normal dir
 REAL (kind=8), INTENT(IN), OPTIONAL ::  k0     !thermal curvature
 REAL (kind=8), INTENT(OUT), OPTIONAL :: stran(6), & !1st & 2nd fund. forms
                                         gamma(3) !distorsion angles
 REAL (kind=8), INTENT(OUT) :: t(3,15),& !derivatives of the element configuration
                               lb       !thickness ratio
 INTEGER(kind=4), INTENT(IN) :: sides(3) !1 = side elements exists
 LOGICAL, INTENT(IN) :: is(3),     & !True = clamped side
                        mems(3)      !True = use quadratic approach at that side
 LOGICAL, INTENT(IN), OPTIONAL :: first        !True = First time computation

 INTEGER (kind=4) :: i,j,l,n,n1,m,m1,ns,ss(2)
 REAL (kind=8) :: d1,d2,knn,r11,r12,r22,ml,mn(3,3),ls(3),nu, &
                  no(2,2),nd(2,2),t11,t12,t22,det,s(3),ni(3)


 t(1:3,1) = -MATMUL(x(1:3,1:3),b(1:3))     ! x(1)
 t(1:3,2) = +MATMUL(x(1:3,1:3),a(1:3))     ! x(2)
 t(1,3) = t(2,1)*t(3,2) - t(3,1)*t(2,2)    ! normal * rA
 t(2,3) = t(3,1)*t(1,2) - t(1,1)*t(3,2)
 t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
 !computes thickness ratio  t/t0
 lb = 1d0/SQRT(t(1,3)*t(1,3)+t(2,3)*t(2,3)+t(3,3)*t(3,3))
 t(1:3,3) = t(1:3,3)*lb                    ! normalizes normal vector

 n = 4           !first position for side element vectors in array T
 m = 10          !first position for grad vectors in array T
 ns = 0          !initializes number of s.s. or free sides
 DO i=1,3        !for each triangle in the patch
   s = x(:,hh(2,i))-x(:,hh(3,i))                !side vector  j =>k
   ls(i) = SQRT(s(1)*s(1)+s(2)*s(2)+s(3)*s(3))  !length
   s = s/ls(i)                                  !unit side vector
   ls(i) = 4d0*a1/lb/ls(i)                      !4A/Ls   = 4A(0)/Ls(0)/Lambda/Lambda_s
   t(1,n) = s(2)*t(3,3) - s(3)*t(2,3)    ! normal * rA
   t(2,n) = s(3)*t(1,3) - s(1)*t(3,3)
   t(3,n) = s(1)*t(2,3) - s(2)*t(1,3)
   knn = SQRT(t(1,n)*t(1,n)+t(2,n)*t(2,n)+t(3,n)*t(3,n))  !length
   t(:,n) = t(:,n)/knn                          !unit outward in-plane normal
   n1 = n + 1                                   !pointer to normal at side element
   m1 = m + 1                                   !position of second derivative
   IF(sides(i) == 1)THEN                !if side exists
     ni  = MATMUL(x(:,hh(:,i)),c(:,i,2))        !in-plane normal derivative at adjacent element
     t(1,n1) = s(2)*ni(3) - s(3)*ni(2)          ! element normal * rA
     t(2,n1) = s(3)*ni(1) - s(1)*ni(3)
     t(3,n1) = s(1)*ni(2) - s(2)*ni(1)
     knn = SQRT(t(1,n1)*t(1,n1)+t(2,n1)*t(2,n1)+t(3,n1)*t(3,n1))  !unit normal
     t(:,n1) = t(:,n1)/knn                      !unit normal at adjacent element
     ! tangent vectors at mid-side
     l = i+3                                    !extra node
     t(:,m:m1) = MATMUL(x(:,1:3),cd(1:3,:,i))   !x(1:2)(i) both derivatives
     t(:,m)    = t(:,m)  + x(:,l)*cd(4,1,i)     !x(1)(i)  add fourth node
     t(:,m1)   = t(:,m1) + x(:,l)*cd(4,2,i)     !x(2)(i)
   ELSE IF(sides(i) == 2)THEN           !side element is solid
     ni  = x(:,i+3)                     !normal vector
     c(1:3,i,2) = -ni                   !initial normal derivative
     t(1,n1) = s(2)*ni(3) - s(3)*ni(2)    ! normal * rA
     t(2,n1) = s(3)*ni(1) - s(1)*ni(3)
     t(3,n1) = s(1)*ni(2) - s(2)*ni(1)
     knn = SQRT(t(1,n1)*t(1,n1)+t(2,n1)*t(2,n1)+t(3,n1)*t(3,n1))  !unit normal
     t(:,n1) = t(:,n1)/knn                !unit normal at ficticious adjacent element
   !ELSE
   !  t(:,n1) = 0d0
   !  t(:,m:m1) = 0d0
   END IF
   n = n+2               !pointer to next side vectors in array T
   m = m+2               !pointer to next grad vectors in array T
 END DO

 IF( PRESENT(stran) ) THEN ! strains required

   ! Bending part (2nd fundamental form multiplied by 2)
   stran(1:6) = 0d0   !initializes
   n = 4              !pointer to vectors on side elements in T array
   m = 10             !pointer to vectors on side elements in T array
   DO i=1,3           !for each side element
     ml = ls(i)       !factor for sides
     IF(sides(i) == 1 )THEN      !if side element exist
       n1 = n + 1               !pointer no normal at side element
       m1 = m+1                               !position of second derivative
       d1 = DOT_PRODUCT(t(:,3),t(:,n1))       !cosine angle tM . ti
       d2 = DOT_PRODUCT(t(:,n),t(:,n1))       !sin    angle nM . ti
       gamma(i) = gamma(i)/ci(i)
       gamma(i) = ci(i)*(ATAN4(d2,d1,a0(i),gamma(i)))  !distorsion
       knn = gamma(i)*ml                      !contribution to curvature
       stran(4) = stran(4) +     b(i)*b(i)*knn    ! k11
       stran(5) = stran(5) +     a(i)*a(i)*knn    ! k22
       stran(6) = stran(6) - 2d0*a(i)*b(i)*knn    !2k12

       ! Membrane part (1st fundamental form, metric tensor at mid-surface)
       IF( mems(i) )THEN
         d1 = ABS(a0(i)+2d0*gamma(i))           !angle at boundary
         IF( d1 <= alp1 )THEN            !smooth  use side (cuadratic) metric tensor
           stran(1) = stran(1) + DOT_PRODUCT(t(:,m),t(:,m))   !g11
           stran(2) = stran(2) + DOT_PRODUCT(t(:,m1),t(:,m1)) !g22
           stran(3) = stran(3) + DOT_PRODUCT(t(:,m),t(:,m1))  !g12

         ELSE IF( d1 >= alp2 )THEN       !folded, use main triangle metric tensor
           stran(1) = stran(1) + DOT_PRODUCT(t(:,1),t(:,1))   !g11
           stran(2) = stran(2) + DOT_PRODUCT(t(:,2),t(:,2))   !g22
           stran(3) = stran(3) + DOT_PRODUCT(t(:,1),t(:,2))   !g12

         ELSE                            !interpolate
           d1 = (d1-alp1)/alp3           !function for main triangle
           d2 = 1d0 - d1                 !function for cuadratic interpolation
           stran(1) = stran(1) + d1*DOT_PRODUCT(t(:,1),t(:,1)) &
                               + d2*DOT_PRODUCT(t(:,m),t(:,m))   !g11
           stran(2) = stran(2) + d1*DOT_PRODUCT(t(:,2),t(:,2)) &
                               + d2*DOT_PRODUCT(t(:,m1),t(:,m1)) !g22
           stran(3) = stran(3) + d1*DOT_PRODUCT(t(:,1),t(:,2)) &
                               + d2*DOT_PRODUCT(t(:,m),t(:,m1))  !g12
         END IF

       ELSE
         stran(1) = stran(1) + DOT_PRODUCT(t(:,1),t(:,1))   !g11
         stran(2) = stran(2) + DOT_PRODUCT(t(:,2),t(:,2))   !g22
         stran(3) = stran(3) + DOT_PRODUCT(t(:,1),t(:,2))   !g12

       END IF

     ELSE IF( sides(i) == 0 .OR. sides(i) == 2 )THEN           !boundary node
       IF( is(i) )THEN  !If clamped side
         !normal vector to plane of symmetry (outward)
         IF(first .AND. sides(i) == 0) c(1:3,i,2) = t(1:3,n)  ! initial normal derivative
         !proyection of normal to the triangle over the normal to the plane
         d1 = DOT_PRODUCT(t(:,n),c(1:3,i,2))   !cos
         d2 =-DOT_PRODUCT(t(:,3),c(1:3,i,2))   !sin
         gamma(i) = gamma(i)/ci(i)
         gamma(i) = ci(i)*(ATAN4(d2,d1,a0(i),gamma(i))) !distorsion
         knn = gamma(i)*ml                     !contribution to curvature
         stran(4) = stran(4) +     b(i)*b(i)*knn    ! k11
         stran(5) = stran(5) +     a(i)*a(i)*knn    ! k22
         stran(6) = stran(6) - 2d0*a(i)*b(i)*knn    !2k12
       ELSE                  !simple supported or free side

         ns = ns+1           !increase counter of SS or Free sides
         ss(ns) = i          !keep side number
       END IF
       !non existing side, use main triangle metric tensor
       stran(1) = stran(1) + DOT_PRODUCT(t(:,1),t(:,1))   !g11
       stran(2) = stran(2) + DOT_PRODUCT(t(:,2),t(:,2))   !g22
       stran(3) = stran(3) + DOT_PRODUCT(t(:,1),t(:,2))   !g12

     ELSE  ! sides(i) <  0 !branching node
       knn = gamma(i)*ml                      !contribution to curvature
       stran(4) = stran(4) +     b(i)*b(i)*knn    ! k11
       stran(5) = stran(5) +     a(i)*a(i)*knn    ! k22
       stran(6) = stran(6) - 2d0*a(i)*b(i)*knn    !2k12
       ! use main triangle metric tensor
       stran(1) = stran(1) + DOT_PRODUCT(t(:,1),t(:,1))   !g11
       stran(2) = stran(2) + DOT_PRODUCT(t(:,2),t(:,2))   !g22
       stran(3) = stran(3) + DOT_PRODUCT(t(:,1),t(:,2))   !g12
     END IF
     n = n+2               !increase pointer in T array
     m = m+2               !pointer to next grad vectors in array T
   END DO

   IF( ANY(sides < 0)) ns = 0  !if a branching side exists
   SELECT CASE(ns)      !if SS or Free sides exist

   CASE (1)          !one side free or simple supported
     i = ss(1)                              !side number
     ml  = SQRT(b(i)**2+a(i)**2)        !side pseudo length
     d1 =  b(i)/ml                        !side normal
     d2 = -a(i)/ml
     r11 = d1*d1                            !auxiliar factors
     r12 = d1*d2
     r22 = d2*d2
!     mn = RESHAPE ((/ 1d0-r11*r11,    -r22*r11,    -2*r12*r11,    &
!                         -r11*r22, 1d0-r22*r22,    -2*r12*r22,    &
!                         -r11*r12,    -r22*r12,    1d0-2*r12  /), &
!                   (/3,3/))
     nu = ci(i)      !Poisson ratio
   mn = RESHAPE ((/ 1d0-r11*r11-nu*r11*r22,    -r22*r11-nu*r22*r22,    -2*r12*r11-nu*r12*r22,    &
&                      -r11*r22-nu*r11*r11, 1d0-r22*r22-nu*r11*r22,    -2*r12*r22-nu*r11*r12,    &
&                      -r11*r12+nu*r11*r12,    -r22*r12+nu*r12*r22,    1d0-2*r12+nu*r11*r22  /), &
&                (/3,3/))
     stran(4:6) = MATMUL(mn,stran(4:6))  !modification to nullify normal curvature
     IF(PRESENT(k0)) stran(4:6) = stran(4:6) + k0* (/ r11, r22, 2d0*r12 /) !modification to include thermal curvature

   CASE (2)         !Two sides free or simple supported
     DO j=1,2       ! compute the normal vectors
       i = ss(j)    ! free side
       ml  = SQRT(b(i)**2+a(i)**2)         !side pseudo length
       no(1:2,j) = (/ b(i), -a(i) /)/ml    !side normal
     END DO
     ! second compute the dual base
     det = no(1,1)*no(2,2) - no(1,2)*no(2,1)   !determinant
     nd(1,1) =  no(2,2)/det                    !inverse
     nd(2,1) = -no(1,2)/det
     nd(1,2) = -no(2,1)/det
     nd(2,2) =  no(1,1)/det

     r11 = nd(1,1)*nd(1,2)                     !auxiliar factors
     r22 = nd(2,1)*nd(2,2)
     r12 = nd(1,1)*nd(2,2)+nd(2,1)*nd(1,2)
     t11 = 2d0*no(1,1)*no(1,2)
     t22 = 2d0*no(2,1)*no(2,2)
     t12 = no(1,1)*no(2,2)+no(2,1)*no(1,2)
     mn = RESHAPE ((/ t11*r11, t11*r22, t11*r12,  &
                      t22*r11, t22*r22, t22*r12,  &
                      t12*r11, t12*r22, t12*r12 /), (/3,3/))
     stran(4:6) = MATMUL(mn,stran(4:6)) !modification to nullify normal curvatures
     IF(PRESENT(k0)) THEN
       knn = k0/(1d0+(r11+r22)**2)
       r11 = no(1,1)*no(1,1)+no(1,2)*no(1,2)
       r22 = no(2,1)*no(2,1)+no(2,2)*no(2,2)
       r12 =(no(1,1)*no(2,1)+no(1,2)*no(2,2))*2d0
       stran(4:6) = stran(4:6) + knn*(/ r11,r22,r12  /)
     END IF
   END SELECT

   stran(1:3) = stran(1:3)/3d0               !average metric tensor
   stran(4) = 2d0*stran(4)                   !2 k11
   stran(5) = 2d0*stran(5)                   !2 k22
 END IF
 t(:,8) = ls  !keep factors
 RETURN
 END SUBROUTINE stra24
