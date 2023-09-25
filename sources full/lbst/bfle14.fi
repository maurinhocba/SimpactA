 SUBROUTINE bfle14(lb,a,b,sides,is,bb,t,t3,h,ci,nonrg)
 !
 !     Computes bending B matrix for element type 14
 !
 IMPLICIT NONE

 LOGICAL, INTENT(IN) :: sides(3), & !True = Side element exists
                        is(3),    & !True = Clamped side
                        nonrg       !True = non-regular meshes

 REAL (kind=8), INTENT(IN) :: lb,                &!A0/A == t/t0
                              a(3,0:3),b(3,0:3), &!Element side proyections
                              ci(3),             &!coeff for curvature comput.
                              t(3,2,0:3),        &!configuration derivatives
                              t3(3),             &!normal vector
                              h(3,3)              !h vectors

 REAL (kind=8), INTENT(OUT) ::  bb(3,6,3) !Bending B matrix

 INTEGER (kind=4) i,j,k,ns,ss(2)
 REAL (kind=8) d1(3),d2(3),n1,n2,c1,c2,mn(3,3)  ,r(3,2) !!
 REAL (kind=8) ::  bl(6,3) !local Bending B matrix

 !***  bending matrix evaluation

 bl = 0d0                      !zero B matrix
 ns = 0                        !initializes number of s.s. or free sides

 ! Contravariant base vectors of the main element (may be disregarded) !!
 d1(1) = (t(2,2,0)*t3(3) - t(3,2,0)*t3(2))*lb                          !!
 d1(2) = (t(3,2,0)*t3(1) - t(1,2,0)*t3(3))*lb                          !!
 d1(3) = (t(1,2,0)*t3(2) - t(2,2,0)*t3(1))*lb                          !!
 d2(1) = (t3(2)*t(3,1,0) - t3(3)*t(2,1,0))*lb                          !!
 d2(2) = (t3(3)*t(1,1,0) - t3(1)*t(3,1,0))*lb                          !!
 d2(3) = (t3(1)*t(2,1,0) - t3(2)*t(1,1,0))*lb                          !!
 ! ro components                                                       !!
 r(1,1) = DOT_PRODUCT(h(:,1) ,d1)                                      !!
 r(2,1) = DOT_PRODUCT(h(:,2) ,d1)                                      !!
 r(3,1) = DOT_PRODUCT(h(:,3) ,d1)                                      !!
 r(1,2) = DOT_PRODUCT(h(:,1) ,d2)                                      !!
 r(2,2) = DOT_PRODUCT(h(:,2) ,d2)                                      !!
 r(3,2) = DOT_PRODUCT(h(:,3) ,d2)                                      !!

 DO i=1,3                      ! for each side
   !  Contributions from normal displacements of the Main element (may be disregarded) !!
   bl(i,1) = bl(i,1) + r(1,1)*b(i,0) - r(1,2)*a(i,0)          !K11                     !!
   bl(i,2) = bl(i,2) + r(2,1)*b(i,0) - r(2,2)*a(i,0)          !K22                     !!
   bl(i,3) = bl(i,3) + r(3,1)*b(i,0) - r(3,2)*a(i,0)          !2K12                    !!
   IF (sides(i)) THEN                 ! side element exists
     !  Contributions from normal displacements of the Adj elements
     DO j=1,3                        ! for each node
       k = kk(j,i)                   ! local node (adjacent element)
       bl(k,1)= bl(k,1) +  b(i,0)*b(j,i)                    !K11
       bl(k,2)= bl(k,2) +  a(i,0)*a(j,i)                    !K22
       bl(k,3)= bl(k,3) - (a(i,0)*b(j,i) + b(i,0)*a(j,i))   !2K12
     END DO
     IF( nonrg )THEN
       bl(i,1)= bl(i,1) +  b(i,0)*b(i,0)*ci(i)                  !K11
       bl(i,2)= bl(i,2) +  a(i,0)*a(i,0)*ci(i)                  !K22
       bl(i,3)= bl(i,3) - (a(i,0)*b(i,0) + b(i,0)*a(i,0))*ci(i) !2K12
     END IF

   ELSE IF( is(i) )THEN        !clamped side

     !  Contributions from displacements of the main element

     n1 =  b(2,i)                   !side normal components n1
     n2 =  b(3,i)                   !                       n2
     c1 = n1*n1 - n2*n2             ! n1^2 - n2^2
     c2 = 2D0*n1*n2                 ! 2 n1 n2
     d1 = c1*b(:,0) - c2*a(:,0)  !corrected cartesian derivatives X_n
     d2 = c2*b(:,0) + c1*a(:,0)  !corrected cartesian derivatives X_s
     DO j=1,3                        ! for each node
       bl(j,1)= bl(j,1) - b(i,0)*d1(j)                 !K11
       bl(j,2)= bl(j,2) + a(i,0)*d2(j)                 !K22
       bl(j,3)= bl(j,3) + a(i,0)*d1(j) - b(i,0)*d2(j)  !2K12
     END DO

   ELSE ! simple supported or free
     !  Contributions from normal displacements of the Main elements
     DO j=1,3                        ! for each node
       bl(j,1)= bl(j,1) +  b(i,0)*b(j,0)                    !K11
       bl(j,2)= bl(j,2) +  a(i,0)*a(j,0)                    !K22
       bl(j,3)= bl(j,3) - (a(i,0)*b(j,0) + b(i,0)*a(j,0))   !2K12
     END DO
     !**** comment next two lines to avoid case below
     ns = ns+1                 !increase number of s.s. or free sides
     ss(ns) = i                !keep side order
   END IF
 END DO

 ! correct B matrix to nullify normal curvatures for free or simple-supported sides

 !ns = 0               !to avoid next case
 SELECT CASE(ns)      !according to the number of free or ss sides

 CASE (1)             !one simple supported side
   i = ss(1)                           !side order
   !    IF nu included
   !mn = RESHAPE ((/ 1d0-r11*r11-nu*r11*r22,    -r22*r11-nu*r22*r22,    -2*r12*r11-nu*r12*r22,    &
   !                    -r11*r22-nu*r11*r11, 1d0-r22*r22-nu*r11*r22,    -2*r12*r22-nu*r11*r12,    &
   !                    -r11*r12+nu*r11*r12,    -r22*r12+nu*r12*r22,    1d0-2*r12+nu*r11*r22  /), &
   !              (/3,3/))
   !    without nu
   !mn =  1d0-r11*r11     -r22*r11     -2*r12*r11
   !         -r11*r22  1d0-r22*r22     -2*r12*r22
   !         -r11*r12     -r22*r12     1d0-2*r12
   mn = RESHAPE ((/ a(1,i),   a(2,i),   a(3,i),    &
                    a(2,i),   b(1,i),   b(2,i),    &
                  a(3,i)/2, b(2,i)/2,   b(3,i) /), (/3,3/))
   CALL proma2(bl(1,1),bl(1,1),mn(1,1),6,3,3)

 CASE (2)         !Two sides free or simple supported
   i = ss(1)
   j = ss(2)
   !mn =  r11*t11, r11*t22, r11*t12,
   !      r22*t11, r22*t22, r22*t12,
   !      r12*t11, r12*t22, r12*t12
   mn = RESHAPE ((/ a(1,i),   a(2,i),   a(3,i),  &
                    b(1,i),   b(2,i),   b(3,i),  &
                    a(1,j),   a(2,j),   a(3,j) /), (/3,3/))
   CALL proma2(bl(1,1),bl(1,1),mn(1,1),6,3,3)

 END SELECT

 DO j=1,3        !for each curvature
   DO i=1,6      !for each node in the patch
     bb(:,i,j) = bl(i,j)*t3
   END DO
 END DO

 RETURN
 END SUBROUTINE bfle14
