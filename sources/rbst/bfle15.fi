 SUBROUTINE bfle15(a,b,c,ci,sides,is,bb,t)
 !
 !     Computes bending B matrix for element type 15
 !
 IMPLICIT NONE

 INTEGER(kind=4), INTENT(IN) :: sides(3)    !1 = Side element exists
 LOGICAL, INTENT(IN) ::  is(3)       !True = Clamped side

 REAL (kind=8), INTENT(IN) :: a(3),b(3),c(3,3,2), & !Element side proyections
                              ci(3),              & !factors
                              t(3,9)                !configuration derivatives

 REAL (kind=8), INTENT(OUT) ::  bb(3,6,3) !Bending B matrix

 INTEGER (kind=4) i,j,kk,n,ns,ss(2)
 REAL (kind=8) c1,c2,c11,c12,c22,r11,r12,r22,ml,nu
 REAL (kind=8) xb(3),mn(3,3),no(2,2),nd(2,2)

 INTEGER (kind=4), SAVE ::  &  !side element connectivities
           k(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )

 !***  bending matrix evaluation

 bb = 0d0                      !zero B matrix
 ns = 0                        !initializes number of s.s. or free sides
 n = 4                         ! position of side derivative 1
 DO i=1,3                      ! for each side
   ml = t(i,4)*ci(i)           !2A/Ls  = 2A(0)/Ls(0)/Lambda/Lambda_s *ci
   r11 =      b(i)*b(i)*ml
   r22 =      a(i)*a(i)*ml
   r12 = -2d0*b(i)*a(i)*ml
   IF (sides(i) == 1) THEN                 ! side element exists
     !  Contributions when side element exist
     DO j=1,3                        ! for each node
       !   contributions from main triangle to distorsion
       xb = c(j,i,1)*t(1:3,3)                !N(j)'n * t3(M)
       bb(1:3,j,1)= bb(1:3,j,1)   + r11*xb   !K11
       bb(1:3,j,2)= bb(1:3,j,2)   + r22*xb   !K22
       bb(1:3,j,3)= bb(1:3,j,3)   + r12*xb   !K12
       !  contributions from side triangle to distorsion
       kk = k(j,i)                   ! local node (adjacent element)
       xb = c(j,i,2)*t(1:3,n+1)              !N(j)'n * t3(i)
       bb(1:3,kk,1)= bb(1:3,kk,1) + r11*xb   !K11
       bb(1:3,kk,2)= bb(1:3,kk,2) + r22*xb   !K22
       bb(1:3,kk,3)= bb(1:3,kk,3) + r12*xb   !K12
     END DO

   ELSE IF(sides(i) == 0) THEN       ! boundary side
     IF( is(i) )THEN        !clamped side

       DO j=1,3                        ! for each node
         !  contributions from main triangle only
         xb = c(j,i,1)*t(1:3,3)                 !   to distorsion
         bb(1:3,j,1)= bb(1:3,j,1)   + r11*xb   !K11
         bb(1:3,j,2)= bb(1:3,j,2)   + r22*xb   !K22
         bb(1:3,j,3)= bb(1:3,j,3)   + r12*xb   !K12
       END DO

     ELSE ! simple supported or free
       ns = ns+1                 !increase number of s.s. or free sides
       ss(ns) = i                !keep side order
     END IF
   END IF
   n = n+2                     ! position of side derivative 1
 END DO
 ! this is done only when SIDES(1:3) >= 0
 SELECT CASE(ns)      !according to the number of free or ss sides

 CASE (1)             !one simple supported side
   i = ss(1)                           !side order
   ml  = SQRT(b(i)**2+a(i)**2)     !side pseudo length
   c1 =  b(i)/ml                     !normalizes side to compute normal
   c2 = -a(i)/ml
   r11 = c1*c1                         !auxiliar factors
   r12 = c1*c2
   r22 = c2*c2
!   mn = RESHAPE ((/ 1d0-r11*r11,    -r22*r11,    -2*r12*r11,    &
!                       -r11*r22, 1d0-r22*r22,    -2*r12*r22,    &
!                       -r11*r12,    -r22*r12,     1d0-2*r12 /), &
!                 (/3,3/))
   nu = ci(i)           !Poisson ratio
   mn = RESHAPE ((/ 1d0-r11*r11-nu*r11*r22,    -r22*r11-nu*r22*r22,    -2*r12*r11-nu*r12*r22,    &
&                      -r11*r22-nu*r11*r11, 1d0-r22*r22-nu*r11*r22,    -2*r12*r22-nu*r11*r12,    &
&                      -r11*r12+nu*r11*r12,    -r22*r12+nu*r12*r22,    1d0-2*r12+nu*r11*r22  /), &
                 (/3,3/))
   CALL proma2(bb(1,1,1),bb(1,1,1),mn(1,1),18,3,3)

 CASE (2)         !Two sides free or simple supported
   DO j=1,2       ! compute the normal vectors
     i = ss(j)    ! free side
     ml  = SQRT(b(i)**2+a(i)**2)         !side pseudo length
     no(1:2,j) = (/ b(i), -a(i) /)/ml    !side normal
   END DO
   ! second compute the dual base
   ml = no(1,1)*no(2,2) - no(1,2)*no(2,1)    !determinant of the base
   nd(1,1) =  no(2,2)/ml                     !inverse of matrix NO
   nd(2,1) = -no(1,2)/ml
   nd(1,2) = -no(2,1)/ml
   nd(2,2) =  no(1,1)/ml
                                             !auxiliar factors
   r11 = nd(1,1)*nd(1,2)
   r22 = nd(2,1)*nd(2,2)
   r12 = nd(1,1)*nd(2,2)+nd(2,1)*nd(1,2)
   c11 = 2d0*no(1,1)*no(1,2)
   c22 = 2d0*no(2,1)*no(2,2)
   c12 = no(1,1)*no(2,2)+no(2,1)*no(1,2)
   mn = RESHAPE ((/ c11*r11, c11*r22, c11*r12,  &
                    c22*r11, c22*r22, c22*r12,  &
                    c12*r11, c12*r22, c12*r12 /), (/3,3/))
   CALL proma2(bb(1,1,1),bb(1,1,1),mn(1,1),18,3,3)

 END SELECT

 RETURN
 END SUBROUTINE bfle15
