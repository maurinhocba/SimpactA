 SUBROUTINE stra14(a,b,x,t3,t,h,sides,is,ci,nonrg,stran,lb)

 !     Compute first and second fundamental forms for element CST-BST (TLF)
 IMPLICIT NONE

 REAL (kind=8), INTENT(IN) :: x(3,6),    &  !nodal cooordinates
                              b(3,0:3),  &  !cartesian derivatives local x1
                              a(3,0:3),  &  !cartesian derivatives local x2
                              ci(3)         !coeff for curvature computation
 REAL (kind=8), INTENT(OUT) :: t(3,2,0:3), & !derivatives of the element configuration
                               t3(3),      & !element normal
                               h(3,3)        !h vectors
 LOGICAL, INTENT(IN) :: sides(3),  & !True = side elements exists
                        is(3),     & !True = clamped side
                        nonrg        !True = non-regular meshes
 REAL (kind=8), INTENT(OUT), OPTIONAL :: stran(6), & !1st & 2nd fund. forms
                                         lb          !thickness ratio


 INTEGER (kind=4) :: i,j,ns,ss(2)
 REAL (kind=8) :: ln,mn(3,3),ts(3),tn(3),lbb



 !for main triangle
 t(:,1,0) = -MATMUL(x(:,1:3),b(:,0))     ! x(1)
 t(:,2,0) = +MATMUL(x(:,1:3),a(:,0))     ! x(2)
 CALL vecpro(t(1,1,0),t(1,2,0),t3(1))    ! normal * rA
 CALL vecuni(3,t3(1),lbb)                 ! normalizes normal vector

 ns = 0          !initializes number of s.s. or free sides
 h  = 0d0        !initializes h vectors
 DO i=1,3        !for each side
   IF(sides(i))THEN    !if side exists
     t(:,1,i) = -MATMUL(x(:,kk(:,i)),b(:,i))     ! x(1)
     t(:,2,i) = +MATMUL(x(:,kk(:,i)),a(:,i))     ! x(2)
     IF( nonrg )THEN  !add the influence of the opposite node
       t(:,1,i) = t(:,1,i) - ci(i)*b(i,0)*x(:,i) ! x(1)
       t(:,2,i) = t(:,2,i) + ci(i)*a(i,0)*x(:,i) ! x(2)
     END IF
   ELSE IF (is(i)) THEN  !clamped side
     ts = (x(:,hh(1,i)) - x(:,hh(2,i)))*b(1,i)   !side derivative
     tn = b(2,i)*t(:,1,0) + b(3,i)*t(:,2,0)      !normal derivative (main element)
     ln = SQRT(DOT_PRODUCT(tn,tn))               !normal stretch
     tn = 2d0*ln*a(:,i) - tn                     !normal derivative (side)
     t(:,1,i) =  b(2,i)*tn - b(3,i)*ts           !final gradient at the side
     t(:,2,i) =  b(3,i)*tn + b(2,i)*ts
   ELSE                  !simple supported or free
     ns = ns+1             !increase counter of SS or Free sides
     ss(ns) = i            !keep side number
     t(:,:,i) = t(:,:,0)   !use main triangle
   END IF
   ! vectors to compute Bending part (2nd fundamental form)
   h(:,1) = h(:,1) - b(i,0)*t(:,1,i)                    ! h11
   h(:,2) = h(:,2) + a(i,0)*t(:,2,i)                    ! h22
   h(:,3) = h(:,3) + a(i,0)*t(:,1,i) - b(i,0)*t(:,2,i)  !2h12
 END DO

 IF( .NOT.PRESENT(stran) )RETURN   ! Only vectors T and H required

 lb = 1d0/lbb                             !computes thickness ratio  t/t0
 ! Membrane part (1st fundamental form, metric tensor at mid-surface)
 stran(1) = DOT_PRODUCT(t(:,1,0),t(:,1,0))
 stran(2) = DOT_PRODUCT(t(:,2,0),t(:,2,0))
 stran(3) = DOT_PRODUCT(t(:,1,0),t(:,2,0))

 ! Bending part (2nd fundamental form multiplied by 2)
 stran(4) = DOT_PRODUCT(h(:,1),t3)
 stran(5) = DOT_PRODUCT(h(:,2),t3)
 stran(6) = DOT_PRODUCT(h(:,3),t3)

 !ns = 0  !to avoid next case
 SELECT CASE(ns)      !if SS or Free sides exist

 CASE (1)          !one side free or simple supported
   i = ss(1)
   mn = RESHAPE ((/ a(1,i),   a(2,i),   a(3,i),    &
                    a(2,i),   b(1,i),   b(2,i),    &
                  a(3,i)/2, b(2,i)/2,   b(3,i) /), (/3,3/))
   stran(4:6) = MATMUL(mn,stran(4:6)) !modification to nullify normal curvature

 CASE (2)         !Two sides free or simple supported
   i = ss(1)
   j = ss(2)
   mn = RESHAPE ((/ a(1,i),   a(2,i),   a(3,i),  &
                    b(1,i),   b(2,i),   b(3,i),  &
                    a(1,j),   a(2,j),   a(3,j) /), (/3,3/))
   stran(4:6) = MATMUL(mn,stran(4:6)) !modification to nullify normal curvatures

 END SELECT

 stran(4) = 2d0*stran(4)                             !2 k11
 stran(5) = 2d0*stran(5)                             !2 k22

 RETURN
 END SUBROUTINE stra14
