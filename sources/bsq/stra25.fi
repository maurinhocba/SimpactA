 SUBROUTINE stra25(cd,cn,ns,t,t0,tb,tns,a0,ci,x,sides,is,h,hb,dt, &
                   area,gamma,lb,stran,strag,fff)

 !     Compute first and second fundamental forms for element NBSQ (TLF)

 IMPLICIT NONE

 REAL (kind=8), INTENT(IN) :: ns(2,4),  &   !side proyections on local x1
                              x(3,12), &   !nodal cooordinates of the patch
                              a0(4),   &   !original angles between elements
                              ci(4)        !stiffness factors
 REAL (kind=8), INTENT(IN OUT) :: cd(4,2,0:8), &!cartesyan derivatives of shape functions
                                  cn(4,2,1:8)      !normal derivatives N_n
 REAL (kind=8), INTENT(OUT) :: t(3,2,0:8), & !derivatives of the element configuration
                               tns(3,2,8), & !normal derivatives of the element configuration
                               h(3,3),     & !h vectors
                               hb(3),    & !h side  vectors
                               t0(3),      & !normal at the center (main element)
                               tb(3,8)      !normal at the side (adyacent element)
 INTEGER(kind=4), INTENT(IN) :: sides(4) !1 = side elements exists
 LOGICAL, INTENT(IN) :: is(4)        !True = clamped side
 REAL (kind=8), INTENT(OUT), OPTIONAL :: fff(3,4),stran(9),strag(4), & !1st & 2nd fund. forms
                                         lb,dt(4)                   !thickness ratio
 REAL (kind=8), INTENT(IN), OPTIONAL ::  area(0:4)            !areas
 REAL (kind=8), INTENT(OUT), OPTIONAL :: gamma(4) !distorsion angles
 !LOGICAL, OPTIONAL, INTENT(IN) :: first        !True = First time computation

 INTEGER (kind=4) :: i,ii
 REAL (kind=8) :: lbl,lbi,lbn,knn(4),ls,fs(4),ml,t3s(3),  &
                  n1,n2,nl,c1,c2,c3,c4,c5,cosa,sina
 REAL (kind=8) :: dk(3,4),d0(3),strak(4)
 ! ******** COMPUTE GRADIENTS AT EACH POINT **************
 tns = 0d0
 knn = 0d0
 h = 0d0
 hb = 0d0
 !  element center
 t(:,:,0) = MATMUL(x(:,1:4),cd(:,:,0))                !x_1_0, x_2_0
 CALL vecpro(t(:,1,0),t(:,2,0),t0)                    !normal direction
 CALL vecuni(3,t0,lbl)                                !computes thickness ratio  t0/t

 d0 = (x(:,1)+x(:,2)+x(:,3)+x(:,4))/4d0
 DO i=1,4        !for each side
   t(:,:,i) = MATMUL(x(:,1:4),cd(:,:,i))              ! x_1_i, x_2_i
   ls = SQRT(ns(1,i)**2+ns(2,i)**2)                   !side length
   fs(i) =  lbl!/2d0/area(0)                         !Ls/2A   = Ls(0)*Lambda*Lambda_s/2A(0)
   tns(:,:,i) = MATMUL(x(:,1:4),cn(:,:,i))                    !t_n
   ! normal at the side
   CALL vecpro(tns(:,1,i),tns(:,2,i),t3s)         !normal direction
   CALL vecuni(3,t3s,lbi)                           !computes thickness ratio  ti0/t
   tb(:,i) = t3s                                   !keep normal element
   dk(:,i) = d0 - (x(:,nextn(i))+x(:,nextn(i+3)))/2d0
   dt(i) = SQRT(dk(1,i)*dk(1,i)+dk(2,i)*dk(2,i)+dk(3,i)*dk(3,i))
   !
   nl = SQRT(ns(1,i)**2+ns(2,i)**2)           !side length
   n1 = ns(1,i)/nl           !normalize normal vector
   n2 = ns(2,i)/nl
   !    side rotation
   c1 = n1*n2*nl                ! n1 n2
   c2 = (n1**2-n2**2)*nl          ! n1^2-n2^2
   !    normal to the side rotation
   c3 = 2d0*c1              !2 n1 n2
   c4 = n1**2*nl             ! n1^2
   c5 = n2**2*nl             ! n2^2
   !I-Part:
   h(:,1) = h(:,1) + c4*tns(:,1,i)  - c1*tns(:,2,i)      ! h11
   h(:,2) = h(:,2) + c5*tns(:,1,i)  + c1*tns(:,2,i)      ! h22
   h(:,3) = h(:,3) + c3*tns(:,1,i)  + c2*tns(:,2,i)      !2h12

   ii = i+4      !side position in array "t"
   ml = fs(i)    !factor for sides
   IF(sides(i) == 1)THEN                              !side element exists
     tns(:,:,ii) = MATMUL(x(:,ln(:,i)),cn(:,:,ii))    !t_n
     ! normal at the side
     CALL vecpro(tns(:,2,ii),tns(:,1,ii),t3s)         !normal direction
     CALL vecuni(3,t3s,lbi)                           !computes thickness ratio  ti0/t
     tb(:,ii) = t3s                                   !keep normal element
     !
     cosa = +DOT_PRODUCT(tns(:,1,i),tns(:,1,ii))               !angle between normals
     sina = -DOT_PRODUCT(tns(:,1,i),t3s)
     !gamma(i) = ci(i)*ATAN2(cosa,sina)                      !angle
     gamma(i) = gamma(i)/ci(i)
     gamma(i) = ci(i)*atan4(sina,cosa,a0(i),gamma(i))         !angle change
     knn(i) = 2d0*gamma(i)*ml                       !contribution to curvature
     strak(i) = knn(i)/dt(i)/ml !/2d0

     !I-Part:
     h(:,1) = h(:,1) + c4*tns(:,1,i) - c1*tns(:,2,i)      ! h11
     h(:,2) = h(:,2) + c5*tns(:,1,i) + c1*tns(:,2,i)      ! h22
     h(:,3) = h(:,3) + c3*tns(:,1,i) + c2*tns(:,2,i)      !2h12
     !II-Part:
     hb(1) = hb(1) + c4*knn(i)
     hb(2) = hb(2) + c5*knn(i)
     hb(3) = hb(3) + c3*knn(i)
   ELSE IF( sides(i) == 0 )THEN    !boundary node
     IF( is(i) )THEN   !clamped side
       lbn = SQRT(DOT_PRODUCT(tns(:,1,i),tns(:,1,i)))                 !stretching in normal direction
       tns(:,1,ii) = 2d0*cn(1:3,2,ii)*lbn - tns(:,1,i)    !use original direction
       tns(:,2,ii) = tns(:,2,i)

     ELSE             !simple supported or free
       tns(:,:,ii) = tns(:,:,i)                       !same values
     END IF
     !I-Part:
     h(:,1) = h(:,1) + c4*tns(:,1,ii) - c1*tns(:,2,ii)      ! h11
     h(:,2) = h(:,2) + c5*tns(:,1,ii) + c1*tns(:,2,ii)      ! h22
     h(:,3) = h(:,3) + c3*tns(:,1,ii) + c2*tns(:,2,ii)      !2h12
   ELSE !IF(sides(i) < 0 ) THEN  branching side
     knn(i) = 2d0*gamma(i)*ml                             !contribution to curvature
     strak(i) = knn(i)/dt(i)/ml !/2d0
     !I-Part:
     h(:,1) = h(:,1) + c4*tns(:,1,i) - c1*tns(:,2,i)      ! h11
     h(:,2) = h(:,2) + c5*tns(:,1,i) + c1*tns(:,2,i)      ! h22
     h(:,3) = h(:,3) + c3*tns(:,1,i) + c2*tns(:,2,i)      !2h12
     !II-Part:
     hb(1) = hb(1) + c4*knn(i)
     hb(2) = hb(2) + c5*knn(i)
     hb(3) = hb(3) + c3*knn(i)
   END IF
 END DO

 ! ******** COMPUTE STRAIN MEASURES  **************

 IF( PRESENT(lb) )lb = 1d0/lbl        !thickness ratio
 IF( .NOT.PRESENT(fff) )RETURN
 DO i=1,4           !for each side element
   !   membrane strain measure, metric tensor C at each side
   fff(1,i) = DOT_PRODUCT(t(:,1,i),t(:,1,i))   !g11
   fff(2,i) = DOT_PRODUCT(t(:,2,i),t(:,2,i))   !g22
   fff(3,i) = DOT_PRODUCT(t(:,1,i),t(:,2,i))   !g12
 END DO
 ! 1st fundamental form, average of values computed at the sides
 stran(1:3) = MATMUL(fff,area(1:4)) !average C
 !fff(:,1) = fff(:,1) - stran(1:3)   !  hourglass C
 !fff(:,2) = fff(:,2) - stran(1:3)   !
 !fff(:,3) = fff(:,3) - stran(1:3)   !
 !fff(:,4) = fff(:,4) - stran(1:3)   !
 !fff(1:2,1:4) = fff(1:2,1:4)/2d0    !  Green-Lagrange
 ! 2nd fundamental form, t3_0 proyections of vectors hij
 stran(4) = -(DOT_PRODUCT(h(:,1),t0))/area(0)     !2k11
 stran(5) = -(DOT_PRODUCT(h(:,2),t0))/area(0)     !2k22
 stran(6) = -(DOT_PRODUCT(h(:,3),t0))/area(0)/2d0 !2k12
 stran(7) = -hb(1)/area(0)                      !2k11
 stran(8) = -hb(2)/area(0)                      !2k22
 stran(9) = -hb(3)/area(0)/2d0                  !2k12

 DO i=1,4
   nl = SQRT(ns(1,i)**2+ns(2,i)**2)       !side length
   n1 = ns(1,i)/nl           !normalize normal vector
   n2 = ns(2,i)/nl
   !    normal to the side rotation
   c1 = n1**2 !/2d0         ! n1^2
   c2 = n2**2 !/2d0         ! n2^2
   c3 = n1*n2         !2 n1 n2
   strag(i) = -strak(i) - (c1*stran(7)+c2*stran(8)+c3*stran(9))
 END DO

 RETURN
 END SUBROUTINE stra25
