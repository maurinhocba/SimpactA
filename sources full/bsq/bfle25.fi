 SUBROUTINE bfle25(cn,ci,ns,sides,is,solsi,bb,bbk,bbg,dt,t0,tb,area)
 !
 !     Computes bending B matrix for element BSQ
 !
 IMPLICIT NONE

 INTEGER(kind=4), INTENT(IN) :: sides(4)    !1 = Side element exists
 LOGICAL, INTENT(IN) :: is(4), &   !True = Clamped side
                        solsi(4)

 REAL (kind=8), INTENT(IN) :: cn(4,2,1:8),  & !normal derivatives
                              ci(4),        & !factors
                              ns(2,4),      & !Element side proyections
                              dt(4),        &
                              t0(3),        & !normal at the center (main element)
                              tb(3,8),      & !normal at the side
                              area(0:4)       !element area

 REAL (kind=8), INTENT(OUT) ::  bb(3,12,4), & !Bending B matrix
                                bbk(3,12,4),&
                                bbg(3,12,4)
 INTEGER (kind=4) :: i,ii,j,k,l,ll
 REAL (kind=8) :: bm(12,3),bnn(12,3,4),bmn(12,3,4),bbn1(3,12,3),bbn2(3,4,3), &
                  bmg(12,4),bng(12,4),aux(3)
 REAL (kind=8) :: ca,c(3),ml,n1,n2,nl,c1,c2,c3
 LOGICAL :: branch

 bb = 0d0                      !zero local B matrix
 bm = 0d0                      !zero local B matrix
 bmn = 0d0                     !zero local B matrix
 bnn = 0d0                     !zero local B matrix
 bbn1 = 0d0                    !zero local B1 matrix
 bbn2 = 0d0
 bbk = 0d0
 bbg = 0d0
 bmg = 0d0
 bng = 0d0
 branch = .FALSE.

 DO i=1,4                      ! for each side
   ii = i+4                          ! position of adjacent element
   IF (sides(i) == 1) THEN           ! side element exists
   ! I - P A R T:  C O N T R I B U T I O N  O F  M A I N  E L E M E N T
       nl = SQRT(ns(1,i)**2+ns(2,i)**2)                      !side length
       n1 = ns(1,i)/nl           !normalize normal vector
       n2 = ns(2,i)/nl
       !     side rotation
       c1 = n1*n2*nl                ! n1 n2
       c2 = (n1**2-n2**2)*nl        ! n1^2-n2^2
       DO j=1,4
         bm(j,1)= bm(j,1) - c1*cn(j,2,i)        !K11
         bm(j,2)= bm(j,2) + c1*cn(j,2,i)        !K22
         bm(j,3)= bm(j,3) + c2*cn(j,2,i)        !2K12
       END DO
       !     normal to the side rotation
       c3 = 2d0*c1    !2 n1 n2
       c1 = n1**2*nl   ! n1^2
       c2 = n2**2*nl   ! n2^2
       DO j=1,4                        ! for each node
         bm(j,1)= bm(j,1) + c1*cn(j,1,i)       !K11
         bm(j,2)= bm(j,2) + c2*cn(j,1,i)       !K22
         bm(j,3)= bm(j,3) + c3*cn(j,1,i)       !2K12
       END DO
   ! II - P A R T:
       ml = ci(i)          !1/2A*ci
       DO j=1,4                        ! for each node
         bmn(j,1,i)= -c1*cn(j,1,i)*ml
         bmn(j,2,i)= -c2*cn(j,1,i)*ml
         bmn(j,3,i)= -c3*cn(j,1,i)*ml

         bmg(j,i)= -cn(j,1,i)*ml

         k = ln(j,i)
         bnn(k,1,i)= +c1*cn(j,1,ii)*ml
         bnn(k,2,i)= +c2*cn(j,1,ii)*ml
         bnn(k,3,i)= +c3*cn(j,1,ii)*ml

         bng(k,i)= +cn(j,1,ii)*ml
       END DO
   ELSE IF(sides(i) == 0) THEN       ! boundary side
   ! I - P A R T:  C O N T R I B U T I O N  O F  M A I N  E L E M E N T
     IF( is(i) )THEN        !clamped side
       nl = SQRT(ns(1,i)**2+ns(2,i)**2)                      !side length
       n1 = ns(1,i)/nl           !normalize normal vector
       n2 = ns(2,i)/nl
       !     side rotation
       c1 = n1*n2*nl                ! n1 n2
       c2 = (n1**2-n2**2)*nl          ! n1^2-n2^2
       DO j=1,4
         bm(j,1)= bm(j,1) - c1*cn(j,2,i)        !K11
         bm(j,2)= bm(j,2) + c1*cn(j,2,i)        !K22
         bm(j,3)= bm(j,3) + c2*cn(j,2,i)        !2K12
       END DO
       IF( solsi(i) )THEN
         !     normal to the side rotation
         c3 = 2d0*c1    !2 n1 n2
         c1 = n1**2*nl   ! n1^2
         c2 = n2**2*nl   ! n2^2
         aux = cn(1:3,2,ii)/area(0)
         !ml = ci(i)          !1/2A*ci
         DO j=1,4                        ! for each node
           bmn(j,1,i)= +c1*cn(j,1,i) !*ml
           bmn(j,2,i)= +c2*cn(j,1,i) !*ml
           bmn(j,3,i)= +c3*cn(j,1,i) !*ml
         END DO

         !DO j=1,4                        ! for each node
         !  bm(j,1)= bm(j,1) + c1*cn(j,1,i)       !K11
         !  bm(j,2)= bm(j,2) + c2*cn(j,1,i)       !K22
         !  bm(j,3)= bm(j,3) + c3*cn(j,1,i)       !2K12
         !END DO

       END IF
     ELSE ! simple supported or free
       nl = SQRT(ns(1,i)**2+ns(2,i)**2)                      !side length
       n1 = ns(1,i)/nl           !normalize normal vector
       n2 = ns(2,i)/nl
       !     side rotation
       c1 = n1*n2*nl                ! n1 n2
       c2 = (n1**2-n2**2)*nl          ! n1^2-n2^2
       DO j=1,4
         bm(j,1)= bm(j,1) - c1*cn(j,2,i)        !K11
         bm(j,2)= bm(j,2) + c1*cn(j,2,i)        !K22
         bm(j,3)= bm(j,3) + c2*cn(j,2,i)        !2K12
       END DO
       !     normal to the side rotation
       c3 = 2d0*c1    !2 n1 n2
       c1 = n1**2*nl   ! n1^2
       c2 = n2**2*nl   ! n2^2
       DO j=1,4                        ! for each node
         bm(j,1)= bm(j,1) + c1*cn(j,1,i)  !K11
         bm(j,2)= bm(j,2) + c2*cn(j,1,i)  !K22
         bm(j,3)= bm(j,3) + c3*cn(j,1,i)  !2K12
       END DO
     END IF
   ELSE !IF(sides(i) < 0 ) THEN       ! branching side
     ! I - P A R T:  C O N T R I B U T I O N  O F  M A I N  E L E M E N T
     branch = .TRUE.
     nl = SQRT(ns(1,i)**2+ns(2,i)**2)                      !side length
     n1 = ns(1,i)/nl           !normalize normal vector
     n2 = ns(2,i)/nl
     !     side rotation
     c1 = n1*n2*nl                ! n1 n2
     c2 = (n1**2-n2**2)*nl        ! n1^2-n2^2
     DO j=1,4
       bm(j,1)= bm(j,1) - c1*cn(j,2,i)        !K11
       bm(j,2)= bm(j,2) + c1*cn(j,2,i)        !K22
       bm(j,3)= bm(j,3) + c2*cn(j,2,i)        !2K12
     END DO
     !     normal to the side rotation
     c3 = 2d0*c1    !2 n1 n2
     c1 = n1**2*nl   ! n1^2
     c2 = n2**2*nl   ! n2^2
     DO j=1,4                        ! for each node
       bm(j,1)= bm(j,1) + c1*cn(j,1,i)      !K11
       bm(j,2)= bm(j,2) + c2*cn(j,1,i)      !K22
       bm(j,3)= bm(j,3) + c3*cn(j,1,i)      !2K12
     END DO
     ! II - P A R T:
     ml = 1d0 !ci(i)          !1/2A*ci
     DO j=1,4                        ! for each node
       bmn(j,1,i)= -c1*cn(j,1,i)*ml
       bmn(j,2,i)= -c2*cn(j,1,i)*ml
       bmn(j,3,i)= -c3*cn(j,1,i)*ml

       bmg(j,i)= -cn(j,1,i)*ml
       k = ln(3,i)
     !signs changed because they will be used to substract average rotation B matrix
       bb(1,k,4)= -c1*ml/area(0)
       bb(2,k,4)= -c2*ml/area(0)
       bb(3,k,4)= -c3*ml/area(0)
     END DO
   END IF
 END DO

 !I-Part:
 c = -t0/area(0)
 DO j=1,3        !for each curvature
   DO i=1,12     !for each node in the patch
     bb(:,i,j) = bm(i,j)*c
   END DO
 END DO

 !II Part:
 ca = -1d0/area(0)
 DO l=1,4            !for each side
   ll = l+4

   !STANDARDT ELEMENT
   IF( sides(l)==1 )THEN
     DO j=1,3        !for each curvature
       DO i=1,12     !for each node in the patch
         bbn1(:,i,j) = bbn1(:,i,j) + bmn(i,j,l)*ca*tb(:,l) + bnn(i,j,l)*ca*tb(:,ll)
       END DO
     END DO
     DO i=1,12
       bbk(:,i,l) = (bmg(i,l)*tb(:,l) +  bng(i,l)*tb(:,ll))/dt(l)
       bbg(:,i,l) = bmg(i,l)*tb(:,l) +  bng(i,l)*tb(:,ll)
     END DO

   !BRANCHING ELEMENT
   ELSE IF( sides(l)==-1 )THEN
     DO j=1,3        !for each curvature
       DO i=1,4      !for each node in the patch
         bbn2(:,i,j) = bbn2(:,i,j) + bmn(i,j,l)*ca*tb(:,l)
       END DO
     END DO
     DO i=1,4
       bbg(:,i,l) = bmg(i,l)*tb(:,l)!/dt(l)
     END DO
   END IF
 END DO

 bb(:,:,1:3) = bb(:,:,1:3) + bbn1(:,:,1:3)
 bb(:,1:4,1:3) = bb(:,1:4,1:3) + bbn2(:,1:4,1:3)

 DO l=1,4
   nl = SQRT(ns(1,l)**2+ns(2,l)**2)      !side length
   n1 = ns(1,l)/nl                       !normalize normal vector
   n2 = ns(2,l)/nl
   !     normal to the side rotation
   c1 = n1**2     ! n1^2
   c2 = n2**2     ! n2^2
   c3 = n1*n2     ! n1 n2
   IF( sides(l)==1 )THEN  !STANDART ELEMENT
     bbk(:,:,l) = -bbk(:,:,l) - (c1*bbn1(:,:,1) + c2*bbn1(:,:,2) + c3*bbn1(:,:,3))
   END IF
 END DO

 RETURN
 END SUBROUTINE bfle25
