 SUBROUTINE bran25( nsb, bhead, coora)
 !
 !  computes B matrix for branching sides
 !  B relates average rotation with incremental displacements
 !
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nsb
 REAL (kind=8), INTENT(IN) :: coora(:,:)
 TYPE (sideb), POINTER :: bhead

 INTEGER (kind=4) :: i,j,k,nn,ni,mi
 REAL (kind=8) :: ls,ang(8),sina,cosa,x(3,4),t(3,9),ti(3,8),n(3,8),lb
 TYPE (sideb), POINTER :: si

 !-----------------================------------------

 si => bhead   !point to first branching side
 DO i=1,nsb    !for each branching side
   nn = si%nn  !number of elements joined to the side
   x(:,1:2) = coora(:,si%lnods(-1:0)) !branching side coordinates j k
   t(:,9)   = x(:,2) - x(:,1)         !side vector j ==> k  (s)
   CALL vecuni(3,t(1,9),ls)           !unit side vector(t9) and length (ls)
   ls = ls * si%alph0(nn)             !side stretching ls/ls(0)
   si%bbr(:,-1:0) = 0d0                !initializes first two nodes (side nodes) in B matrix
   k = 0
   DO j=1,nn                          !for each connected element
     ni = ABS(si%lnods(j+k))            !extra node number
     mi = ABS(si%lnods(j+k+1))
     x(:,3) = coora(:,ni)             !coordinates of the first node
     x(:,4) = coora(:,mi)             !coordinates of the second node
     n(:,j) = MATMUL(x(:,1:4),si%c(:,j,nn))       !normal derivative
     CALL vecpro(n(1,j),t(1,9),ti(1,j))  !normal vector to the side tn x s = t
     CALL vecuni(3,ti(1,j),lb)         !unit normal vector (tj) and det/ls = A/A0/lambda_s
     CALL vecpro(t(1,9),ti(1,j),n(1,j))!normal vector to the side n = s x t (outward)
     lb = si%fc(j,nn)/lb     !stiffness / heigth
     !     contributions to average rotation B matrix
     si%bbr(:,-1)   = si%bbr(:,-1) + lb*ti(:,j)*si%c(1,j,nn)  !node j
     si%bbr(:, 0)   = si%bbr(:, 0) + lb*ti(:,j)*si%c(2,j,nn)  !node i
     IF (j == 1) THEN !first element
       si%bbr(:, j)     =            + lb*ti(:,j)*si%c(3,j,nn)  !node k
       si%bbr(:, j+1)   =            + lb*ti(:,j)*si%c(4,j,nn)  !node l
     ELSE             !second element upto the last element
       si%bbr(:, j+k)   =            + lb*ti(:,j)*si%c(3,j,nn)  !node k
       si%bbr(:, j+k+1) =            + lb*ti(:,j)*si%c(4,j,nn)  !node l
     END IF

     k = k + 1
   END DO
   !  compute present angles
   DO j=1,nn-1                 !for each consecutive pair (not ordered)
     !cosa = -DOT_PRODUCT(n(:,j),n(:,j+1)) !n_i . n_i+1
     !sina = +DOT_PRODUCT(n(:,j+1),ti(:,j))!n_i+1 . t_i
     !ang(j) = atan4(sina,cosa,si%alph0(j),0d0)         !angle change

     cosa = +DOT_PRODUCT(n(:,j),n(:,j+1)) !n_i . n_i+1
     sina = +DOT_PRODUCT(n(:,j+1),ti(:,j))!n_i+1 . t_i
     ang(j) = atan4(sina,cosa,si%alph0(j),0d0)         !angle change
   END DO
   si%gamma = MATMUL(si%fc(:,1:nn-1),ang(1:nn-1))      !distorsions
   si => si%next      !point to next side
 END DO

 RETURN
 END SUBROUTINE bran25
