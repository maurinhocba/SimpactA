 SUBROUTINE bran13( nsb, bhead, coora)
 !
 !  computes B matrix for branching sides
 !  B relates average rotation with incremental displacements
 !
 !$ USE omp_lib
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nsb
 REAL (kind=8), INTENT(IN) :: coora(:,:)
 TYPE (sideb), POINTER :: bhead

 INTEGER (kind=4) :: i,j,nn,ni
 REAL (kind=8) :: ls,ang(8),cosa,sina,x(3,3),t(3,9),n(3,8),lb
 TYPE (pside), POINTER :: si(:)  !array of pointers to sides

 !-----------------================------------------
 ALLOCATE(si(nsb))
 si(1)%p => bhead   !point to first branching side
 DO i=2,nsb      !for each branching side
   si(i)%p => si(i-1)%p%next
 END DO
 !$OMP DO  &
 !$OMP PRIVATE(j,nn,ni,ls,ang,cosa,sina,x,t,n,lb)
 DO i=1,nsb    !for each branching side
   nn = si(i)%p%nn  !number of elements joined to the side
   x(:,1:2) = coora(:,si(i)%p%lnods(-1:0)) !branching side coordinates j k
   t(:,9)   = x(:,2) - x(:,1)              !side vector j => k
   CALL vecuni(3,t(1,9),ls)                !unit side vector(t9) and length (ls)
   ls = ls * si(i)%p%alph0(nn)             !side stretching ls/ls(0)
   si(i)%p%bb(:,-1:0) = 0d0                !initializes first two nodes (side nodes) in B matrix
   DO j=1,nn                          !for each connected element
     ni = ABS(si(i)%p%lnods(j))            !extra node number
     x(:,3) = coora(:,ni)                  !coordinates of the node
     n(:,j) = MATMUL( x(:,:),si(i)%p%c(:,j)) !normal derivative
     CALL vecpro(n(1,j),t(1,9),t(1,j))     !normal vector to the plane tn x s = t
     CALL vecuni(3,t(1,j),lb)              !unit normal vector (tj) and det/ls = A/A0/lambda_s
     lb = si(i)%p%fc(j,nn)/lb/ls           !ri * lambda_j (lambda_j = 1/det
     CALL vecpro(t(1,9),t(1,j),n(1,j))     !normal vector to the side n = s x t (outward)
     !     contributions to average rotation B matrix
     si(i)%p%bb(:,-1) = si(i)%p%bb(:,-1) + lb*t(:,j)*si(i)%p%c(1,j)  !node j
     si(i)%p%bb(:, 0) = si(i)%p%bb(:, 0) + lb*t(:,j)*si(i)%p%c(2,j)  !node k
     si(i)%p%bb(:, j) =                  + lb*t(:,j)*si(i)%p%c(3,j)  !node i
   END DO
   !  compute present angles
   DO j=1,nn-1                 !for each consecutive pair (not ordered)
     cosa = DOT_PRODUCT(t(:,j),t(:,j+1))               !t_i . t_i+1
     sina = DOT_PRODUCT(n(:,j),t(:,j+1))               !n_i . t_i+1
     ang(j) = atan4(sina,cosa,si(i)%p%alph0(j),0d0)    !angle change
   END DO
   si(i)%p%gamma = MATMUL(si(i)%p%fc(:,1:nn-1),ang(1:nn-1))      !distorsions
 END DO
 !$OMP END DO
 DEALLOCATE(si)

 RETURN
 END SUBROUTINE bran13
