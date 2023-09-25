 SUBROUTINE jacob05(deriv,cds,dvol,x,ielem,istop, &
                    cdq,lb,sides,nfdas,jacin,quad)
 !*********************************************************************
 !
 !**** this SUBROUTINE evaluates the jacobian matrix and the cartesian
 !     shape FUNCTION derivatives
 !
 !*********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT(IN) :: x(:,:),        & !(3,nnode)patch coordinates
                             lb(:,:),       & !(3,3) Local cartesian system at element center
                             deriv(:,:,:)     !(6,3,3)Local natural derivatives at faces and center
 REAL(kind=8), INTENT(OUT) :: dvol(:),      & !(ngaud=3) gauss point volume at faces and center
                              cds(:),       & !(nnb)   standard cartesian derivative (y3) at center
                              cdq(:,:,:,:), & !(4,2,3,2) cartesian derivatives for in-plane
                              nfdas(:,:,:), & !(nnb,3,2) N_3 at gauss points for assumed shear
                              jacin(:,:,:)    !(2,2,2) in-plane inverse jacobian at gauss points
 LOGICAL, INTENT(IN) :: quad,sides(3)         !if QUAD aproach and if side elements exist
 INTEGER(kind=4) :: ielem,istop               !in case of error
 ! local variables
 INTEGER (kind=4) ii,jj,i,j,k,g
 REAL (kind=8) l1(3),l2(3),t1(3),t2(3),t3(3), &
               lt,ji(2,2),jac(2,2),a(3),b(3),aux,j0, &
               shape(nnb),der(nnb,3),jaci(3,3),cd(nnb,3),gcova(3,3)
 ! natural derivatives (2) of the 6-node patch shape functions at (3) mid side points (QUAD)
 REAL (kind=8), PARAMETER :: nd(nnb,2,3) =          &
   (/ -0.5d0, 0.5d0,-0.5d0, 0.5d0, 0.0d0, 0.0d0 ,   &
      -0.5d0,-0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0 ,   &
      -0.5d0, 0.5d0, 0.5d0, 0.0d0,-0.5d0, 0.0d0 ,   &
      -1.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ,   &
      -1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ,   &
      -0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0,-0.5d0 /)
 ! mid-side gauss points position (SHELL)
 REAL (kind=8), PARAMETER :: gpa(2,3) =   (/ 0.5d0, 0.5d0, 0.0d0, 0.5d0, 0.5d0, 0.0d0  /)
 REAL (kind=8), PARAMETER :: gpz(2)   =   (/ -1d0, 1d0  /)

 DO jj=1,2  !for each face
    i = 3*jj-2      !point to first node in the face (1 or 4)

    !   in-plane behavior

    IF( quad ) THEN   !for quadratic approach
      k = i + nnb       !point to fourth node in the face
      DO ii=1,3       !for each side
        IF( sides(ii) )THEN                   !adjacent element exist
          ! compute cartesian derivatives for membrane part
          l1 =      MATMUL(x(:,i:i+2),nd(1:3,1,ii))   ! natural derivative (xi)   first 3 nodes
          l2 =      MATMUL(x(:,i:i+2),nd(1:3,2,ii))   ! natural derivative (eta)  first 3 nodes
          l1 = l1 + MATMUL(x(:,k:k+2),nd(4:6,1,ii))   ! natural derivative (xi)   add extra node
          l2 = l2 + MATMUL(x(:,k:k+2),nd(4:6,2,ii))   ! natural derivative (eta)  add extra node

          CALL vecpro(l1,l2,t3)          ! normal direction at gauss point
          CALL vecuni(3,t3,lt)           ! unit vector
          CALL vecpro(t3,lb(1,1),t2)     !  t2 direction at face side
          CALL vecuni(3,t2,lt)           !     normalizes t2
          CALL vecpro(t2,t3,t1)          !  t1 direction at face side

          !     Jacobian Matrix
          jac(1,1) = DOT_PRODUCT(t1,l1)  !x1,xi
          jac(1,2) = DOT_PRODUCT(t1,l2)  !x1,eta
          jac(2,1) = DOT_PRODUCT(t2,l1)  !x2,xi
          jac(2,2) = DOT_PRODUCT(t2,l2)  !x2,eta
          !     Jacobian Inverse
          lt = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)  !determinant
          ji(1,1) = jac(2,2)/lt         !xi,x1
          ji(1,2) =-jac(1,2)/lt         !xi,x2
          ji(2,1) =-jac(2,1)/lt         !eta,x1
          ji(2,2) = jac(1,1)/lt         !eta,x2
          !compute cartesian derivatives
          DO j=1,3             !main nodes
            cdq(j,1,ii,jj) = ji(1,1)*nd(j,1,ii) + ji(2,1)*nd(j,2,ii)   !N(j),x1
            cdq(j,2,ii,jj) = ji(1,2)*nd(j,1,ii) + ji(2,2)*nd(j,2,ii)   !N(j),x2
          END DO
          j = ii+3              !extra node
          cdq(4,1,ii,jj) = ji(1,1)*nd(j,1,ii) + ji(2,1)*nd(j,2,ii)     !N(4),x1
          cdq(4,2,ii,jj) = ji(1,2)*nd(j,1,ii) + ji(2,2)*nd(j,2,ii)     !N(4),x2
        ELSE
          !*** evaluate the side vectors

          l1 = x(1:3,i+2) - x(1:3,i+1)                             !side 1
          l2 = x(1:3,i  ) - x(1:3,i+2)                             !side 2

          !*** evaluate the cross product => plane normal

          CALL vecpro(l1,l2,t3)       !  t3 direction
          CALL vecuni(3,t3,lt)        !     normalizes t3 lt=twice the area

          CALL vecpro(t3,lb(1,1),t2)        !  t2 direction at face center
          CALL vecuni(3,t2,aux)         !     normalizes t1
          CALL vecpro(t2,t3,t1)         !  t1 direction at face center

          !*** find the local derivatives por central triangle

          a(1) = DOT_PRODUCT(l1,t1) /lt      ! l1 . t1
          a(2) = DOT_PRODUCT(l2,t1) /lt      ! l2 . t1
          a(3) = -a(1)-a(2)
          b(1) = DOT_PRODUCT(l1,t2) /lt      ! l1 . t2
          b(2) = DOT_PRODUCT(l2,t2) /lt      ! l2 . t2
          b(3) = -b(1)-b(2)
          cdq(1:3,1,ii,jj) = -b               !N(j),x1
          cdq(1:3,2,ii,jj) = a                !N(j),x2
          cdq(4,1:2,ii,jj) = 0d0
        END IF
      END DO
    ELSE    !standard Constant Strain Triangle
      !*** evaluate the side vectors

      l1 = x(1:3,i+2) - x(1:3,i+1)                             !side 1
      l2 = x(1:3,i  ) - x(1:3,i+2)                             !side 2

      !*** evaluate the cross product => plane normal

      CALL vecpro(l1,l2,t3)       !  t3 direction
      CALL vecuni(3,t3,lt)        !     normalizes t3 lt=twice the area

      CALL vecpro(t3,lb(1,1),t2)        !  t2 direction at face center
      CALL vecuni(3,t2,aux)         !     normalizes t1
      CALL vecpro(t2,t3,t1)         !  t1 direction at face center

      !*** find the local derivatives por central triangle

      a(1) = DOT_PRODUCT(l1,t1) /lt      ! l1 . t1
      a(2) = DOT_PRODUCT(l2,t1) /lt      ! l2 . t1
      a(3) = -a(1)-a(2)
      b(1) = DOT_PRODUCT(l1,t2) /lt      ! l1 . t2
      b(2) = DOT_PRODUCT(l2,t2) /lt      ! l2 . t2
      b(3) = -b(1)-b(2)
      cdq(1:3,1,1,jj) = -b               !N(j),x1
      cdq(1:3,2,1,jj) = a                !N(j),x2
    END IF

    !  compute derivatives at assumed strain points for transverse shear

    DO g=1,3   !shear assumed strain points
      CALL shape4(nnb,shape,der,gpa(1,g),gpa(2,g),gpz(jj)) !shape function derivatives at mid-surface

      gcova = MATMUL(x(:,1:nnb),der)  ! covariant base

      j0 =  gcova(1,1)*gcova(2,2)*gcova(3,3) +gcova(1,3)*gcova(2,1)*gcova(3,2)  &
           +gcova(3,1)*gcova(1,2)*gcova(2,3) -gcova(3,1)*gcova(2,2)*gcova(1,3)  &
           -gcova(3,3)*gcova(1,2)*gcova(2,1) -gcova(1,1)*gcova(2,3)*gcova(3,2)

      IF( j0   <= 0d0 ) THEN
         WRITE(*,600,ERR=9999) ielem
         WRITE(*,900,ERR=9999)
         WRITE(*,910,ERR=9999) x(1:3,1:nnb)
         WRITE(lures,600,ERR=9999) ielem
         WRITE(lures,900,ERR=9999)
         WRITE(lures,910,ERR=9999) x(1:3,1:nnb)
         istop = 1
         RETURN
      END IF
      ! determine a cartesyan local system (y coord) based on center system (lb)
      ! here gcova(:,3) must be set as direction 3
      CALL vecpro(gcova(1,1),gcova(1,2),jaci(1,3)) !tl_3 = g_xita x g_eta
      CALL vecuni(3,jaci(1,3),aux)                 !unit vector, normal to xita-eta plane
      CALL vecpro(jaci(1,3),lb(1,1),jaci(1,2))     !tl_2 = tl_3 x lb_1
      CALL vecuni(3,jaci(1,2),aux)                 !unit vector
      CALL vecpro(jaci(1,2),jaci(1,3),jaci(1,1))   !tl_1 = tl_2 x tl_3
      !compute components of gcova respect to the local system
      gcova = MATMUL(TRANSPOSE(jaci),gcova)        !d_y/d_xita
      !WRITE(58,"(i5,2i2,9f10.6)")ielem,jj,g,jaci
      ! compute inverse jacobian (third column) d_(xita,eta,zeta)/d_y3
      t3(1) =  (gcova(1,2)*gcova(2,3)-gcova(1,3)*gcova(2,2))/j0
      t3(2) = -(gcova(1,1)*gcova(2,3)-gcova(1,3)*gcova(2,1))/j0
      t3(3) =  (gcova(1,1)*gcova(2,2)-gcova(2,1)*gcova(1,2))/j0
      ! compute cartesyan (y3) derivatives of the shape functions
      ! necessary to compute f_3
      DO i=1,nnb
        nfdas(i,g,jj) = DOT_PRODUCT(der(i,:),t3)
      END DO

    END DO
 END DO

 !  standard cartesian derivatives at element center

 DO g=1,3  ! for each gauss point
    gcova = MATMUL(x(:,1:nnb),deriv(:,:,g)) !create jacobian matrix gcova
    !     calculate determinant
    j0 =  gcova(1,1)*gcova(2,2)*gcova(3,3) + gcova(1,3)*gcova(2,1)*gcova(3,2)  &
         +gcova(3,1)*gcova(1,2)*gcova(2,3) - gcova(3,1)*gcova(2,2)*gcova(1,3)  &
         -gcova(3,3)*gcova(1,2)*gcova(2,1) - gcova(1,1)*gcova(2,3)*gcova(3,2)
    IF(j0 <= 0) THEN                        !check element is not distorted
      WRITE(*,600,ERR=9999) ielem
      WRITE(*,900,ERR=9999)
      WRITE(*,910,ERR=9999) x(1:3,1:nnb)
      WRITE(lures,600,ERR=9999) ielem
      WRITE(lures,900,ERR=9999)
      WRITE(lures,910,ERR=9999) x(1:3,1:nnb)
      istop = 1
      RETURN
    END IF
    dvol(g) = j0      !keep gauss point volume

    ! and inverse of jacobian matrix
    jaci(1,1) =  (gcova(2,2)*gcova(3,3)-gcova(2,3)*gcova(3,2))/j0
    jaci(2,1) = -(gcova(2,1)*gcova(3,3)-gcova(3,1)*gcova(2,3))/j0
    jaci(3,1) =  (gcova(2,1)*gcova(3,2)-gcova(2,2)*gcova(3,1))/j0
    jaci(1,2) = -(gcova(1,2)*gcova(3,3)-gcova(1,3)*gcova(3,2))/j0
    jaci(2,2) =  (gcova(1,1)*gcova(3,3)-gcova(3,1)*gcova(1,3))/j0
    jaci(3,2) = -(gcova(1,1)*gcova(3,2)-gcova(1,2)*gcova(3,1))/j0
    jaci(1,3) =  (gcova(1,2)*gcova(2,3)-gcova(1,3)*gcova(2,2))/j0
    jaci(2,3) = -(gcova(1,1)*gcova(2,3)-gcova(1,3)*gcova(2,1))/j0
    jaci(3,3) =  (gcova(1,1)*gcova(2,2)-gcova(2,1)*gcova(1,2))/j0
    IF( g == 2 )THEN
      ! compute cartesyan (y3) derivatives of the shape functions
      ! necessary to compute f_3
      cd = MATMUL(deriv(:,:,g),jaci)             !  calculate cartesian derivatives
      DO i=1,nnb
        cds(i) = DOT_PRODUCT(cd(i,:),lb(:,3))
      END DO
    ELSE
      !  keep in plane inverse jacobian
      i = g
      IF( i == 3 ) i = 2
      jacin(1,1,i) = DOT_PRODUCT(jaci(1,:),lb(:,1))  !dxi/dx1
      jacin(2,1,i) = DOT_PRODUCT(jaci(2,:),lb(:,1))  !deta/dx1
      jacin(1,2,i) = DOT_PRODUCT(jaci(1,:),lb(:,2))  !dxi/dx2
      jacin(2,2,i) = DOT_PRODUCT(jaci(2,:),lb(:,2))  !deta/dx2
    END IF
 END DO

 RETURN
  600 FORMAT(//,' PROGRAM halted in SUBROUTINE jacob05',/,11x,           &
     &' zero or negative volume',/,10x,' element number ',i5)
  900 FORMAT(//,5x,'coordinates of element nodes')
  910 FORMAT(5x,3e15.8)
 9999 CALL runen2('')
 END SUBROUTINE jacob05
