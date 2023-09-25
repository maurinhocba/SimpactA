 SUBROUTINE jacob04(cds,dvol,x,ielem,istop,cdq,nfdas,lb,jacin,area)
 !*********************************************************************
 !
 !**** this SUBROUTINE evaluates the jacobian matrix and the cartesian
 !     shape FUNCTION derivatives
 !
 !*********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT(IN) :: x(:,:),        & !(3,nnode)patch coordinates
                             lb(:,:)          !(3,3) Local cartesian system at element center
 REAL(kind=8), INTENT(OUT) :: dvol(:),      & !(ngaud=3) gauss point volume at faces and center
                              cds(:),       & !(nnode) standard cartesian derivative (y3) at center
                              cdq(:,:,:,:), & !(4,2,4,2) cartesian derivatives for in-plane
                              nfdas(:,:,:), & !(8,4,2) cartesian derivatives y3 direction
                              jacin(:,:,:), & !(2,2,2) in-plane inverse jacobian at gauss points
                              area(:,:)
 INTEGER(kind=4) :: ielem,istop               !in case of error
 ! local variables
 INTEGER (kind=4) ii,jj,i,j
 REAL (kind=8) jaci(3,3),cd(nnode,3),gcova(3,3),t(3,3),qk(3), &
               ca,sa,lt,qn(3),aux,ji(2,2),j3(3),j0,area0

 ! natural derivatives (3) of the 8-node shape functions at (4) mid side points at two faces
 REAL (kind=8), PARAMETER :: nd(nnode,3,0:4,2) =          &
    (/ -0.250D0,  0.250D0,  0.250D0, -0.250D0,  0.000D0,  0.000D0,  0.000D0,  0.000D0, &
       -0.250D0, -0.250D0,  0.250D0,  0.250D0,  0.000D0,  0.000D0,  0.000D0,  0.000D0, &
       -0.125D0, -0.125D0, -0.125D0, -0.125D0,  0.125D0,  0.125D0,  0.125D0,  0.125D0, &
       -0.25D0,  0.25D0,  0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
       -0.50D0,  0.00D0,  0.00D0,  0.50D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
       -0.25D0,  0.00D0,  0.00D0, -0.25D0,  0.25D0,  0.00D0,  0.00D0,  0.25D0,   &
       -0.50D0,  0.50D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
       -0.25D0, -0.25D0,  0.25D0,  0.25D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
       -0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.25D0,  0.25D0,  0.00D0,  0.00D0,   &
       -0.25D0,  0.25D0,  0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
        0.00D0, -0.50D0,  0.50D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
        0.00D0, -0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.25D0,  0.25D0,  0.00D0,   &
        0.00D0,  0.00D0,  0.50D0, -0.50D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
       -0.25D0, -0.25D0,  0.25D0,  0.25D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,   &
        0.00D0,  0.00D0, -0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.25D0,  0.25D0,   &
        0.000D0,  0.000D0,  0.000D0,  0.000D0, -0.250D0,  0.250D0,  0.250D0, -0.250D0, &
        0.000D0,  0.000D0,  0.000D0,  0.000D0, -0.250D0, -0.250D0,  0.250D0,  0.250D0, &
       -0.125D0, -0.125D0, -0.125D0, -0.125D0,  0.125D0,  0.125D0,  0.125D0,  0.125D0, &
        0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.25D0,  0.25D0,  0.25D0, -0.25D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.50D0,  0.00D0,  0.00D0,  0.50D0,   &
       -0.25D0,  0.00D0,  0.00D0, -0.25D0,  0.25D0,  0.00D0,  0.00D0,  0.25D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.50D0,  0.50D0,  0.00D0,  0.00D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.25D0, -0.25D0,  0.25D0,  0.25D0,   &
       -0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.25D0,  0.25D0,  0.00D0,  0.00D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.25D0,  0.25D0,  0.25D0, -0.25D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.50D0,  0.50D0,  0.00D0,   &
        0.00D0, -0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.25D0,  0.25D0,  0.00D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.50D0, -0.50D0,   &
        0.00D0,  0.00D0,  0.00D0,  0.00D0, -0.25D0, -0.25D0,  0.25D0,  0.25D0,   &
        0.00D0,  0.00D0, -0.25D0, -0.25D0,  0.00D0,  0.00D0,  0.25D0,  0.25D0  /),     &
   ndc(8,3) =                                                                          &
    (/ -0.125D0,  0.125D0,  0.125D0, -0.125D0, -0.125D0,  0.125D0,  0.125D0, -0.125D0, &
       -0.125D0, -0.125D0,  0.125D0,  0.125D0, -0.125D0, -0.125D0,  0.125D0,  0.125D0, &
       -0.125D0, -0.125D0, -0.125D0, -0.125D0,  0.125D0,  0.125D0,  0.125D0,  0.125D0 /)

 INTERFACE
   INCLUDE 'vecpro.h'
   INCLUDE 'vecuni.h'
 END INTERFACE

 DO jj=1,2  !for each face
    i = 4*jj-3      !point to first node in the face (1 or 5)
    DO ii=0,4       !for each mid-point side
      gcova = MATMUL(x,nd(:,:,ii,jj))  !covariant base
      ! compute local system
      CALL vecpro(gcova(:,1),gcova(:,2),t(:,3)) !normal
      CALL vecuni(3,t(:,3),aux)             !unit-normal and jacobian
      CALL vecpro(lb(:,3),t(:,3),qk)        !intersection
      CALL vecuni(3,qk,lt)                  !unit direction
      IF( lt > 0d0 )THEN                  ! t3 /= t3c
        ca = DOT_PRODUCT(qk,lb(:,1))
        sa = DOT_PRODUCT(qk,lb(:,2))
        CALL vecpro(t(:,3),qk,qn)           !normal to intersection on tangent plane
        t(:,1) = ca*qk - sa*qn          !first direction on tangent plane at side
        t(:,2) = sa*qk + ca*qn          !second direction on tangent plane at side
      ELSE                                ! t3 = t3c
        t(:,1) = lb(:,1)
        t(:,2) = lb(:,2)
      END IF
      ! rotate gcova to local system t1-t2-t3
      gcova = MATMUL(TRANSPOSE(t),gcova)
      ! compute 3D jacobian
      j0 =  aux*gcova(3,3)               !(3,1)=(3,2)=0
      !j0 =  gcova(1,1)*gcova(2,2)*gcova(3,3) +gcova(1,3)*gcova(2,1)*gcova(3,2)  &
      !     +gcova(3,1)*gcova(1,2)*gcova(2,3) -gcova(3,1)*gcova(2,2)*gcova(1,3)  &
      !     -gcova(3,3)*gcova(1,2)*gcova(2,1) -gcova(1,1)*gcova(2,3)*gcova(3,2)
      IF( j0   <= 0d0 ) THEN !check
         WRITE(*,600,ERR=9999) ielem
         WRITE(*,900,ERR=9999)
         WRITE(*,910,ERR=9999) x(1:3,1:nnode)
         WRITE(lures,600,ERR=9999) ielem
         WRITE(lures,900,ERR=9999)
         WRITE(lures,910,ERR=9999) x(1:3,1:nnode)
         istop = 1
         RETURN
      END IF
      ! in plane jacobian for in-plane strains
      !jac(1,1) = DOT_PRODUCT(t(:,1),gcova(:,1))   !x1,xi
      !jac(1,2) = DOT_PRODUCT(t(:,1),gcova(:,2))   !x1,eta
      !jac(2,1) = DOT_PRODUCT(t(:,2),gcova(:,1))   !x2,xi
      !jac(2,2) = DOT_PRODUCT(t(:,2),gcova(:,2))   !x2,eta
      !area(ii,jj) = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1) !=aux
      !ji(1,1) = jac(2,2)/area(ii,jj)
      !ji(1,2) =-jac(1,2)/area(ii,jj)
      !ji(2,1) =-jac(2,1)/area(ii,jj)
      !ji(2,2) = jac(1,1)/area(ii,jj)
      ji(1,1) = gcova(2,2)/aux
      ji(1,2) =-gcova(1,2)/aux
      ji(2,1) =-gcova(2,1)/aux
      ji(2,2) = gcova(1,1)/aux
      IF( ii == 0 )THEN
        jacin(:,:,jj) = ji
        dvol(jj) = j0
      ELSE
        area(ii,jj) =aux
        ! in-plane cartesian derivatives
        cdq(:,:,ii,jj) = MATMUL(nd(i:i+3,1:2,ii,jj),ji)
        ! transverse inverse jacobian
        j3(1) =  (gcova(1,2)*gcova(2,3)-gcova(1,3)*gcova(2,2))/j0
        j3(2) = -(gcova(1,1)*gcova(2,3)-gcova(1,3)*gcova(2,1))/j0
        j3(3) =  (gcova(1,1)*gcova(2,2)-gcova(2,1)*gcova(1,2))/j0
        DO j=1,nnode
          nfdas(j,ii,jj) = DOT_PRODUCT(nd(j,:,ii,jj),j3)
        END DO
      END IF
    END DO
    area0 = SUM(area(1:4,jj)) !total face area
    area(:,jj) = area(:,jj)/area0
 END DO

 dvol(3) = dvol(2)  !move to third position
 !  standard cartesian derivatives at element center
 gcova = MATMUL(x,ndc(:,:))  !covariant base
 !     calculate determinant
 j0 =  gcova(1,1)*gcova(2,2)*gcova(3,3) +gcova(1,3)*gcova(2,1)*gcova(3,2)  &
      +gcova(3,1)*gcova(1,2)*gcova(2,3) -gcova(3,1)*gcova(2,2)*gcova(1,3)  &
      -gcova(3,3)*gcova(1,2)*gcova(2,1) -gcova(1,1)*gcova(2,3)*gcova(3,2)

 IF( j0   <= 0d0 ) THEN
    WRITE(*,600,ERR=9999) ielem
    WRITE(*,900,ERR=9999)
    WRITE(*,910,ERR=9999) x(1:3,1:nnode)
    WRITE(lures,600,ERR=9999) ielem
    WRITE(lures,900,ERR=9999)
    WRITE(lures,910,ERR=9999) x(1:3,1:nnode)
    istop = 1
    RETURN
 END IF

 dvol(2) = j0      !keep gauss point volume at element center

 ! compute inverse of jacobian matrix
 jaci(1,1) =  (gcova(2,2)*gcova(3,3)-gcova(2,3)*gcova(3,2))/j0
 jaci(2,1) = -(gcova(2,1)*gcova(3,3)-gcova(3,1)*gcova(2,3))/j0
 jaci(3,1) =  (gcova(2,1)*gcova(3,2)-gcova(2,2)*gcova(3,1))/j0
 jaci(1,2) = -(gcova(1,2)*gcova(3,3)-gcova(1,3)*gcova(3,2))/j0
 jaci(2,2) =  (gcova(1,1)*gcova(3,3)-gcova(3,1)*gcova(1,3))/j0
 jaci(3,2) = -(gcova(1,1)*gcova(3,2)-gcova(1,2)*gcova(3,1))/j0
 jaci(1,3) =  (gcova(1,2)*gcova(2,3)-gcova(1,3)*gcova(2,2))/j0
 jaci(2,3) = -(gcova(1,1)*gcova(2,3)-gcova(1,3)*gcova(2,1))/j0
 jaci(3,3) =  (gcova(1,1)*gcova(2,2)-gcova(2,1)*gcova(1,2))/j0

 cd = MATMUL(ndc,jaci)             !  calculate cartesian derivatives
 ! compute cartesyan (y3) derivatives of the shape functions
 ! necessary to compute f_3
 DO i=1,nnode
   cds(i) = DOT_PRODUCT(cd(i,:),lb(:,3))
 END DO

 RETURN
  600 FORMAT(//,' PROGRAM halted in SUBROUTINE jacob04',/,11x,           &
     &' zero or negative volume',/,10x,' element number ',i5)
  900 FORMAT(//,5x,'coordinates of element nodes')
  910 FORMAT(5x,3e15.8)
 9999 CALL runen2('')
 END SUBROUTINE jacob04
