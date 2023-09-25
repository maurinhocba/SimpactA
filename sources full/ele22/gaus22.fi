      SUBROUTINE gaus22(nelem, ngaus, nnode, lnods, coord,         &
     &                  posgp, weigp, shape, deriv, istop)
!********************************************************************
!
!*** evaluates gauss points constants for 3D 4/8 node element
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4) :: nelem,ngaus,nnode,lnods(:,:),istop
      REAL    (kind=8) coord(:,:),posgp(2),weigp(2),         &
     &                 shape(:,:),deriv(:,:,:)

      INTEGER (kind=4) kg,ig,jg,lg,elem,ng
      REAL    (kind=8) xita,eta,zeta,jac,         &
     &                 cartd(nnode,3,ngaus*ngaus*ngaus),         &
     &                 elcod(3,nnode)

!     gauss points in local coordinates and weigths

      ng = ngaus*ngaus*ngaus
      IF (nnode == 8) THEN
        CALL gaussq (ngaus ,posgp ,weigp )
      ELSE !IF (nnode == 4) THEN
        CALL gauste (ngaus ,posgp ,weigp )
      END IF

!     gauss points shape and derivatives of nodal functions

      kg=0
      DO lg=1,ngaus
        zeta = posgp(lg)
        DO jg=1,ngaus
          eta = posgp(jg)
          DO ig=1,ngaus
            xita = posgp(ig)
            kg = kg+1
            CALL shape5(deriv(1,1,kg),shape(1,kg),xita,eta,zeta,nnode)
          END DO
        END DO
      END DO

      !   check jacobian

      DO elem=1,nelem
        elcod(1:3,1:nnode) = coord(1:3,lnods(1:nnode,elem))
        DO kg = 1,ng
          CALL jacob5(cartd(1,1,kg),deriv(1,1,kg),jac,         &
     &                elcod,nnode,elem,istop)
          IF(istop == 1) STOP

        END DO
      END DO

      END SUBROUTINE gaus22
