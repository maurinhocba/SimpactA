SUBROUTINE gaus21(nelem, ngaus, nnode, lnods, coord, posgp, weigp, &
                  shape, deriv, ntype, istop)

!*** evaluates gauss points constants
  IMPLICIT NONE

  INTEGER (kind=4) :: nelem,ngaus,nnode,lnods(:,:),ntype,istop
  REAL    (kind=8) coord(:,:),posgp(2),weigp(2),shape(:,:),deriv(:,:,:)

  INTEGER (kind=4) kg,ig,jg,ielem,gg
  REAL    (kind=8) exisp,etasp,derx(2,nnode),dvolu(ngaus*ngaus), &
                   elcod(2,nnode)

! gauss points in local coordinates and weigths

  IF (nnode == 4) THEN
    CALL gaussq (ngaus ,posgp ,weigp )
  ELSE !IF (nnode == 3) THEN
    !CALL gaustr (ngaus ,posgp ,weigp )
    posgp(1) = 1.0d0/3.0d0
    weigp(1) = 1.0d0/2.0d0
  END IF

! gauss points shape and derivatives of nodal functions

  kg = 0
  DO ig=1,ngaus
    exisp = posgp(ig)
    DO jg=1,ngaus
      kg = kg+1
      etasp = posgp(jg)
      CALL shape3(deriv(1,1,kg),shape(1,kg),exisp,etasp,nnode)
    END DO
  END DO

  gg = 0
  DO ielem=1,nelem
    elcod(1:2,1:nnode) = coord(1:2,lnods(1:nnode,ielem))

    ! verify coordinates for axisymmetric model
    IF( ntype == 3 .AND. ANY(elcod(1,:) < 0d0))THEN
      WRITE(55,"('For axisymmetric problems all X1 (r) coordinates' &
                 ' must be possitive ',/, &
                 ' ELEMENT ',i6,/,(2E15.6) )",ERR=9999) ielem,elcod
      CALL runen3('GAUS21: negative coord. r detected ')
    END IF

    kg=0
    DO ig=1,ngaus
      exisp = posgp(ig)
      DO jg=1,ngaus
        gg = gg+1
        kg = kg+1
        etasp = posgp(jg)
        ! check jacobian
        CALL jacob3(derx,deriv(1,1,kg),dvolu(kg),elcod,ielem,nnode,istop)
      END DO
    END DO

  END DO

RETURN
 9999 CALL runen2('')
END SUBROUTINE gaus21
