 SUBROUTINE heat22(ngaus,nelem,matno,lnods,weigp,shape,deriv,lcur,       &
                   hgs,nnode,temp,coord,ndime,q,ttime)
 !     calculates Thermal residual
 !     for three-dimensional geometry
 USE npo_db, ONLY: iftmp
 IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: ndime,ngaus,lcur,nelem,nnode,matno(:),lnods(:,:)
  REAL(kind=8)::  hgs,ttime,weigp(2),temp(:),coord(:,:),q(:), shape(:,:),deriv(:,:,:)
  !--- Function
  REAL(kind=8):: functs
  !--- Local variables
  INTEGER(kind=4):: i,ig,jg,kg,lg,iel,istop,isec,osec
  REAL(kind=8):: hf(nnode),sk(nnode,nnode),tnod(nnode),elcod(ndime,nnode),   &
                 xkt(nnode),cartd(nnode,ndime)
  REAL(kind=8):: wt,djacb,dvolu,ek11,tip, qg
  TYPE(section),POINTER:: sec

  istop=0 !initializes

  ! internal source
  qg = hgs         !constant value
  IF (lcur > 0) qg = functs (lcur,ttime)*hgs  !curve factor

  osec = -1          !initializes
  DO iel=1,nelem     !for each element

    isec = matno(iel)        !associated section
    IF( isec /= osec )THEN
      sec => psecs(isec)%p
      ! temperature dependent conductivities can be used in the future
      ! isotropic conductivity assumed, but prepared for anisotropy
      ek11 = sec%mtbas%prope(8)           !conductivity
      osec = isec
    END IF

    !     zero element conductivity matrix sk(nnode,nnode)

    sk = 0d0
    hf = 0d0

    tnod (    1:nnode) = temp (    lnods(1:nnode,iel))  !nodal temperatures
    elcod(1:3,1:nnode) = coord(1:3,lnods(1:nnode,iel))  !nodal coordinates

    !  Loop over all integration points

    kg=0  !kg=1,ngaus^3
    DO ig=1,ngaus
      DO jg=1,ngaus
        DO lg=1,ngaus
          kg = kg+1
          ! compute Gauss-point Weigth
          IF (nnode == 8) wt = weigp(ig)*weigp(jg)*weigp(lg)  !brick
          IF (nnode == 4) wt = weigp(1)                       !tetrahedra
          !  Find jacobian, its determinant (djacb),
          !  jacobian inverse and cartesian derivative operator (cartd)
          CALL jacob5(cartd,deriv(1,1,kg),djacb,elcod,nnode,iel,istop)
          dvolu = djacb*wt                                    !

          !  Compute temperature at integration point
          tip = DOT_PRODUCT(tnod(1:nnode),shape(1:nnode,kg))

          IF (qg /= 0d0) THEN
            DO i=1,nnode
              hf(i) = hf(i) + shape(i,kg)*dvolu
            END DO
          END IF

          !  Form  SK = B(TRANSPOSE)*EK*B for integration point

          sk = sk + MATMUL ( cartd, TRANSPOSE(cartd) ) *ek11 * dvolu

        END DO
      END DO
    END DO


    CALL proma3(xkt,sk,tnod,nnode,nnode,.TRUE.)  !     XKT = SK *TNOD

    hf = hf * qg                                 ! internal source

    DO i=1,nnode                          !for each element node
      jg=lnods(i,iel)                     !node number
      IF (jg /= 0) THEN                   ! if node exists
        kg = iftmp(1,jg)                   !DOF
        IF( kg > 0 )THEN                   ! if DOF exists
          q(kg) = q(kg) + xkt(i) - hf(i)
        END IF
      END IF
    END DO
  END DO

RETURN
END SUBROUTINE heat22
