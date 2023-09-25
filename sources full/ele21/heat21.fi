SUBROUTINE heat21 (ntype,ngaus,nelem,matno,lnods,weigp,shape,  &
                   deriv,lcur,hgs,nnode,temp,coord,ndime,q,ttime)

  !     calculates element conductivity matrix
  !     for two dimensional planar and axisymmetric geometry

  USE npo_db, ONLY: iftmp
  IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: ngaus,ntype,lcur,nelem,matno(:),nnode,lnods(:,:),ndime
  REAL(kind=8)::  hgs,ttime,weigp(2),temp(:),coord(:,:),q(:), &
                    shape(:,:),deriv(:,:,:)
  !--- Function
  REAL(kind=8):: functs
  !--- Local variables
  INTEGER(kind=4):: j,ii,iel,istop,ndof,i,nd,k,ig,jg,kg,nodm, &
                    isec,osec
  REAL(kind=8):: jacm(ndime,ndime),sk(nnode,nnode),tnod(nnode),  &
                  elcod(ndime,nnode),xkt(nnode),b(ndime,nnode),hf(nnode)
  REAL(kind=8):: wt,detj,fac,rad,tip,ek(3),ek12,ek11,ek22,bt1,bt2,  &
                  db1,db2, qg
  REAL(kind=8),PARAMETER:: twopi=6.283185307179586d0
  TYPE(section),POINTER:: sec

  ndof  = nnode*nnode

  qg = hgs
  IF (lcur > 0) qg = functs(lcur,ttime)*hgs

  osec = -1
  DO iel=1,nelem

    isec = matno(iel)
    IF( isec /= osec )THEN
      sec => psecs(isec)%p
      ek(1) = sec%mtbas%prope(8)           !conductivity
      !  temperature dependent conductivities can be used in the future
      !  isotropic conductivity assumed, but prepared for anisotropy
      ek11 = ek(1)
      ek22 = ek(1)
      ek12 = 0d0

      osec = isec
    END IF

    !     zero element conductivity matrix sk(nnode,nnode)

    sk = 0d0
    hf = 0d0

    !     determine element nodal temperature vector
    !     and element nodal coordinates

    DO i=1,nnode
      nd=lnods(i,iel)
      IF (nd /= 0) THEN
        tnod(i)=temp(nd)
        DO k=1,2
          elcod(k,i)=coord(k,nd)
        END DO
      END IF
    END DO

    !  Loop over all integration points

    kg=0
    DO ig=1,ngaus
      DO jg=1,ngaus
        kg=kg+1

        IF (nnode == 4) wt = weigp(ig)*weigp(jg)
        IF (nnode == 3) wt = weigp(1)

        !  Find jacobian (xj) and its determinant (detj).

        CALL shape1t (elcod,deriv(:,:,kg),jacm,detj,iel,nnode,ndime,istop)

        !  Evaluate jacobian inverse (xji) and global derivative operator (b)
        !  at each integration point (r,s) within the element

        CALL  deriv1 (elcod,shape(:,kg),deriv(:,:,kg),b,jacm,detj,  &
                      rad,ntype,iel,nnode,ndime)

        fac = wt*rad*detj
        IF (ntype == 3) fac = fac * twopi


        !  Compute temperature at integration point

        !IF  (ntp /= 1) THEN
        tip = 0d0
        DO i=1,nnode
          tip = tip + shape(i,kg)*tnod(i)
        END DO
        !END IF

        IF (lcur > 0) THEN
          DO i=1,nnode
            hf(i) = hf(i) + shape(i,kg)*fac
          END DO
        END IF

        !  Form  SK = B(TRANSPOSE)*EK*B for integration point

        DO i=1,nnode
          bt1 = b(1,i)
          bt2 = b(2,i)
          DO j=i,nnode
            db1 = ek11*b(1,j) + ek12*b(2,j)
            db2 = ek12*b(1,j) + ek22*b(2,j)
            sk(i,j) = sk(i,j) + (bt1*db1 + bt2*db2)*fac
          END DO
        END DO

      END DO
    END DO

    nodm = nnode - 1
    DO i=1,nodm
      ii = i + 1
      DO j=ii,nnode
        sk(j,i) = sk(i,j)
      END DO
    END DO

    !     XKT = SK *TNOD

    CALL proma3(xkt,sk,tnod,nnode,nnode,.TRUE.)


    hf = hf * qg

    DO i=1,nnode                          !for each element node
      nd=lnods(i,iel)                     !node number
      IF (nd /= 0) THEN
        k = iftmp(1,nd)                   !DOF
        IF( k > 0 )THEN
          q(k) = q(k) + xkt(i) - hf(i)
        END IF
      END IF
    END DO
  END DO

  RETURN

END SUBROUTINE heat21
