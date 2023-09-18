SUBROUTINE luma21 (ntype,ngaus,nelem,matno,lnods,weigp, &
                   shape,deriv,nnode,coord,ndime,c)

  !     calculates lumped heat capacity
  !     for two dimensional planar and axisymmetric geometry

  USE npo_db, ONLY : iftmp

  IMPLICIT NONE
  INTEGER(kind=4) :: ntype,ndime,nelem,ngaus,nnode, &
                     matno(:),lnods(:,:)
  REAL(kind=8) :: weigp(2), shape(:,:), deriv(:,:,:), &
                  coord(:,:),c(:)

  INTEGER (kind=4) :: i,k,ig,jg,kg,iel,nd,istop,isec,osec
  REAL (kind=8) :: wt, detj, fac, rad, dens, sp,  &
                   jacm(ndime,ndime),sc(nnode),elcod(ndime,nnode)
  REAL (kind=8),PARAMETER :: twopi=6.283185307179586
  TYPE (section), POINTER :: sec

  istop=0

  osec = -1
  DO iel=1,nelem

    isec = matno(iel)
    IF( isec /= osec )THEN
      sec => psecs(isec)%p
      dens = sec%mtbas%prope(5)   !material density (not scaled)
      sp   = sec%mtbas%prope(7)   !material capacity
      osec = isec
    END IF

    !  determine element nodal temperature vector
    !  and element nodal coordinates

    DO i=1,nnode
      nd=lnods(i,iel)
      IF (nd /= 0) THEN
        DO k=1,2
          elcod(k,i)=coord(k,nd)
        END DO
      END IF
    END DO

    sc = 0d0

    kg=0
    DO ig=1,ngaus
      DO jg=1,ngaus
        kg=kg+1

        IF (nnode == 4) wt = weigp(ig)*weigp(jg)
        IF (nnode == 3) wt = weigp(1)

        !  Find jacobian (xj) and its determinant (detj).

        CALL shape1t (elcod,deriv(:,:,kg),jacm,detj,iel,nnode,ndime,istop)

        rad = 1.0d0
        IF (ntype == 3) THEN

          !  Compute the radius at point (r,s) for axisymmetric solid

          rad = 0.0d0
          DO i=1,nnode
            rad = rad + shape(i,kg)* elcod(1,i)
          END DO
        END IF
        fac = wt*rad*detj
        IF (ntype == 3) fac = fac * twopi
        !  this type of lumped matrix has been chosen to be consistent with
        !  integration scheme of heat generated in the element
        sc(1:nnode) = sc(1:nnode) + shape(1:nnode,kg)*fac

      END DO
    END DO

    sc(1:nnode) = sc(1:nnode) * dens * sp

    !  assemble heat capacity vector (c)

     DO i=1,nnode                          !for each element node
       nd=lnods(i,iel)                     !node number
       IF (nd /= 0) THEN
         k = iftmp(1,nd)                   !DOF
         IF( k > 0 )c(k) = c(k) + sc(i)
       END IF
     END DO

  END DO

END SUBROUTINE luma21

