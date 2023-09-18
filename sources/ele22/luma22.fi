      SUBROUTINE luma22 (ngaus,nelem,matno,lnods,weigp, &
     &                   shape,deriv,nnode,coord,ndime,c)

      !     calculates lumped heat capacity
      !     for three dimensional 4/8-node element

      USE npo_db, ONLY : iftmp

      IMPLICIT NONE
      INTEGER(kind=4) :: ndime,nelem,ngaus,nnode,matno(:),lnods(:,:)
      REAL(kind=8) :: weigp(2), shape(:,:), deriv(:,:,:), &
     &                coord(:,:),c(:)

!      INTERFACE
!        INCLUDE 'addc.h'
!      END INTERFACE

      INTEGER (kind=4) :: ig,jg,kg,lg,ng,iel,inode,istop,isec,osec
      REAL (kind=8) :: djacb, dvolu, fac, dens, shapi, sp, sumas, tarea
      REAL(kind=8),ALLOCATABLE :: cartd(:,:),elcod(:,:),sc(:), &
     &                            jacm(:,:)
      TYPE (section), POINTER :: sec

      ng = ngaus*ngaus*ngaus
      ALLOCATE (cartd(nnode,ndime), elcod(3,nnode),  &
     &          sc(nnode), jacm(ndime,ndime))

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

        !  determine element nodal coordinates

        elcod(1:3,1:nnode) = coord(1:3,lnods(1:nnode,iel))

        sc = 0d0
        tarea = 0.0
        kg=0

        DO ig=1,ngaus
          DO jg=1,ngaus
            DO lg=1,ngaus
              kg = kg+1
              CALL jacob5(cartd,deriv(1,1,kg),djacb,elcod,nnode, &
     &                    iel,istop)
              IF(istop == 1)STOP
              IF (nnode == 8) THEN
                dvolu = djacb*weigp(ig)*weigp(jg)*weigp(lg)
                DO inode = 1,nnode
                  shapi = shape(inode,kg)
                  sc(inode) = sc(inode) + shapi*shapi*shapi*dvolu
                END DO
              ELSE IF (nnode == 4) THEN ! is valid only for ngaus = 1
                dvolu = djacb*weigp(ig)
                DO inode = 1,nnode
                  sc(inode) = dvolu/nnode
                END DO
              ENDIF
              tarea = tarea+dvolu
            END DO
          END DO
        END DO

        !  generates lumped heat capacity matrix proportional to diagonal

        sumas = 0.
        DO inode=1,nnode
          sumas = sumas + sc(inode)
        END DO
        fac  = tarea/sumas
        sc   = sc * fac * dens * sp

        !  assemble heat capacity vector (c)

        !call  addc (c,sc,lnods(1:nnode,iel),nnode)
        DO ig=1,nnode                          !for each element node
          ng=lnods(ig,iel)                     !node number
          IF (ng /= 0) THEN
            kg = iftmp(1,ng)                   !DOF
            IF( kg > 0 )c(kg) = c(kg) + sc(ig)
          END IF
        END DO

      END DO

      DEALLOCATE (cartd, elcod, sc, jacm)

      RETURN

      END SUBROUTINE luma22
