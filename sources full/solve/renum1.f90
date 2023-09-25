      SUBROUTINE renum1(nnode,npoin,nelem,lnods,nodad,nposi)
!***********************************************************************
!*
!*****computes the connections between degrees of freedom (stored as
!*    a linked list)
!*
!***********************************************************************
      IMPLICIT none
!     routine parameters
      INTEGER (kind=4),INTENT(IN) :: nnode,npoin,nelem,lnods(:,:)
      INTEGER (kind=4),INTENT(OUT):: nodad(:),nposi
!     local variables
      INTEGER (kind=4) :: ipoin,nfree,ielem,inode,jnode,jpoin,nadre,    &
     &                    nadro,kpoin

      nodad = 0
      nfree = npoin+1

      DO ielem = 1,nelem
        DO inode = 1,nnode
          ipoin = lnods(inode,ielem)
          IF(ipoin /= 0) THEN
            DO jnode = 1,nnode
              IF(jnode /= inode) THEN
                jpoin = lnods(jnode,ielem)
                IF(jpoin /= 0) THEN
                   nadre = nodad(ipoin)
                   nadro = ipoin
                   DO
                     IF(nadre <= 0)EXIT
                     kpoin = nodad(nadre)
                     IF(kpoin == jpoin) EXIT
                     nadro = nadre + 1
                     nadre = nodad(nadro)
                   END DO
                   IF(kpoin /= jpoin)THEN
                     nodad(nadro) = nfree
                     nodad(nfree) = jpoin
                     nodad(nfree+1) = 0
                     nfree = nfree + 2
                   END IF
                 END IF
              END IF
            END DO
          END IF
        END DO
      END DO
      nposi = nfree

      RETURN

      END SUBROUTINE renum1
