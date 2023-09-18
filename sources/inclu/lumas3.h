      SUBROUTINE lumas3 (npoin,ndofn,nelem,ntype,ngaus,nnode,lnods, &
     &                   matno,coord,emass,posgp,weigp, &
     &                   shape,deriv,iwrit,sumat,ispli)
!******************************************************************
!
! *** calculates lumped mass for 3/4 nodes element
!
!******************************************************************
      USE split_db
      IMPLICIT NONE
      INTEGER (kind=4) nelem,npoin,ndofn,ntype,ngaus,nnode,iwrit, &
     &                 lnods(:,:),matno(:)
      INTEGER (kind=4), OPTIONAL :: ispli
      REAL    (kind=8) coord(2,npoin),emass(:,:), &
     &                 posgp(2),weigp(2),shape(:,:),deriv(:,:,:),sumat

      END SUBROUTINE lumas3
