      SUBROUTINE lumas5(ndofn,npoin,nelem,ngaus,nnode,lnods,matno,coord, &
     &                  emass,weigp,shape,deriv,iwrit,sumat)
!******************************************************************
!
! *** calculates lumped mass for 3d 4- or 8-node element
!
!******************************************************************
      IMPLICIT NONE
      INTEGER (kind=4) ndofn,nelem,ngaus,nnode,npoin,iwrit, &
     &                 lnods(:,:),matno(:)
      REAL    (kind=8) coord(3,npoin),emass(ndofn,*), &
     &                 weigp(:),shape(:,:),deriv(:,:,:),sumat

      END SUBROUTINE lumas5
