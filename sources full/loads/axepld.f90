      SUBROUTINE axepld(coorn,trans,x,y,area2,flag)
!***********************************************************************
!
!****this routine compute the element local axes system
!****for the 3 node element
!****(local x-axis is aligned with a side 1-2
!***********************************************************************
      IMPLICIT NONE

      REAL (kind=8),INTENT(IN) :: coorn(3,3)
      REAL (kind=8),INTENT(OUT) :: trans(3,3),x(3),y(3),area2
      LOGICAL, INTENT(IN) :: flag

!
      INTEGER (kind=4) i,j
      REAL (kind=8) xno1,xno3,coorb(2,3)

      INTERFACE
        INCLUDE 'vecpro.h'
        INCLUDE 'vecuni.h'
      END INTERFACE
!*** evaluate the local axes

      trans(1:3,1) = coorn(1:3,2) - coorn(1:3,1)    !parallel to side 1
      trans(1:3,2) = coorn(1:3,3) - coorn(1:3,1)    !parallel to side 3

!*** evaluate the cross product  T3 = T1 x T2

!      CALL vecpro(trans(1:3,1),trans(1:3,2),trans(1:3,3))
      trans(1,3) = trans(2,1)*trans(3,2) - trans(3,1)*trans(2,2)
      trans(2,3) = trans(3,1)*trans(1,2) - trans(1,1)*trans(3,2)
      trans(3,3) = trans(1,1)*trans(2,2) - trans(2,1)*trans(1,2)

!      CALL vecuni(3,trans(1:3,1),xno1)
      xno1 = SQRT(DOT_PRODUCT(trans(:,1),trans(:,1)))
      trans(:,1) = trans(:,1)/xno1
!      CALL vecuni(3,trans(1:3,3),xno3)
      xno3 = SQRT(DOT_PRODUCT(trans(:,3),trans(:,3)))
      trans(:,3) = trans(:,3)/xno3

      IF (xno3 == 0.0D0) THEN
         WRITE (*,*) 'ERROR: AREA2 == 0 '
         write(*,"(3e15.5)")coorn
         WRITE (55,*,ERR=9999) 'ERROR: AREA2 == 0 '
         write(55,"(3e15.5)",ERR=9999) coorn
         CALL runen3('AXEPLD: ELEMENT DISTORTION OCCURRED')
      END IF

!      CALL vecpro(trans(1:3,3),trans(1:3,1),trans(1:3,2))
      trans(1,2) = trans(2,3)*trans(3,1) - trans(3,3)*trans(2,1)
      trans(2,2) = trans(3,3)*trans(1,1) - trans(1,3)*trans(3,1)
      trans(3,2) = trans(1,3)*trans(2,1) - trans(2,3)*trans(1,1)

      area2 = xno3

      IF(flag)THEN

       !*** find the local coordinates

!        coorb = MATMUL(TRANSPOSE(trans(1:3,1:2)),coorn)
        DO i=1,2
          DO j=1,3
            coorb(i,j) = DOT_PRODUCT(trans(:,i),coorn(:,j))
          END DO
        END DO

        !side proyections
        x(1) = coorb(1,1) - coorb(1,2)
        y(1) = coorb(2,1) - coorb(2,2)
        x(2) = coorb(1,2) - coorb(1,3)
        y(2) = coorb(2,2) - coorb(2,3)
        x(3) = coorb(1,3) - coorb(1,1)      ! -x(1) -x(2)
        y(3) = coorb(2,3) - coorb(2,1)      ! -y(1) -y(2)

      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE axepld
