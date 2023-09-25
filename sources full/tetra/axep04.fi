 SUBROUTINE axep04(x,cd,dvol,lb)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 4 node tetrahedra, and for the adjacent elements
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN) :: x(:,:), &  !nodal coordinates
                             lb(3,3)    !rotation between standard X1 and local X1
 REAL (kind=8),INTENT(OUT) :: cd(4,3),dvol

 ! local variables
 INTEGER (kind=4) j

 REAL (kind=8) :: cda(4,3)

 !*** compute jacobian at the main triangle

 CALL jacob5(cda(1,1),nd0(1,1),dvol,x(1,1),4,0,j)

 IF (dvol <= 0.0D0 .OR. j == 1)THEN
   WRITE (55,*) 'ERROR: VOLUME <= 0 '
   WRITE(55,"(3e15.5)")x(1:3,1:4)
   CALL runend('AXEP04: WRONG ELEMENT DEFINITION   ')
 END IF

 ! change to local system
 cd(:,1) = cda(:,1)*lb(1,1) + cda(:,2)*lb(2,1) + cda(:,3)*lb(3,1)
 cd(:,2) = cda(:,1)*lb(1,2) + cda(:,2)*lb(2,2) + cda(:,3)*lb(3,2)
 cd(:,3) = cda(:,1)*lb(1,3) + cda(:,2)*lb(2,3) + cda(:,3)*lb(3,3)

 dvol = dvol/6d0                 !element volume

 RETURN
 END SUBROUTINE axep04
