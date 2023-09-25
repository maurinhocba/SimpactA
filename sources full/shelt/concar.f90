SUBROUTINE concar(con,car,j,flag)
!*********************************************************
!     transform a tensor from convected to cartesian & viceversa
!     covariant (strain, note that the third element is twice the tensor value)
!         FLAG = 1  from convected to cartesian
!         FLAG = 2  from cartesian to convected
!     contravariant (stress)
!         FLAG = 3  from convected to cartesian
!         FLAG = 4  from cartesian to convected
!     DATA
!         CON in convected system
!         CAR in cartesian system
!         J   inverse of the jacobian matrix
!*********************************************************
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: flag
REAL (kind=8),INTENT(IN) :: j(2,2)
REAL (kind=8),INTENT(IN OUT) :: con(3),car(3)

REAL (kind=8) :: value

SELECT CASE (flag)
CASE (1)
  car(1) =  con(1)*j(1,1)*j(1,1)+                                 &
            con(2)*j(1,2)*j(1,2)+                                 &
            con(3)*j(1,1)*j(1,2)
  car(2) =  con(1)*j(2,1)*j(2,1)+                                 &
            con(2)*j(2,2)*j(2,2)+                                 &
            con(3)*j(2,1)*j(2,2)
  car(3)=2*(con(1)*j(2,1)*j(1,1)+                                 &
            con(2)*j(2,2)*j(1,2))+                                &
            con(3)*(j(1,1)*j(2,2)+j(1,2)*j(2,1))
CASE (2)
  value  =  1d0/(j(1,1)*j(2,2)-j(1,2)*j(2,1))**2
  con(1) = (car(1)*j(2,2)*j(2,2)+                                 &
            car(2)*j(1,2)*j(1,2)-                                 &
            car(3)*j(2,2)*j(1,2))*value
  con(2) = (car(1)*j(2,1)*j(2,1)+                                 &
            car(2)*j(1,1)*j(1,1)-                                 &
            car(3)*j(2,1)*j(1,1))*value
  con(3) = (-2*(car(1)*j(2,1)*j(2,2)+                             &
                car(2)*j(1,1)*j(1,2))+                            &
                car(3)*(j(1,1)*j(2,2)+j(1,2)*j(2,1)))*value
CASE (3)
  value  =  1d0/(j(1,1)*j(2,2)-j(1,2)*j(2,1))**2
  car(1) = (con(1)*j(2,2)*j(2,2)+                                 &
            con(2)*j(2,1)*j(2,1)-                                 &
            con(3)*j(2,2)*j(2,1)*2d0)*value
  car(2) = (con(1)*j(1,2)*j(1,2)+                                 &
            con(2)*j(1,1)*j(1,1)-                                 &
            con(3)*j(1,2)*j(1,1)*2d0)*value
  car(3) =(-con(1)*j(2,2)*j(1,2)-                                 &
            con(2)*j(2,1)*j(1,1)+                                 &
            con(3)*(j(1,1)*j(2,2)+j(1,2)*j(2,1)))*value
CASE (4)
  con(1) =  car(1)*j(1,1)*j(1,1)+                                 &
            car(2)*j(2,1)*j(2,1)+                                 &
            car(3)*j(1,1)*j(2,1)*2d0
  con(2) =  car(1)*j(1,2)*j(1,2)+                                 &
            car(2)*j(2,2)*j(2,2)+                                 &
            car(3)*j(1,2)*j(2,2)*2d0
  con(3) =  car(1)*j(1,1)*j(1,2)+                                 &
            car(2)*j(2,1)*j(2,2)+                                 &
            car(3)*(j(1,1)*j(2,2)+j(1,2)*j(2,1))
END SELECT
RETURN

END SUBROUTINE concar
