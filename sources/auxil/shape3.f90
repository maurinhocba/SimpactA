      SUBROUTINE shape3(deriv,shape,s,t,nnode )
      !********************************************************************
      !
      !***  calculates shape functions and their derivatives for 2d elements
      !
      !********************************************************************
      IMPLICIT NONE
      INTEGER (kind=4) :: nnode
      REAL (kind=8)    deriv(nnode,2),shape(nnode),s,t
      REAL (kind=8)    st

      st = s*t

      IF (nnode == 4) THEN

         !  shape functions for 4 noded element

         shape(1) = (1D0 - s - t + st)/4D0
         shape(2) = (1D0 + s - t - st)/4D0
         shape(3) = (1D0 + s + t + st)/4D0
         shape(4) = (1D0 - s + t - st)/4D0

         ! and derivatives

         deriv(1,1) = (-1d0 + t)/4d0
         deriv(2,1) = -deriv(1,1)
         deriv(3,1) = ( 1d0 + t)/4d0
         deriv(4,1) = -deriv(3,1)
         deriv(1,2) = (-1d0 + s)/4d0
         deriv(2,2) = (-1d0 - s)/4d0
         deriv(3,2) = -deriv(2,2)
         deriv(4,2) = -deriv(1,2)

      ELSE IF (nnode == 3) THEN

         !  shape functions for 3-node element

         shape(1) = 1.d0 - s - t
         shape(2) = s
         shape(3) = t

         !  and derivatives

         deriv(1,1) = -1.d0
         deriv(1,2) = -1.d0
         deriv(2,1) =  1.d0
         deriv(2,2) =  0.d0
         deriv(3,1) =  0.d0
         deriv(3,2) =  1.d0

      END IF

      END SUBROUTINE shape3
