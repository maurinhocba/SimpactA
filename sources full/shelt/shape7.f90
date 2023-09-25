 SUBROUTINE shape7(a2,a3,shape,deriv)
 IMPLICIT NONE
 REAL (kind=8) a2,a3,shape(6),deriv(6,2)
 REAL (kind=8) a1

 a1 = 1d0-a2-a3

 shape(1) = (2d0*a1-1d0)*a1
 shape(2) = (2d0*a2-1d0)*a2
 shape(3) = (2d0*a3-1d0)*a3
 shape(4) = 4d0*a1*a2
 shape(5) = 4d0*a2*a3
 shape(6) = 4d0*a1*a3

 deriv(1,1) =  1d0-4d0*a1
 deriv(2,1) =  4d0*a2-1d0
 deriv(3,1) =  0d0
 deriv(4,1) =  4d0*(a1-a2)
 deriv(5,1) =  4d0*a3
 deriv(6,1) = -4d0*a3
 deriv(1,2) =  1d0-4d0*a1
 deriv(2,2) =  0d0
 deriv(3,2) =  4d0*a3-1d0
 deriv(4,2) = -4d0*a2
 deriv(5,2) =  4d0*a2
 deriv(6,2) =  4d0*(a1-a3)
 RETURN
 END SUBROUTINE shape7
