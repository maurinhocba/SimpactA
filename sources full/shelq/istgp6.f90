 SUBROUTINE istgp6(ngaus,nst,s,e,d,ambda,matty,b,sf)
 !*****************************************************************************
 !
 !*****evaluates total resultant stresses for shell element
 !     for linear elastic isotropic material
 !
 !****************************************************************************
 IMPLICIT NONE

 !                        routine parameters

 INTEGER (kind=4), INTENT(IN) :: ngaus,nst,matty

 REAL (kind=8), INTENT(IN) :: e(nst,ngaus),ambda(2,ngaus),d(15),b(21),sf(3)
 REAL (kind=8), INTENT(OUT) :: s(nst,ngaus)

 !                        local variables
 INTEGER (kind=4) :: g  !Gauss point Index
 REAL    (kind=8) :: h  !present thickness ratio

 DO g = 1,ngaus
   h = ambda(2,g)
   !h = 1d0
   !                     compute effective stresses
   SELECT CASE (matty)
   CASE (1:2)
     s(1,g) = d(1) * e(1,g) + d(2) * e(2,g)
     s(2,g) = d(2) * e(1,g) + d(1) * e(2,g)
     s(3,g) = d(3) * e(3,g)
     s(4,g) =(d(4) * e(4,g) + d(5) * e(5,g))*h**2
     s(5,g) =(d(5) * e(4,g) + d(4) * e(5,g))*h**2
     s(6,g) = d(6) * e(6,g) * h**2
     s(7,g) = d(7) * e(7,g) * h
     s(8,g) = d(7) * e(8,g) * h
   CASE (3)
     s(1,g) = d(1) * e(1,g) + d(2) * e(2,g)
     s(2,g) = d(2) * e(1,g) + d(3) * e(2,g)
     s(3,g) = d(4) * e(3,g)
     s(4,g) =(d(5) * e(4,g) + d(6) * e(5,g))*h**2
     s(5,g) =(d(6) * e(4,g) + d(7) * e(5,g))*h**2
     s(6,g) = d(8) * e(6,g) * h**2
     s(7,g) = d(9) * e(7,g) * h
     s(8,g) =d(10) * e(8,g) * h
   CASE (4)
     s(1,g) =  d(1)*e(1,g)+ d(2)*e(2,g)+ d(3)*e(3,g)
     s(2,g) =  d(2)*e(1,g)+ d(4)*e(2,g)+ d(5)*e(3,g)
     s(3,g) =  d(3)*e(1,g)+ d(5)*e(2,g)+ d(6)*e(3,g)
     s(4,g) =( d(7)*e(4,g)+ d(8)*e(5,g)+ d(9)*e(6,g))*h**2
     s(5,g) =( d(8)*e(4,g)+d(10)*e(5,g)+d(11)*e(6,g))*h**2
     s(6,g) =( d(9)*e(4,g)+d(11)*e(5,g)+d(12)*e(6,g))*h**2
     s(7,g) =(d(13)*e(7,g)+d(14)*e(8,g))*h
     s(8,g) =(d(14)*e(7,g)+d(15)*e(8,g))*h
   CASE (5)
     s(1,g) =  d(1)*e(1,g)+ d(2)*e(2,g)+ d(3)*e(3,g)             &
            +  b(1)*e(4,g)+ b(2)*e(5,g)+ b(3)*e(6,g)             !&
            !+ b(10)*e(7,g)+b(11)*e(8,g)
     s(2,g) =  d(2)*e(1,g)+ d(4)*e(2,g)+ d(5)*e(3,g)             &
            +  b(4)*e(4,g)+ b(5)*e(5,g)+ b(6)*e(6,g)             !&
            !+ b(12)*e(7,g)+b(13)*e(8,g)
     s(3,g) =  d(3)*e(1,g)+ d(5)*e(2,g)+ d(6)*e(3,g)             &
            +  b(7)*e(4,g)+ b(8)*e(5,g)+ b(9)*e(6,g)             !&
            !+ b(14)*e(7,g)+b(15)*e(8,g)
     s(4,g) =( d(7)*e(4,g)+ d(8)*e(5,g)+ d(9)*e(6,g))*h**2       &
            +  b(1)*e(1,g)+ b(4)*e(2,g)+ b(7)*e(3,g)             !&
            !+ b(16)*e(7,g)+b(17)*e(8,g)
     s(5,g) =( d(8)*e(4,g)+d(10)*e(5,g)+d(11)*e(6,g))*h**2       &
            +  b(2)*e(1,g)+ b(5)*e(2,g)+ b(8)*e(3,g)             !&
            !+ b(18)*e(7,g)+b(19)*e(8,g)
     s(6,g) =( d(9)*e(4,g)+d(11)*e(5,g)+d(12)*e(6,g))*h**2       &
            +  b(3)*e(1,g)+ b(6)*e(2,g)+ b(9)*e(3,g)             !&
            !+ b(20)*e(7,g)+b(21)*e(8,g)
     s(7,g) =(d(13)*sf(1)*e(7,g)+d(14)*sf(2)*e(8,g))*h           !&
            !+ b(11)*e(1,g)+b(13)*e(2,g)+b(15)*e(3,g)             &
            !+ b(16)*e(4,g)+b(18)*e(5,g)+b(20)*e(6,g)
     s(8,g) =(d(14)*sf(2)*e(7,g)+d(15)*sf(3)*e(8,g))*h           !&
            !+ b(10)*e(1,g)+b(12)*e(2,g)+b(14)*e(3,g)             &
            !+ b(17)*e(4,g)+b(19)*e(5,g)+b(21)*e(6,g)
   END SELECT
 END DO
 RETURN
 END SUBROUTINE istgp6
