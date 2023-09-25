 SUBROUTINE stra02(ndime,gtype,l0,e,dx,l)
 !***********************************************************************
 !
 !**** this routine compute the strain measure for element no. 2
 !     2-3d 2-node truss element
 !***********************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4), INTENT(IN) :: ndime, & !spatial dimension
                                 gtype    !type of strain-stress pair
 REAL (kind=8), INTENT(IN) :: l0,       & !original length
                              dx(ndime)   !present distance vector
 REAL (kind=8), INTENT(OUT) :: l,       & !present length ratio
                               e          !present strain measure

 !     compute present length and strain
 l = SQRT(DOT_PRODUCT(dx,dx))                ! length
 l = l/l0                                    ! lambda

 SELECT CASE (gtype)
 CASE (0)
   e = 0.5d0*(l**2 - 1d0)         !Green-Lagrange strain
 CASE (1,2)
   e = LOG(l)                     !Natural (logarithmic) strain
 CASE (3)
   e = l - 1d0                    !Engng. strain
 END SELECT

 RETURN
 END SUBROUTINE stra02
