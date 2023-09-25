 SUBROUTINE bmem11(ntype,l0,r0,lb,t,b,sh)
   ! compute membrame B matrix for 2-d beam/shell rotation free elements
   IMPLICIT NONE
   INTEGER (kind=4), INTENT(IN) :: ntype  !problem type
   REAL (kind=8), INTENT(IN) :: l0, &     !initial length
                                r0, &     !initial radius
                                lb(3), &  !stretchings
                                t(2),  &  !tangent vector
                                sh(2)     !shape functions
   REAL (kind=8), INTENT(OUT) :: b(2,2,2) !Bm
   ! local variable
   REAL (kind=8) :: fac

   ! axial strain
   fac = lb(1)/l0              !Ls/L0
   b(1:2,1,1) = -fac*t         !first node
   b(1:2,2,1) =  fac*t         !second node
   ! circunferencial strain
   IF( ntype == 2 )THEN
     b(:,:,2) = 0d0            !zero
   ELSE IF( ntype == 3 )THEN   !axilsymmetryc problems
     fac = lb(3)/r0            !lt/r0
     b(1,1:2,2) =  fac*sh
     b(2,1:2,2) =  0d0
   END IF

   RETURN
 END SUBROUTINE bmem11
