 SUBROUTINE smog19(nvar,Area1,strsg,strsn)
 !=======================================================================
 ! This routine computes the elemental contribution to nodal values
 !=======================================================================
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4),INTENT(IN):: nvar
   REAL(kind=8),INTENT(IN):: Area1, strsg(nvar)
   REAL(kind=8),INTENT(OUT):: strsn(nvar,3)
   !--- Local variables
   INTEGER(kind=4):: i, j
   REAL(kind=8):: auxva, emass(3,3)

   !Evaluates the element smoothing matrix and inverts it (emass)
   CALL smom19(Area1,emass)

   !Integrate nodal weighted stresses and compute nodal contributions
   strsn(1:nvar,1:3) = 0d0
   DO j=1,3
     auxva = Area1/3d0
     DO i=1,nvar
       strsn(i,j) = strsn(i,j) + auxva*strsg(i)
     END DO
   END DO

   strsn(1:nvar,1:3) = MATMUL(strsn(1:nvar,1:3),TRANSPOSE(emass(1:3,1:3)))

 RETURN
 END SUBROUTINE smog19
