 SUBROUTINE resv04e( heade, coora, resid, stab)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF INTERNAL FORCES ELEMENT T-L 3-D SOLID
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 ! dummy arguments
 REAL(kind=8),INTENT (IN) :: coora(:,:), &       !present coordinates
                             stab

 REAL(kind=8),INTENT (IN OUT) :: resid(:,:)   !internal nodal forces
 TYPE (ele04), POINTER  :: heade              !pointer to first element

 ! local variables

 INTEGER (kind=4)  i,n      !different indexes

 REAL (kind=8) strel(7)     !local deviatoric 2nd Piola streses

 REAL (kind=8) bmat(12,6),  & !B matrix
               x(3,nnode),  & !nodal coordinates of the element
               g(3,3),      & !deformation gradient
               rl(3,nnode), &!nodal internal equivalent forces
               p

 TYPE (ele04), POINTER :: e      !pointer to an element data

 !***  Computes nodal internal forces

 e => heade                          !point to first element in the list
 DO                                  !loop over all the elements
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   !*** computes actual configuration
   x = coora(:,e%lnods(1:4))
   g = MATMUL(x,e%cd) !x(1:3) Main element derivatives (deformation gradient)
   !       actual configuration B Matrices
   CALL bmat04(e%cd,bmat,g,nnode)

   strel = 0d0
   DO i=1,nnode
     n = e%lnods(i)
     strel = strel + varn(1:7,n)
   END DO
   ! for post-process only
   p = strel(7)/4d0
   e%stint = strel(1:6)/4d0
   e%stint(1:3) = e%stint(1:3) + p
   ! ---------
   strel(1:6) = strel(1:6)*e%dvol/4d0*stab

   CALL proma3(rl(1,1),bmat(1,1),strel(1),12,6,.TRUE.)  !nodal contributions

   DO i=1,nnode             !for each node in the element
     n = e%lnods(i)     !node
     resid(1:3,n) = resid(1:3,n)+rl(1:3,i) !sum internal force
   END DO

   e => e%next          !point to next node in the list
 END DO

 RETURN

 END SUBROUTINE resv04e
