 SUBROUTINE gaus04n(npoin,coord)

 ! Generate nodal connectivities and space for B matrices

 IMPLICIT NONE
 INTEGER (kind=4) :: npoin
 REAL (kind=8) :: coord(:,:)

 ! local variables
 INTEGER (kind=4) i
 REAL (kind=8), ALLOCATABLE ::  x(:,:)         !nodal coordinates
 REAL (kind=8) :: g(3,3),stran(6)


 !_____________________________


 DO i=1,npoin
   IF( ncc(1,i) == 1 )CYCLE
   ncard(i)%p = ncard(i)%p/varn(0,i)
   ALLOCATE( x(3,ncc(1,i)) )
   x = coord(:,nlnod(i)%p)
   g = MATMUL(x,ncard(i)%p) !Main derivatives (deformation gradient)
   stran(1) = g(1,1)*g(1,1) + g(2,1)*g(2,1) + g(3,1)*g(3,1)  - 1d0
   stran(2) = g(1,2)*g(1,2) + g(2,2)*g(2,2) + g(3,2)*g(3,2)  - 1d0
   stran(3) = g(1,3)*g(1,3) + g(2,3)*g(2,3) + g(3,3)*g(3,3)  - 1d0
   stran(4) = g(1,1)*g(1,2) + g(2,1)*g(2,2) + g(3,1)*g(3,2)
   stran(5) = g(1,1)*g(1,3) + g(2,1)*g(2,3) + g(3,1)*g(3,3)
   stran(6) = g(1,2)*g(1,3) + g(2,2)*g(2,3) + g(3,2)*g(3,3)
   varn(8:13,i) = stran
   DEALLOCATE( x )
 END DO

 RETURN
 END SUBROUTINE gaus04n
