 SUBROUTINE bbra13(bmat,lb,sides,nbn,esi,bb,ln)

 ! assembles standard bending B matrix and branchin side B matrix

 IMPLICIT NONE
 INTEGER(Kind=4), INTENT(IN) :: nbn,      & !number of nodes in the extended patch
                                sides(3), & !side types
                                ln(3)       !position in side conns
 REAL(kind=8), INTENT(IN) :: bmat(:,:,:), & !standard bending B matrix
                             lb             !h/h0 at main triangle
 REAL(kind=8), INTENT(OUT) :: bb(:,:,:)     !assembled B matrix
 TYPE(pside), POINTER :: esi(:)             !array of pointers to branching sides

 ! local variables
 INTEGER(kind=4) :: i,j,n,jj,kk,k,nn
 REAL(kind=8) :: r11,r22,r12        !factors for curvature contributions
 REAL(kind=8) :: bs(3,-1:16)           !average rotation B matrix at branching side
 TYPE(sideb), POINTER :: si         !pointer to the branching side

 !-----------------------++++++++++++++++++++++++++++++++++++++++++

 j = 4                                  !first extra node
 bb (1:3,4:nbn,1:6) = 0d0               !initializes part of extra nodes
 bb (1:3,1:3,1:6)   = bmat(1:3,1:3,1:6) !main triangle contributions
 DO i=1,3                     !for each side
   SELECT CASE ( sides(i) )   !according to side type
   CASE (1)                   !adyacent element exist
     bb (1:3,j,1:6) = bmat(1:3,i+3,1:6)   !adjacent triangle contribution
     j = j+1                 !update node position
   CASE (0)                  !boundary side
     CYCLE                   !nothing
   CASE (-1)                 !branching side
     si => esi(i)%p          !change pointer
     nn = si%nn              !number of nodes in the branching side
     bs(:,-1:nn) = si%bb/lb  !B/lambda_M
     n = -ln(i)              !position in side connectivites
     r11 = bmat(1,i+3,4)
     r22 = bmat(1,i+3,5)
     r12 = bmat(1,i+3,6)
     IF( si%lnods(n) > 0 )THEN  !same order
       jj = hh(3,i)  !first node
       kk = hh(2,i)  !second node of the side
     ELSE                       !reverse
       jj = hh(2,i)  !first node
       kk = hh(3,i)  !second node of the side
       r11 = -r11    ! change sign to coefficients
       r22 = -r22
       r12 = -r12
     END IF
     bb (1:3,jj,4) = bb(1:3,jj,4) + r11*bs(1:3,-1)   !first node of the side
     bb (1:3,jj,5) = bb(1:3,jj,5) + r22*bs(1:3,-1)   !
     bb (1:3,jj,6) = bb(1:3,jj,6) + r12*bs(1:3,-1)   !
     bb (1:3,kk,4) = bb(1:3,kk,4) + r11*bs(1:3, 0)   !second node of the side
     bb (1:3,kk,5) = bb(1:3,kk,5) + r22*bs(1:3, 0)   !
     bb (1:3,kk,6) = bb(1:3,kk,6) + r12*bs(1:3, 0)   !

     DO k=1,nn
       IF( k == n )THEN
         bb (1:3,i,4) = bb(1:3,i,4) + r11*bs(1:3,k)   !main triangle node
         bb (1:3,i,5) = bb(1:3,i,5) + r22*bs(1:3,k)   !
         bb (1:3,i,6) = bb(1:3,i,6) + r12*bs(1:3,k)   !
       ELSE
         bb (1:3,j,4) = bb(1:3,j,4) + r11*bs(1:3,k)   !other node
         bb (1:3,j,5) = bb(1:3,j,5) + r22*bs(1:3,k)   !
         bb (1:3,j,6) = bb(1:3,j,6) + r12*bs(1:3,k)   !
         j = j+1
       END IF
     END DO
   END SELECT
 END DO
 RETURN
 END SUBROUTINE bbra13
