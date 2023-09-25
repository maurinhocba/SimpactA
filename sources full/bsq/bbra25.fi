 SUBROUTINE bbra25(bmat,bb,lb,sides,nbn,esi,bbr,bbrk,lm,dt,ns)

 ! assembles standard bending B matrix and branchin side B matrix

 IMPLICIT NONE
 INTEGER(Kind=4), INTENT(IN) :: nbn,      & !number of nodes in the extended patch
                                sides(4), & !side types
                                lm(8)       !position in side conns
 REAL(kind=8), INTENT(IN) :: bmat(:,:,:), & !standard bending B matrix
                             bb(:,:,:),   & !standard stabilization bending matrix
                             lb,          & !h/h0 at main triangle
                             dt(4),       &
                             ns(2,4)        !Element side proyections
 REAL(kind=8), INTENT(OUT) :: bbr(:,:,:), & !assembled B matrix
                              bbrk(:,:,:)
 TYPE(pside), POINTER :: esi(:)             !array of pointers to branching sides

 ! local variables
 INTEGER(kind=4) :: i,j,l,ll,ii,n,m,jj,kk,k,nn,mm,mn
 REAL(kind=8) :: r11,r22,r12,n1,n2,nl,c1,c2,c3,c4   !factors for curvature contributions
 REAL(kind=8) :: bs(3,-1:18)           !average rotation B matrix at branching side
 REAL(kind=8) :: bm(3,20,3)
 REAL(kind=8) :: bmg(3,20,4)
 TYPE(sideb), POINTER :: si         !pointer to the branching side

 !-----------------------++++++++++++++++++++++++++++++++++++++++++

 m = 5                                   !first extra node
 bm = 0d0
 bmg = 0d0
 bbr (1:3,5:nbn,1:6) = 0d0               !initializes part of extra nodes
 bbr (1:3,1:4,1:6)   = bmat(1:3,1:4,1:6) !main quadrilateral contributions
 bbrk(1:3,5:nbn,1:4) = 0d0               !initializes part of extra nodes
 bbrk(1:3,1:4,1:4) = bb(1:3,1:4,1:4)     !main quadrilateral contributions

 DO i=1,4                     !for each side
   SELECT CASE ( sides(i) )   !according to side type
   CASE (1)                   !adyacent element exist
     DO j=3,4
       k = ln(j,i)
       bbr (1:3,m,1:6) = bmat(1:3,k,1:6)    !adjacent quadrilateral contribution
       bbrk(1:3,m,i) = bb(1:3,k,i)          !adjacent quadrilateral contribution
       m = m+1
     END DO
   CASE (0)                  !boundary side
     CYCLE                   !nothing
   CASE (-1)                 !branching side
     ! I - P A R T:  C O N T R I B U T I O N  O F  M A I N  E L E M E N T
     si => esi(i)%p          !change pointer
     nn = si%nn              !2*nn number of nodes in the branching side
     bs(:,-1:2*nn) = si%bbr/lb  !B/lambda_M
     n = -lm(ln(3,i)-4)              !position in side connectivites
     k = ln(3,i)
     r11 = bmat(1,k,7)
     r22 = bmat(2,k,7)
     r12 = bmat(3,k,7)
     c4  = +1d0

     IF( si%lnods(2*n) > 0 )THEN  !same order
       jj = hh(3,i)  !first node
       kk = hh(2,i)  !second node of the side
       ii = nextn(i+1)
       ll = nextn(i+2)
       mm = m
       mn = m+1
     ELSE                       !reverse
       jj = hh(2,i)  !first node
       kk = hh(3,i)  !second node of the side
       ii = nextn(i+1)
       ll = nextn(i+2)
       mm = m
       mn = m+1
       r11 = -r11    ! change sign to coefficients
       r22 = -r22
       r12 = -r12
       c4  = -c4
     END IF

     ! II - P A R T:
       !     normal to the side rotation
       bm(1:3,jj,1) = bm(1:3,jj,1) + r11*bs(1:3,-1)   !first node of the side
       bm(1:3,jj,2) = bm(1:3,jj,2) + r22*bs(1:3,-1)   !
       bm(1:3,jj,3) = bm(1:3,jj,3) + r12*bs(1:3,-1)   !

       bmg(1:3,jj,i) = bmg(1:3,jj,i) + c4*bs(1:3,-1)   !first node of the side

       !bbrk(1:3,jj,i) = bbrk(1:3,jj,i) - (c1*bm(1:3,jj,1)+c2*bm(1:3,jj,2)+c3*bm(1:3,jj,3)) + &
       !                 bmg(1:3,jj)/dt(i)

       bm(1:3,kk,1) = bm(1:3,kk,1) + r11*bs(1:3, 0)   !second node of the side
       bm(1:3,kk,2) = bm(1:3,kk,2) + r22*bs(1:3, 0)   !
       bm(1:3,kk,3) = bm(1:3,kk,3) + r12*bs(1:3, 0)   !

       bmg(1:3,kk,i) = bmg(1:3,kk,i) + c4*bs(1:3, 0) !second node of the side

       !bbrk(1:3,kk,i) = bbrk(1:3,kk,i) - (c1*bm(1:3,kk,1)+c2*bm(1:3,kk,2)+c3*bm(1:3,kk,3)) + &
       !                 bmg(1:3,kk)/dt(i)

       DO k=1,nn
         j = 2*k-1
         IF( k == n )THEN
           bm (1:3,ii,1) = r11*bs(1:3,j)   !main triangle node
           bm (1:3,ii,2) = r22*bs(1:3,j)   !
           bm (1:3,ii,3) = r12*bs(1:3,j)   !

           bmg (1:3,ii,i) = c4*bs(1:3,j)  !main triangle node

           !bbrk(1:3,ii,i) = -(c1*bm(1:3,ii,1)+c2*bm(1:3,ii,2)+c3*bm(1:3,ii,3)) + &
           !                 bmg(1:3,ii)/dt(i)
           !
           bm (1:3,ll,1) = r11*bs(1:3,j+1)   !
           bm (1:3,ll,2) = r22*bs(1:3,j+1)   !
           bm (1:3,ll,3) = r12*bs(1:3,j+1)   !

           bmg (1:3,ll,i) = c4*bs(1:3,j+1)

           !bbrk(1:3,ll,i) = -(c1*bm(1:3,ll,1)+c2*bm(1:3,ll,2)+c3*bm(1:3,ll,3)) + &
           !                 bmg(1:3,ll)/dt(i)

         ELSE
           bm (1:3,mm,1) = r11*bs(1:3,j)   !other node
           bm (1:3,mm,2) = r22*bs(1:3,j)   !
           bm (1:3,mm,3) = r12*bs(1:3,j)   !
           m = m+2                                    !

           bmg (1:3,mm,i) = c4*bs(1:3,j) !other node

           !bbrk(1:3,mm,i) = -(c1*bm(1:3,mm,1)+c2*bm(1:3,mm,2)+c3*bm(1:3,mm,3)) + &
           !                 bmg(1:3,mm)/dt(i)

           bm (1:3,mn,1) = r11*bs(1:3,j+1)   !
           bm (1:3,mn,2) = r22*bs(1:3,j+1)   !
           bm (1:3,mn,3) = r12*bs(1:3,j+1)   !

           bmg (1:3,mn,i) = c4*bs(1:3,j+1)

           !bbrk(1:3,mn,i) = -(c1*bm(1:3,mn,1)+c2*bm(1:3,mn,2)+c3*bm(1:3,mn,3)) + &
           !                 bmg(1:3,mn)/dt(i)

           mm = mm+2
           mn = mn+2

         END IF
       END DO
   END SELECT
 END DO

 bbr(1:3,1:nbn,4:6) = bbr(1:3,1:nbn,4:6) + bm(1:3,1:nbn,1:3)
 bbrk(1:3,1:nbn,1:4) = bbrk(1:3,1:nbn,1:4) + bmg(1:3,1:nbn,1:4)

 DO l=1,4
   nl = SQRT(ns(1,l)**2+ns(2,l)**2)      !side length
   n1 = ns(1,l)/nl                       !normalize normal vector
   n2 = ns(2,l)/nl
   !     normal to the side rotation
   c1 = n1**2     ! n1^2
   c2 = n2**2     ! n2^2
   c3 = n1*n2     ! n1 n2
   bbrk(:,1:nbn,l) = -bbrk(:,1:nbn,l)/dt(l) - (c1*bbr(:,1:4,4) + c2*bbr(:,1:4,5) + c3*bbr(:,1:4,6))
 END DO

 RETURN
 END SUBROUTINE bbra25
