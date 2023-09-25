 SUBROUTINE projtr(imn,x,lcseg,xita,eta,y)
 !     project over a triangle
 IMPLICIT NONE

 !     arguments
 INTEGER (kind=4), INTENT (IN) :: imn,lcseg(:)
 REAL (kind=8), INTENT (IN) :: x(:,:),y(:)
 REAL (kind=8), INTENT (OUT) :: xita,eta
 !     local variables
 INTEGER (kind=4) inseg,jnseg,knseg,jmn,kmn
 REAL    (kind=8) vj(3),vk(3),auxi(3),vn(3),dot1,dot2,area2


 !.... identify local segment nodes
 inseg = 1
 DO
   IF (lcseg(inseg) == imn) EXIT
   inseg = inseg+1
 END DO
 jnseg = MOD(inseg,3) + 1
 knseg = MOD(inseg + 1, 3) + 1
 !.... identify global master nodes
 jmn = lcseg(jnseg)
 kmn = lcseg(knseg)
 !.... form tangent vectors to master segment
 vj = x(1:3,jmn) - x(1:3,imn)
 vk = x(1:3,kmn) - x(1:3,imn)
 !.... compute normal vector and twice the area of the triangle
 vn(1) = vj(2)*vk(3) - vj(3)*vk(2)
 vn(2) = vj(3)*vk(1) - vj(1)*vk(3)
 vn(3) = vj(1)*vk(2) - vj(2)*vk(1)
 area2 = SQRT(DOT_PRODUCT(vn(:),vn(:)))
 vn = vn/area2
 auxi(1) = vj(2)*y(3) - vj(3)*y(2)
 auxi(2) = vj(3)*y(1) - vj(1)*y(3)
 auxi(3) = vj(1)*y(2) - vj(2)*y(1)
 dot1 = MAX(DOT_PRODUCT(auxi,vn)/area2 , 0d0)
 auxi(1) = y(2)*vk(3) - y(3)*vk(2)
 auxi(2) = y(3)*vk(1) - y(1)*vk(3)
 auxi(3) = y(1)*vk(2) - y(2)*vk(1)
 dot2 = MAX(DOT_PRODUCT(vn,auxi)/area2 , 0d0)
 SELECT CASE (inseg)
 CASE (1)
   xita = dot2
   eta  = dot1
 CASE (2)
   xita = MAX(1d0 - (dot1+dot2), 0d0)
   eta  = dot2
 CASE (3)
   xita = dot1
   eta  = MAX(1d0 - (dot1+dot2), 0d0)
 END SELECT
 RETURN
 END SUBROUTINE projtr
