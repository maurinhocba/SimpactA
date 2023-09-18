 SUBROUTINE fixdpd(x,euler)
 !***********************************************************************
 !
 !     updates configuration and factors of dependant displacements
 !
 !***********************************************************************
 USE ctrl_db, ONLY : ndime,neulr
 USE kinc_db, ONLY : ndepd,npsdf,nndpd,ftsdf,distd
 USE npo_db,  ONLY : naeul,ifpre
 IMPLICIT NONE
 REAL (kind=8),INTENT(IN OUT) :: x(:,:)
 REAL (kind=8), POINTER :: euler(:,:)

 INTEGER (kind=4) :: j,k,n,n1,n2
 REAL    (kind=8) :: d(ndime),t1,t2
 LOGICAL :: eul

 eul = neulr > 0  !consider local systems

 DO n=1,ndepd
   n1 = nndpd(1,n)     !slave node
   n2 = nndpd(2,n)     !master node
   j  = -ifpre(1,n1)   !position in the pointer array
   k  = npsdf(j)       !pointer to ftsdf
   d = distd(1:ndime,n)!default
   IF( eul )THEN
     IF( naeul(n2) ) THEN
       IF(ndime == 2) THEN            !for 2-D problems
         t1 = COS(euler(1,n2))        !Cosine of angle at master node
         t2 = SIN(euler(1,n2))        !sine of angle at master node
         d(1) = t1*distd(1,n) - t2*distd(2,n) !actual distance vector (global coordinates)
         d(2) = t2*distd(1,n) + t1*distd(2,n)
         ftsdf(k+1) = -t2*distd(1,n) - t1*distd(2,n)  !u_1  theta dependency
         ftsdf(k+3) =  t1*distd(1,n) - t2*distd(2,n)  !u_2  theta dependency
         IF(naeul(n1)) euler(1,n1) = euler(1,n2) + distd(3,n) !update angle
       ELSE
         d = MATMUL(RESHAPE(euler(1:9,n2),(/3,3/)),distd(1:3,n))
         !ftsdf(k   ) =  1d0                                            !u_1 u_1 does not change
         ftsdf(k+ 1) = -euler(4,n2)*distd(3,n) + euler(7,n2)*distd(2,n) !u_1 Theta dependency
         ftsdf(k+ 2) = -euler(7,n2)*distd(1,n) + euler(1,n2)*distd(3,n)
         ftsdf(k+ 3) = -euler(1,n2)*distd(2,n) + euler(4,n2)*distd(1,n)
         !ftsdf(k+4 ) =  1d0                                            !u_2 u_2 does not change
         ftsdf(k+ 5) = -euler(5,n2)*distd(3,n) + euler(8,n2)*distd(2,n) !u_2 Theta dependency
         ftsdf(k+ 6) = -euler(8,n2)*distd(1,n) + euler(2,n2)*distd(3,n)
         ftsdf(k+ 7) = -euler(2,n2)*distd(2,n) + euler(5,n2)*distd(1,n)
         !ftsdf(k+8 ) =  1d0                                            !u_3 u_3 does not change
         ftsdf(k+ 9) = -euler(6,n2)*distd(3,n) + euler(9,n2)*distd(2,n) !u_3 Theta dependency
         ftsdf(k+10) = -euler(9,n2)*distd(1,n) + euler(3,n2)*distd(3,n)
         ftsdf(k+11) = -euler(3,n2)*distd(2,n) + euler(6,n2)*distd(1,n)
         IF(naeul(n1)) CALL proma1(euler(1,n1),euler(1,n2),ftsdf(k+12),3,3,3) !update local systems at slave node
       END IF
     END IF
   END IF
   x(1:ndime,n1) = d + x(1:ndime,n2)        !update coordinates
 END DO
 RETURN

 END SUBROUTINE fixdpd
