 SUBROUTINE waterb_ini (numfl,ndime,flparm,headf)
 !Reads follower loads
 USE loa_db,ONLY: foll_seg, flpar
 USE npo_db,ONLY: coora
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: numfl, ndime
   TYPE(flpar):: flparm
   TYPE(foll_seg),POINTER:: headf
   !--- Local variables
   REAL (kind=8), PARAMETER  :: pi=3.141592653589793
   INTEGER(kind=4):: nn
   REAL(kind=8):: vol,x0,y0,xc,h,xl,xr,x(3,3),rr
   REAL (kind=8) l1(3),l2(3),t1(3),t3(3),area2
   TYPE(foll_seg),POINTER:: seg
   LOGICAL :: full,axisy

     axisy = flparm%nodcen == 3     !flag
     y0    = flparm%press            !water surface coordinate
     x0    = flparm%timini           !x-coordinate for volume computation
     vol = 0d0             !Initializes volume
     rr = 1d0              !Initializes radius
     IF( flparm%nodcen > 10 ) THEN
       xr=-1e9
       full = .TRUE.
     ELSE
       xl = 0d0
       xr = 0d0  !initializes coordinates for free surfaces
       full = .FALSE.
     END IF
     seg => headf
     DO nn=1,numfl

       IF (ndime == 3) THEN
         x(1,1:3) = coora(1,seg%lnofl(1:3)) - x0 !distance to reference line
         x(2,1:3) = y0 - coora(2,seg%lnofl(1:3)) !depth at nodes
         !*** evaluate the first two side vectors
         l1 = coora(:,seg%lnofl(1)) - coora(:,seg%lnofl(2))                             !side 1
         l2 = coora(:,seg%lnofl(3)) - coora(:,seg%lnofl(1))                             !side 2
         !*** evaluate the cross product => plane normal
         CALL vecpro(l1,l2,t3)                             !normal * area2  !plano yz
         CALL vecuni(3,t3,area2)                           !computes twice area
         t1 = (/-1d0,0d0,0d0/)
         area2 = area2*DOT_PRODUCT(t3,t1)
       ELSE
         x(1,1:2) = coora(1,seg%lnofl(1:2)) - x0 !distance to reference line
         x(2,1:2) = y0 - coora(2,seg%lnofl(1:2)) !depth at nodes
         IF( axisy ) rr = pi*(x(1,1)+x(1,2))     !2 Pi * r
       END IF

       IF( full )THEN
         IF(ndime == 3) THEN
           xl = -MINVAL(x(2,1:3))   !highest point
           IF( xl > xr ) xr = xl   !update
           vol = vol+((x(1,1)+x(1,2)+x(1,3))/3d0)*area2/2d0
           seg%fload = (x(2,1)+x(2,2)+x(2,3))/3d0         !depth at centroid
         ELSE
           xl = -MINVAL(x(2,1:2))   !highest point
           IF( xl > xr ) xr = xl   !update
           vol = vol+((x(1,1)+x(1,2))/2d0)*(x(2,2)-x(2,1))*rr
           seg%fload = (x(2,1)+x(2,2))/2d0         !depth at centroid
         END IF
       ELSE
         IF(ndime == 3) THEN   !only for triangular element
            !IF( x(2,1) >= 0d0 .AND. x(2,2)>= 0d0 .AND. x(2,3) >= 0d0 )THEN    !all sumerged
            vol = vol+((x(1,1)+x(1,2)+x(1,3))/3d0)*area2/2d0
            seg%fload = (x(2,1)+x(2,2)+x(2,3))/3d0   !depth at centroid
         ELSE
            IF( x(2,1) >= 0d0 .AND. x(2,2) >= 0d0 )THEN    !both sumerged
            vol = vol+((x(1,1)+x(1,2))/2d0)*(x(2,2)-x(2,1))*rr
            seg%fload = (x(2,1)+x(2,2))/2d0         !depth at centroid

           ELSE IF( x(2,1) > 0d0 .OR. x(2,2) > 0d0 )THEN !if one node is below line
            h  = x(2,1)/(x(2,1) - x(2,2))           !normalized position from node 1 (+)
            xc = (1d0-h)*x(1,1) + h*x(1,2)          !x-intersection with free surface
             IF( x(2,2) > 0 )THEN                    !node 1 is above line
              IF( axisy ) rr = pi*(xc+x(1,2))       !2 Pi * r
              vol = vol+(xc    +x(1,2))/2d0*x(2,2)*rr
              seg%fload = x(2,2)*(1d0-h)/2d0        !equivalent pressure at centroid
              IF( xc > xr ) xr = xc
             ELSE                                    !node 2 is above line
              IF( axisy ) rr = pi*(x(1,1)+xc )      !2 Pi * r
              vol = vol-(x(1,1)+xc    )/2d0*x(2,1)*rr
              seg%fload = x(2,1)*h/2d0              !equivalent pressure at centroid(-)
              IF( xc < xl ) xl = xc
             END IF
           ELSE
           seg%fload = 0d0                         !pressure at centroid
           END IF
         END IF
       END IF
       seg => seg%next
     END DO   !loop over SEGMENTS

     flparm%volin = vol
     flparm%vol   = vol
     SELECT CASE (flparm%nodcen)
     CASE  (1)
       flparm%pext  = xr           !symmetric problems
     CASE  (2)
       flparm%pext  = xr - xl      !closed bags
     CASE  (3)
       flparm%pext  = xr*xr*pi*2d0 !axisymmetric problems
     CASE (11:)                    !update initial pressure
       IF( flparm%press == 0d0 ) THEN
         flparm%press = xr
         seg => headf
         DO nn=1,numfl
           seg%fload = seg%fload + xr !depth at centroid
           seg => seg%next
         END DO
       END IF
     END SELECT

 RETURN
 END SUBROUTINE waterb_ini
