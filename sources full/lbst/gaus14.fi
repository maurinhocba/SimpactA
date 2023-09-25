 SUBROUTINE gaus14(heade,coord,coora,iffix,gauss,angdf,nonrg,locax,plstr)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 14
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: iffix(:)
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              coora(:,:), & !present coordinates
                              angdf         !default angle X_1 - Ort_1
 TYPE (ele14), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists
 LOGICAL, INTENT (IN) :: nonrg        !TRUE if non-regular meshes will be considered
 INTEGER(kind=4), INTENT (IN) :: locax      ! local axis
 INTEGER, INTENT (IN OUT) :: plstr    !TRUE if plastic strains are to be computed

 INTEGER (kind=4) i,j,k,m,n,l
 REAL (kind=8) x(3,6),stra0(6),t(3,2,0:3),t3(3),h(3,3)
 REAL (kind=8) angle,factor,ls
 TYPE (ele14), POINTER :: e
 LOGICAL :: sides(3),is(3)

 !     ---------------------------------------------------------
 factor = ASIN(1d0)/90d0         !pi/180
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   ! is = .FALSE.                    !initializes boundary codes
   x(:,1:3) = coord(:,e%lnods(1:3))     !get original coordinates
   DO m =1,3                       !for all nodes in the patch
     i = m+3                       !associated side element
     n = e%lnods(i)                !node number
     IF(n /= 0) THEN               !if node exist
       x(:,i) = coord(:,n)         !get original coordinates
       sides(m) = .TRUE.           !side element exists
     ELSE
       sides(m) = .FALSE.          !side element does not exists
       j = e%lnods(hh(2,m))        !first node of the side
       k = e%lnods(hh(1,m))        !second node of the side
       IF( iffix(j) < 0 .AND. iffix(k) < 0 )THEN !node on a side
         l = -iffix(j)  !compute position of the i+1 node
         t(:,1,1) =  coord(:,nndpd(2,l))  !node 1
         t(:,2,1) =  coord(:,nndpd(3,l))  !node 2
         l = -iffix(k)  !compute position of the i+2 node
         t(:,1,2) =  coord(:,nndpd(3,l))  !node 3
         t(:,2,2) =  coord(:,nndpd(2,l))  !node 4 of the quadrilateral
         t(:,1,3) = -t(:,1,1)+t(:,2,1)+t(:,1,2)-t(:,2,2) !phi,xita at the center
         t(:,2,3) = x(:,hh(1,m))-x(:,hh(2,m))    !side vector
         CALL vecpro(t(1,1,3),t(1,2,3),x(1,i)) !outward (to the solid) normal vector
         CALL vecuni(3,x(1,i),ls)          !unit normal to the quad
         e%a(1:3,m) = -x(:,i)           !unit normal to the solid side
         e%b(1,m) = 1d0                 !flag
       ELSE
         !x(1:3,i) = 0d0         !null values (unnecessary perhaps)
         !e%a(:,m) = 0d0         !null values (unnecessary perhaps)
         IF(.NOT.gauss)e%b(:,m) = 0d0         !null values (1,) is a flag
       END IF
       is(m) = ( iffix(j) /= 0 .AND. iffix(k) /= 0 )  !clamped side ?
     END IF
   END DO

   IF(.NOT.gauss)THEN    !Compute initial constants if requested
     ! compute element local systems, side proyections and element areas
     angle = ( angdf + e%angle )*factor
     CALL axep14(x,e%a,e%b,e%area1,sides,is,angle,e%ci,nonrg,locax)

     ! compute second fundamental form of the undeformed structure
     CALL stra14(e%a,e%b,x,t3,t,h,sides,is,e%ci,nonrg,stra0,e%lb)
     e%stra0 = stra0(4:6)      !keep original curvatures

   END IF
   IF( plstr == 1 )  CALL stst14( e, iffix, coora, .FALSE., nonrg)
   e => e%next          !point to next element
 END DO
 gauss = .TRUE.
 plstr = 0

 RETURN
 END SUBROUTINE gaus14
