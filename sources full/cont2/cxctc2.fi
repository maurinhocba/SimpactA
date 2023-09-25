 SUBROUTINE cxctc2(x,surf)
 !
 ! Compute nodal and segment normals for self contact
 !
 IMPLICIT NONE
 REAL (KIND=8), INTENT(IN) :: x(:,:)
 TYPE (surf2_db) :: surf
 ! INTENT(OUT) :: surf%tn(2,ncnod) :nodal normals
 !                surf%tc(2,nsegm) :segment normals for self contact
 !                surf%xc(2,nsegm) :segment position of center
 !local variables
 INTEGER (KIND=4) :: ncnod,nsegm,i,j,n,n1,n2
 REAL (KIND=8) :: t(2),long

 !INTEGER (kind=4), POINTER :: nr(:)
 REAL (KIND=8), POINTER :: tn(:,:),tc(:,:),xc(:,:)


 ncnod = surf%ncnod            !number of nodes in the surface
 nsegm = surf%nsegm            !number of segments defining surface
 !  point to surface arrays
 xc => surf%xc                 !(2,nsegm)
 tc => surf%tc                 !(2,nsegm)
 tn => surf%tn                 !(2,ncnod)
! nr => surf%nr

 tn = 0d0                      !averaged normals initialization

 ! Compute segment normal and center
 j = 1                      !initializes for an efficient search
 DO n=1,nsegm               !for each segment
   n1 = surf%lcseg(1,n)     !first node
   n2 = surf%lcseg(2,n)     !second node
   xc(:,n) = (x(:,n1) + x(:,n2) )/2d0 !compute center
   t(1) = -x(2,n2) + x(2,n1)          !normal*l
   t(2) =  x(1,n2) - x(1,n1)    !normal*l
   !------------------
   DO i=1,2                    !sums on average normal weighted by the area
     n1 = surf%lcseg(i,n)      !node
     DO                        !loop over list
       IF( surf%lcnod(j) == n1 )EXIT
       j = j+1
       IF( j > ncnod ) j = 1
     END DO
     tn(:,j) = tn(:,j) + t     !sums averaged normal
   END DO
   !n1 = nr(n1)
   !tn(:,n1) = tn(:,n1) + t     !sums averaged normal
   !n2 = nr(n2)
   !tn(:,n2) = tn(:,n2) + t     !sums averaged normal
   !------------------
   long = SQRT(DOT_PRODUCT(t,t))     !length
   tc(:,n) = t/long                  !unit vector at segment
 END DO
 ! compute unitary normals
 DO n=1,ncnod
   t = tn(:,n)                       !weighted normal at nodes
   long = SQRT(DOT_PRODUCT(t,t))     !length
   tn(:,n) = t/long                  !unit vector
 END DO

 RETURN
 END SUBROUTINE cxctc2
