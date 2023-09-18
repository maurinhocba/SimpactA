 SUBROUTINE cxctc3(x,surf)
 !
 ! Compute nodal and segment normals for self contact
 !
 IMPLICIT NONE
 REAL (KIND=8), INTENT(IN) :: x(:,:)
 TYPE (surf3_db) :: surf
 ! INTENT(OUT) :: surf%tn(2,ncnod) :nodal normals
 !                surf%tc(2,nsegm) :segment normals for self contact
 !                surf%xc(2,nsegm) :segment position of center
 !local variables
 INTEGER (KIND=4) :: ncnod,nsegm,i,j,n,nn(3)
 REAL (KIND=8) :: t1(3),t2(3),t3(3),long

 INTEGER (kind=4), POINTER :: nr(:)
 REAL (KIND=8), POINTER :: tn(:,:),tc(:,:),xc(:,:)


 ncnod = surf%ncnod            !number of nodes in the surface
 nsegm = surf%nsegm            !number of segments defining surface
 !  point to surface arrays
 xc => surf%xc                 !(3,nsegm)
 tc => surf%tc                 !(3,nsegm)
 tn => surf%tn                 !(3,ncnod)
 nr => surf%nr

 tn = 0d0                      !averaged normals initialization

 ! Compute segment normal and center
 DO n=1,nsegm               !for each segment
   nn = surf%lcseg(:,n)     !nodes
   xc(:,n) = (x(:,nn(1)) + x(:,nn(2)) + x(:,nn(3)) )/3d0 !compute center
   t1 = x(:,nn(2)) - x(:,nn(1))
   t2 = x(:,nn(3)) - x(:,nn(1))
   t3(1) = t1(2)*t2(3) - t1(3)*t2(2)
   t3(2) = t1(3)*t2(1) - t1(1)*t2(3)
   t3(3) = t1(1)*t2(2) - t1(2)*t2(1)
   DO i=1,3                     !sums on average normal weighted by the area
     j = nr(nn(i))              !node
     tn(:,j) = tn(:,j) + t3     !sums averaged normal
   END DO
   long = SQRT(DOT_PRODUCT(t3,t3))   !length
   tc(:,n) = t3/long                  !unit vector at segment
 END DO
 ! compute unitary normals
 DO n=1,ncnod
   t3 = tn(:,n)                      !weighted normal at nodes
   long = SQRT(DOT_PRODUCT(t3,t3))     !length
   tn(:,n) = t3/long                 !unit vector
 END DO

 RETURN
 END SUBROUTINE cxctc3
