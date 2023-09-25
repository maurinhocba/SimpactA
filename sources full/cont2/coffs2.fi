 SUBROUTINE coffs2(offset,x,surf,np)
 !
 ! move surface position along the normal using offset
 !
 IMPLICIT NONE
 REAL (KIND=8), INTENT(IN) :: offset  !offset value
 INTEGER (KIND=4), INTENT (IN) :: np  !number of nodes in the model
 REAL (KIND=8), POINTER :: x(:,:)     !IN OUT
 TYPE (surf2_db) :: surf              !connectivities

 !local variables
 INTEGER (KIND=4) :: ncnod,nsegm,i,k,n,nn(2)
 REAL (KIND=8), POINTER :: tn(:,:)
 REAL (KIND=8) :: ls(2),t(2),long


 ncnod = surf%ncnod            !number of nodes in the surface
 nsegm = surf%nsegm            !number of segments defining surface

 ALLOCATE(tn(2,np))  !auxiliar array

 tn = 0d0                      !average normals initialization

 ! Compute average normals
 DO n=1,nsegm                  !for each segment
   nn(:) = surf%lcseg(:,n)     !connectivities of the segment
   ls(:) = coora(:,nn(2)) - coora(:,nn(1))  !segment direction
   t  = (/ -ls(2), ls(1) /)    !normal*length
   long = DOT_PRODUCT(t,t)     !squared length
   DO i=1,2                    !sums on average normal weighted by the inverse of the area
     k = nn(i)                 !local numeration of the node
     tn(:,k) = tn(:,k) + t/long     !weigthed sum (inverse of length is the weigthing factor)
   END DO
 END DO
 ! displace nodes along average normal
 DO n=1,np               !for each node in the model
   t = tn(:,n)                       !weighted normal
   long = SQRT(DOT_PRODUCT(t,t))     !offset / length
   IF( long == 0d0 )CYCLE
   long = offset/long
   t = t*long                  !unit vector x offset
   x(:,n) = x(:,n) + t         !modify original coordinates
   coora(:,n) = coora(:,n) + t !and present coordinates
 END DO

 DEALLOCATE (tn)
 RETURN
 END SUBROUTINE coffs2
