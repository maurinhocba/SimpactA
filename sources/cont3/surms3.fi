 SUBROUTINE surms3(emass,x)

!.... compute nodal mass associated to surface

 IMPLICIT NONE
!     dummy arguments
 REAL (kind=8), INTENT (IN) :: x(:,:)
 REAL (kind=8), INTENT (IN OUT) :: emass(:,:)
!     local variables
 INTEGER (kind=4) :: i,j,k,isurf,nsegm

 TYPE (surf3_db), POINTER :: surf
 LOGICAL :: smass(nsurf)

!.... compute mass

 smass = .FALSE.                !initializes mass control array
 surf => shead                 !pointer to first surface
 DO isurf=1,nsurf              !for each surface
   nsegm  = surf%nsegm         !number of segments defining the surface
   IF( surf%density > 0d0)  &  !compute mass
     CALL csmas3(surf%lcseg,x,emass,nsegm,surf%density)

   IF( ASSOCIATED(surf%lcnod) )THEN  !list of nodes exist
     DO i=1,surf%ncnod               !for each node in the surface
       j = surf%lcnod(i)             !global node number
       IF(ALL(emass(1:3,j) == 0d0) ) smass(isurf) = .TRUE. !no mass
       IF(smass(isurf))EXIT !unless one node have no mass
     END DO
   ELSE                              !use connectivities
     DO i=1,surf%nsegm               !for each segment in the surface
       DO k=1,3                      !for each node in segment
         j = surf%lcseg(k,i)         !global node number
         IF(ALL(emass(1:3,j) == 0d0)) smass(isurf) = .TRUE. !no mass
       END DO
       IF(smass(isurf))EXIT !unless one node have no mass
     END DO
   END IF
   surf => surf%next       !pointer to next surface
 END DO
 RETURN
 END SUBROUTINE surms3
