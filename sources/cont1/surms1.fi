      SUBROUTINE surms1(emass)

!.... compute nodal mass associated to surface

      IMPLICIT NONE
!     dummy arguments
      REAL (kind=8), INTENT (IN OUT)  :: emass(:,:)
!     local variables
      INTEGER (kind=4) :: i,j,k,isurf,nsegm,n,n1,n2
      REAL (kind=8) :: nmass,l1(3)
      TYPE (surf1_db), POINTER :: surf
      LOGICAL :: mass(nsurf)

      ! compute mass

      mass = .FALSE.
      surf => shead                 !pointer to first surface
      DO isurf=1,nsurf              !for each surface
        nsegm  = surf%nsegm         !number of segments defining the surface
        IF( surf%density > 0d0)THEN    !compute mass
          DO n=1,nsegm                       ! for each segment
            n1 = surf%lcseg(1,n)              ! first node
            n2 = surf%lcseg(2,n)              ! second node
            l1 = coord(:,n2) - coord(:,n1)            ! side 3
            CALL vecuni(3,l1,nmass)           ! side length
            nmass = nmass*surf%density/2d0    ! nodal mass = Total_mass/2
            emass(1:3,n1) = emass(1:3,n1) + nmass  !add to each node
            emass(1:3,n2) = emass(1:3,n2) + nmass
          END DO
        END IF
        IF( ASSOCIATED(surf%lcnod) )THEN
          DO i=1,surf%ncnod
            j = surf%lcnod(i)
            IF(ALL(emass(1:3,j) == 0d0) ) mass(isurf) = .TRUE.
            IF(mass(isurf))EXIT
          END DO
        ELSE
          DO i=1,surf%nsegm
            DO k=1,2
              j = surf%lcseg(k,i)
              IF(ALL(emass(1:3,j) == 0d0)) mass(isurf) = .TRUE.
            END DO
            IF(mass(isurf))EXIT
          END DO
        END IF
        surf => surf%next       !pointer to next surface
      END DO

      RETURN
      END SUBROUTINE surms1
