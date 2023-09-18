  SUBROUTINE upndc2(surf,oldlb,flag,keep)
!
!     Modify internal denomination of nodes when they have been modified
!     by a new strategy
!
  IMPLICIT NONE
  TYPE (surf2_db):: surf  ! INTENT(IN OUT)
  INTEGER (kind=4), INTENT(IN) :: oldlb(:)
  LOGICAL, INTENT(OUT) :: flag,keep

  ! local variables
  INTEGER (kind=4) :: ncnod,nsegm,i,j,n,lab, chnode,nn,nel
  INTEGER (kind=4), POINTER :: newn(:),news(:,:)
  INTEGER (kind=4), PARAMETER :: nnseg = 2  !No of nodes per segment

  ncnod = surf%ncnod            !No of nodes in the surface
  nsegm = surf%nsegm            !No of segments defining the surface

  flag = .FALSE.                !initializes to NO changes
  IF( ASSOCIATED( surf%lcseg ) )THEN    !if connectivities are given
    ! ... change nodes in connectivities
    nel = 1                     !initializes
    OUTER : DO j=1,nsegm          !for each segment
      DO i=1,nnseg                  !for each node
        n = surf%lcseg(i,j)         !old internal number
        lab = oldlb(n)              !label
        lab = chnode(lab)           !new internal number
        IF( lab > 0 )THEN
          surf%lcseg(i,nel) = lab !new internal number
        ELSE
          CYCLE OUTER
        END IF
      END DO
      nel = nel + 1
    END DO OUTER
    nel = nel - 1
    IF( nel < nsegm )THEN
      flag = .TRUE.
      surf%nsegm = nel
      IF( nel > 0 )THEN
        ALLOCATE( news(nnseg,nel) )
        news = surf%lcseg(:,1:nel)
        DEALLOCATE(surf%lcseg)
        surf%lcseg => news
      ELSE
        DEALLOCATE(surf%lcseg)
      END IF
    END IF
    keep = nel > 0      !does surface exist?
  END IF
  IF( ASSOCIATED( surf%lcnod ) )THEN
     ! ... change list of nodes
    nn = 0                   !initializes
    DO j=1,ncnod                 !for each node
      n = surf%lcnod(j)          !old internal number
      lab = oldlb(n)             !label
      lab = chnode(lab)          !new internal node
      IF( lab > 0 )THEN            !if node exist
        nn = nn+1                  !update counter
        surf%lcnod(nn) = lab       !new internal number
      END IF
    END DO
    IF( nn < ncnod )THEN         !if nodes have been deleted
      flag = .TRUE.                !update flag
      surf%ncnod = nn              !update number of nodes
      IF( nn > 0 )THEN               !if nodes remain
        ALLOCATE( newn(nn) )
        newn = surf%lcnod(1:nn)
        DEALLOCATE(surf%lcnod)
        surf%lcnod => newn
      ELSE
        DEALLOCATE(surf%lcnod)
      END IF
    END IF
    keep = nn > 0      !does surface exist?
  END IF

  RETURN
  END SUBROUTINE upndc2
