  SUBROUTINE upndc1(surf,oldlb,flag)
!
!     Modify internal denomination of nodes when they have been modified
!     by a new strategy
!
  IMPLICIT NONE
  TYPE (surf1_db):: surf  ! INTENT(IN OUT)
  INTEGER (kind=4), INTENT(IN) :: oldlb(:)
  LOGICAL, INTENT(OUT) :: flag

  ! local variables
  INTEGER (kind=4) :: ncnod,nsegm,i,j,n,lab, chnode,nn,nel
  INTEGER (kind=4), POINTER :: newn(:),news(:,:)
  INTEGER (kind=4), PARAMETER :: nnseg = 2  !No of nodes per segment

  ncnod = surf%ncnod            !No of nodes in the surface
  nsegm = surf%nsegm            !No of segments defining the surface

  flag = .FALSE.
  IF( ASSOCIATED( surf%lcseg ) )THEN
    ! ... change nodes in connectivities
    nel = 1
    OUTER : DO j=1,nsegm                !for each segment
      DO i=1,nnseg              !for each node
        n = surf%lcseg(i,j)      !old internal number
        lab = oldlb(n)           !label
        lab = chnode(lab) !new internal number
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
  END IF
  IF( ASSOCIATED( surf%lcnod ) )THEN
     ! ... change list of nodes
    nn = 0
    DO j=1,ncnod                 !for each node
      n = surf%lcnod(j)          !old internal number
      lab = oldlb(n)             !label
      lab = chnode(lab)
      IF( lab > 0 )THEN
        nn = nn+1
        surf%lcnod(nn) = lab !new internal number
      END IF
    END DO
    IF( nn < ncnod )THEN
      flag = .TRUE.
      surf%ncnod = nn
      IF( nn > 0 )THEN
        ALLOCATE( newn(nn) )
        newn = surf%lcnod(1:nn)
        DEALLOCATE(surf%lcnod)
        surf%lcnod => newn
      ELSE
        DEALLOCATE(surf%lcnod)
      END IF
    END IF
  END IF

  RETURN
  END SUBROUTINE upndc1
