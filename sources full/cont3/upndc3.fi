  SUBROUTINE upndc3(surf,oldlb,flag,keep)
!
!     Modify internal denomination of nodes when they have been modified
!     by a new strategy
!
  IMPLICIT NONE
  TYPE (surf3_db):: surf          ! INTENT(IN OUT)
  INTEGER (kind=4), INTENT(IN) :: oldlb(:)  !original labels
  LOGICAL, INTENT(OUT) :: flag, & !changes have occured in segments or nodes
                          keep    !segments or nodes remains in the model

  ! local variables
  INTEGER (kind=4) :: ncnod,nsegm,i,j,n,lab, chnode,nn,nel
  INTEGER (kind=4), POINTER :: newn(:),news(:,:)
  INTEGER (kind=4), PARAMETER :: nnseg = 3  !No of nodes per segment

  ncnod = surf%ncnod            !No of nodes in the surface
  nsegm = surf%nsegm            !No of segments defining the surface

  flag = .FALSE.                !initializes to NO CHANGES
  IF( ASSOCIATED( surf%lcseg ) )THEN    !if connectivities are given
    ! ... change nodes in connectivities
    nel = 1                     !initializes
    OUTER : DO j=1,nsegm          !for each segment
      DO i=1,nnseg                  !for each node in the segment
        n = surf%lcseg(i,j)         !old internal number
        lab = oldlb(n)              !original label
        lab = chnode(lab)           !new internal number
        IF( lab > 0 )THEN           !if node exists
          surf%lcseg(i,nel) = lab   !new internal number
        ELSE
          CYCLE OUTER               !a node does not exists, go to next segment
        END IF
      END DO
      nel = nel + 1                 !increase number of remaining elements
    END DO OUTER
    nel = nel - 1               !update final number
    IF( nel < nsegm )THEN       !if number of elements have changed
      flag = .TRUE.               !change flag to CHANGES
      surf%nsegm = nel            !update number of segments
      IF( nel > 0 )THEN              !if Surface still exists
        ALLOCATE( news(nnseg,nel) )    !get memory for new array
        news = surf%lcseg(:,1:nel)     !pass data
        DEALLOCATE(surf%lcseg)         !release old array
        surf%lcseg => news             !point to new data
      ELSE
        DEALLOCATE(surf%lcseg)       !surface does not exist, release memory
      END IF
    END IF
    keep = nel > 0      !does surface exist?
  END IF
  !
  IF( ASSOCIATED( surf%lcnod ) )THEN  !IF list of nodes is active
     ! ... change list of nodes
    nn = 0                   !initializes
    DO j=1,ncnod                 !for each node
      n = surf%lcnod(j)          !old internal number
      lab = oldlb(n)             !original label
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
        ALLOCATE( newn(nn) )         !get memory for new array
        newn = surf%lcnod(1:nn)      !pass data
        DEALLOCATE(surf%lcnod)       !release old array
        surf%lcnod => newn           !point to new data
      ELSE
        DEALLOCATE(surf%lcnod)     !no data remains, release array
      END IF
    END IF
    keep = nn > 0      !does surface exist?
  END IF

  RETURN
  END SUBROUTINE upndc3
