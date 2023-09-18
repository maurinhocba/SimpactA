 SUBROUTINE cupdl3( )
 !
 !     RENUMBER nodes in surface contact Data_Base
 !     Deletes invalid pairs and surfaces
 !
 IMPLICIT NONE
 ! routine arguments

 ! local variables
 INTEGER (kind=4), PARAMETER :: nnseg = 3  !number of nodes per segment
 LOGICAL :: found,flag,keep                !flagas
 INTEGER (kind=4) :: isurf,j,n,lab,chnode  !index
 TYPE (surf3_db), POINTER :: surf,prev     !auxiliar surfaces
 TYPE (pair3_db), POINTER :: pair,pprev    !auxiliar pairs
 INTEGER (kind=4), POINTER :: lcseg(:,:)   !auxiliar surface connectivities
 REAL(kind=8), POINTER :: wwean(:)         !
 CHARACTER(len=mnam):: name                !just a name

 INTERFACE
   INCLUDE 'elemnt.h'
 END INTERFACE

 ! Updates surface connectivities, internal numeration

 surf => shead                               ! initializes to first surface
 NULLIFY(prev)                               ! initializes previous to none
 DO                                          ! for each existing surface
   NULLIFY( lcseg)                           ! initializes connectivities pointer
   IF( .NOT. ASSOCIATED(surf)) EXIT          ! if all surface processed EXIT
   keep = .FALSE.                            ! initializes not to keep
   CALL new_srf(surfa)                       ! allocates surface in surf_db
   CALL elemnt ('SURFAC', name=surf%sname, flag2=found)  ! get surface from element database
   IF( found )THEN                           ! is an element set
     ALLOCATE( lcseg(nnseg,2*surfa%nelem) )  ! temporary auxiliar space
     CALL store_segs (surfa%head,lcseg,nnseg,surfa%nelem) !surface ==> lcseg
     flag = surfa%nelem /= surf%nsegm        !surface have changed size
     surf%nsegm = surfa%nelem                !store number of present segments
     keep = .TRUE.                           !KEEP this surface
   ELSE                                      ! a surface defined as is (or a deleted surface)
     CALL upndc3(surf,oldlb,flag,keep)       ! update nodal numeration
   END IF

   IF( keep ) THEN   !if surface remains
     CALL chsur3(surf,found,flag,lcseg)        !update all the associated arrays
     prev => surf
     stail => surf
     surf => surf%next                         ! next surface
   ELSE              !if surface has dissapeared
     name = surf%sname                         !surface name
     nsurf = nsurf - 1                         !decrease number of active surfaces
     CALL del3_srf (shead, stail, prev, surf)  !deletes surface from list
     CALL delete_name(name,5)                  !delete surface name from database
     IF( ASSOCIATED(prev))THEN                 !if a previous is active
       surf => prev%next                         !point to next surface
     ELSE IF( ASSOCIATED(shead))THEN           !if a top surface still exist (changed in del3_srf)
       surf => shead                             !point to next surface
     ELSE                                      !if no surface remains
       EXIT                                      !EXIT DO
     END IF
   END IF
   CALL dalloc_srf(surfa)                      !release auxiliar surface
   IF(ASSOCIATED(lcseg)) DEALLOCATE ( lcseg )     ! free auxiliar space

 END DO

 ! Updates pair data in terms of remaining surfaces

 pair => headp                       ! point to first pair
 NULLIFY(pprev)                      ! no previous pair in the list
 DO                                     ! for each existing pair
   IF( .NOT. ASSOCIATED(pair)) EXIT     ! if all pair processed EXIT
   keep = .FALSE.                    ! check is master surface yet exists
   surf => shead                     ! point to first remaining surface
   DO isurf=1,nsurf                     !for each surface
     IF( surf%sname == pair%master )THEN   !compare names
       keep  = .TRUE.                    !keep pair
       pair%imast = isurf                !update position of master surface in the list
       EXIT                              !exit loop
     END IF
     surf => surf%next                  !point to next surface if not found yet
   END DO
   IF( keep )THEN                    ! is master exist
     keep = .FALSE.                  ! check is slave surface yet exists
     surf => shead                     ! point to first surface
     DO isurf=1,nsurf                    !for each surface
       IF( surf%sname == pair%slave )THEN  !compare names
         keep = .TRUE.                     !both surfaces exist, KEEP pare
         pair%islav = isurf                !updates position of slave surface in the list
         EXIT                              !exit loop
       END IF
       surf => surf%next                 !point to next surface if not found yet
     END DO
   END IF
   !
   IF( keep ) THEN  !both master and slave surface still exist
     pprev => pair                             ! keep previous
     tailp => pair                             ! update tail
     pair => pair%next                         ! next pair
   ELSE             !delete pair
     name = pair%pname       !keep name to delete from data base
     CALL del_pair3 (headp, tailp, pprev , pair)   !delete pair from list
     npair = npair - 1   !decrement number of existing pairs
     CALL delete_name(name,4)                  !delete pair name from database
     IF( ASSOCIATED(pprev))THEN                !if a previous pair
       pair => pprev%next                        !point to next pair
     ELSE IF( ASSOCIATED(headp))THEN           !if the first pair is the head
       pair => headp                             !point to the head
     ELSE
       EXIT                                    !no pair remain
     END IF
   END IF

 END DO

 nconp (3) = npair ! put number of contact pairs of cont. type 3

! ! check if active surface do exist in any pair
! ! Surfaces not used in any pair are DELETED, this loop is CANCELLED
!
! surf => shead                               ! initializes
! NULLIFY(prev)                               ! initializes
! DO                                          ! for each existing surface
!   IF( .NOT. ASSOCIATED(surf)) EXIT          ! if all surface processed EXIT
!   keep = .FALSE.                            ! initializes
!   IF( npair > 0 )THEN
!     pair => headp                             ! point to first pair
!     DO
!       IF(.NOT.ASSOCIATED(pair)) EXIT
!       keep = TRIM(surf%sname) == TRIM(pair%slave) .OR. TRIM(surf%sname) == TRIM(pair%master)
!       IF(keep)EXIT
!       pair => pair%next
!     END DO
!     IF( keep ) THEN
!       prev => surf
!       stail => surf
!       surf => surf%next                         ! next surface
!     END IF
!   ELSE
!     ! ERROR, this modifies internal numeration updated above
!     name = surf%sname
!     nsurf = nsurf - 1
!     CALL del3_srf (shead, stail, prev, surf)
!     CALL delete_name(name,5)
!     IF( ASSOCIATED(prev))THEN
!       surf => prev%next
!     ELSE IF( ASSOCIATED(shead))THEN
!       surf => shead
!     ELSE
!       EXIT
!     END IF
!     EXIT
!   END IF
!
! END DO

 ! change node numeration for friction work
 IF(wear)THEN
   IF( npair > 0)THEN
     ALLOCATE ( wwean(npoin) )
     wwean = 0d0
     DO j=1,npoio                     !for each node
       lab = oldlb(j)                 !label
       n = chnode(lab)                !new internal number
       IF(n > 0) wwean(n) = wwear(j)
     END DO
     DEALLOCATE( wwear )
     wwear => wwean
   ELSE
     wear = .FALSE.
     DEALLOCATE( wwear )
   END IF
 END IF

 RETURN
 END SUBROUTINE cupdl3
