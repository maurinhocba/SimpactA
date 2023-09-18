 SUBROUTINE cupdl2( )
 !
 !     RENUMBER nodes in surface contact Data_Base
 !     Deletes invalid pairs and surfaces
 !
 IMPLICIT NONE
 ! routine arguments

 ! local variables
 INTEGER (kind=4), PARAMETER :: nnseg = 2
 LOGICAL :: found,flag,keep
 INTEGER (kind=4) :: isurf,j,n,lab,chnode
 TYPE (surf2_db), POINTER :: surf,prev
 TYPE (pair2_db), POINTER :: pair,pprev
 INTEGER (kind=4), POINTER :: lcseg(:,:)
 REAL(kind=8), POINTER :: wwean(:)
 CHARACTER(len=mnam):: name

 INTERFACE
   INCLUDE 'elemnt.h'
 END INTERFACE

 ! Updates surface connectivities

 surf => shead                               ! initializes
 NULLIFY(prev)                               ! initializes
 DO                                          ! for each existing surface
   IF( .NOT. ASSOCIATED(surf)) EXIT          ! if all surface processed EXIT
   keep = .FALSE.                            ! initializes
   CALL new_srf(surfa)                       ! allocates surface in surf_db
   CALL elemnt ('SURFAC', name=surf%sname, flag2=found)  ! get surface from element database
   IF( found )THEN                    ! is an element set
     ALLOCATE( lcseg(nnseg,2*surfa%nelem) )  ! temporary auxiliar space
     CALL store_segs (surfa%head,lcseg,nnseg,surfa%nelem) !surface ==> lcseg
     flag = surfa%nelem /= surf%nsegm        !surface have changed size
     surf%nsegm = surfa%nelem
     keep = .TRUE.
   ELSE                                      ! a surface defined as is (or a deleted surface)
     CALL upndc2(surf,oldlb,flag,keep)       ! update nodal numeration
   END IF
   IF( keep ) THEN
     CALL chsur2(surf,found,flag,lcseg)
     prev => surf
     stail => surf
     surf => surf%next                         ! next surface
   ELSE
     name = surf%sname
     nsurf = nsurf - 1
     CALL del2_srf (shead, stail, prev, surf)
     CALL delete_name(name,5)
     IF( ASSOCIATED(prev))THEN
       surf => prev%next
     ELSE IF( ASSOCIATED(shead))THEN
       surf => shead
     ELSE
       EXIT
     END IF
   END IF
   CALL dalloc_srf(surfa)
   IF(ASSOCIATED(lcseg)) DEALLOCATE ( lcseg )     ! free auxiliar space

 END DO

 pair => headp                               ! initializes
 NULLIFY(pprev)                              ! initializes
 DO                                          ! for each existing pairace
   IF( .NOT. ASSOCIATED(pair)) EXIT          ! if all pairace processed EXIT
   surf => shead                     ! point to first surface
   keep = .FALSE.                    ! check is master surface yet exists
   DO isurf=1,nsurf
     IF( surf%sname == pair%master )THEN
       keep  = .TRUE.
       pair%imast = isurf
       EXIT
     END IF
     surf => surf%next
   END DO
   IF( keep )THEN                    ! is master exist
     surf => shead                     ! point to first surface
     keep = .FALSE.                  ! check is slave surface yet exists
     DO isurf=1,nsurf
       IF( surf%sname == pair%slave )THEN
         keep = .TRUE.
         pair%islav = isurf
         EXIT
       END IF
       surf => surf%next
     END DO
   END IF
   IF( keep ) THEN
     pprev => pair                             ! keep previous
     tailp => pair                             ! update tail
     pair => pair%next                         ! next pair
   ELSE                            !delete pair
     name = pair%pname
     CALL del_pair2 (headp, tailp, pprev , pair)
     npair = npair - 1   !decrement number of existing pairs
     IF( ASSOCIATED(pprev))THEN
       pair => pprev%next
     ELSE IF( ASSOCIATED(headp))THEN
       pair => headp
     ELSE
       EXIT
     END IF
     CALL delete_name(name,4)
   END IF

 END DO

 nconp (2) = npair ! put number of contact pairs of cont. type 2

 ! check if active surface do exist in any pair

 surf => shead                               ! initializes
 NULLIFY(prev)                               ! initializes
 DO                                          ! for each existing surface
   IF( .NOT. ASSOCIATED(surf)) EXIT          ! if all surface processed EXIT
   keep = .FALSE.                            ! initializes
   IF( npair > 0 )THEN
     pair => headp                             ! point to first pair
     DO
       IF(.NOT.ASSOCIATED(pair)) EXIT
       keep = TRIM(surf%sname) == TRIM(pair%slave) .OR. TRIM(surf%sname) == TRIM(pair%master)
       IF(keep)EXIT
       pair => pair%next
     END DO
     IF( keep ) THEN
       prev => surf
       stail => surf
       surf => surf%next                         ! next surface
     END IF
   ELSE
     ! ERROR, this modifies internal numeration updated above
     name = surf%sname
     nsurf = nsurf - 1
     CALL del2_srf (shead, stail, prev, surf)
     CALL delete_name(name,5)
     IF( ASSOCIATED(prev))THEN
       surf => prev%next
     ELSE IF( ASSOCIATED(shead))THEN
       surf => shead
     ELSE
       EXIT
     END IF
     EXIT
   END IF

 END DO

 ! change node numeration for friction work
 IF( wear)THEN
   IF( npair > 0)THEN
     ALLOCATE ( wwean(npoin) )
     wwean = 0d0
     DO j=1,npoio               !for each node
       lab = oldlb(j)           !label
       n = chnode(lab)          !new internal number
       wwean(n) = wwear(j)
     END DO
     DEALLOCATE( wwear); NULLIFY(wwear)
     wwear => wwean
   ELSE
     wear = .FALSE.
     DEALLOCATE( wwear )
   END IF
 END IF

 RETURN
 END SUBROUTINE cupdl2
