 SUBROUTINE chsur1(surf,found,flag,lcseg)
 !
 !  update surface from an element set
 !
 IMPLICIT NONE
 ! routine arguments
 TYPE (surf1_db), POINTER :: surf
 INTEGER(kind=4), POINTER :: lcseg(:,:)
 LOGICAL :: found,flag

 ! local variables
 INTEGER (kind=4), PARAMETER :: nnseg = 2
 INTEGER (kind=4) :: nsegm,ncnod,ncnxx,ipair
 INTEGER (kind=4), POINTER :: lcnod(:)
 TYPE (pair1_db), POINTER :: pair

 INTERFACE
   INCLUDE 'getnod.h'
 END INTERFACE

   nsegm = surf%nsegm    !number of segments

   !  nodes
   IF( ASSOCIATED(surf%lcnod) .AND. ASSOCIATED(lcseg) )THEN
     DEALLOCATE ( surf%lcnod )  ! free list
     IF( nsegm > 0 )THEN
       ncnxx = 2*nsegm+2                         ! maximum number of nodes
       ALLOCATE ( lcnod( ncnxx ) )               ! temporary auxiliar space
       ncnod = 0                                 ! initializes
       CALL getnod(ncnod,ncnxx,nsegm,nnseg,lcnod,lcseg)
       surf%ncnod = ncnod                        ! new number of nodes
       ALLOCATE( surf%lcnod(ncnod) )             ! get memory for list of nodes
       surf%lcnod = lcnod(1:ncnod)               ! assign list of nodes
       DEALLOCATE ( lcnod )            ! free auxiliar space
     END IF
   END IF

   ! standard connectivities
   IF ( found )THEN
     IF( ASSOCIATED (surf%lcseg) ) THEN
       IF( flag )THEN
         DEALLOCATE ( surf%lcseg )
         IF( nsegm > 0 )ALLOCATE ( surf%lcseg(nnseg,nsegm))
       END IF
       IF( nsegm > 0 )surf%lcseg = lcseg(:,1:nsegm)
     END IF
   END IF

   IF( .NOT.surf%imcod .AND. .NOT.surf%iwrit ) surf%nsegm = 0  !keep as flag

   IF( surf%iscod )THEN  ! if surface acts as slave
     ! Updates pairs associated to SLAVE surface
     pair => headp                               ! initializes
     DO ipair=1,npair                            ! for each pair
       IF( TRIM(pair%slave) == TRIM(surf%sname) )THEN             ! if the right surface
         DEALLOCATE ( pair%issdb, pair%rssdb )   ! release memory
         pair%ncnod = surf%ncnod                 ! new number of nodes
         IF( pair%ncnod == 0 )THEN
           pair%end = 0d0           !cancel pair
         ELSE
           ALLOCATE ( pair%issdb(nisdb,pair%ncnod), pair%rssdb(nrsdb,pair%ncnod) )
           pair%issdb = 0                          ! initializes database
           pair%rssdb = 0d0
           pair%prev = .FALSE.
           IF(pair%press) THEN
             IF(ASSOCIATED(pair%presn))DEALLOCATE( pair%presn )
             ALLOCATE( pair%presn(pair%ncnod) )
             pair%presn = 0d0                      ! initializes
           END IF
         END IF
       END IF
       pair => pair%next                         ! next pair
     END DO
   END IF

   IF( surf%nsegm == 0 .AND.  ASSOCIATED(surf%lcseg))DEALLOCATE(surf%lcseg)

 RETURN
 END SUBROUTINE chsur1
