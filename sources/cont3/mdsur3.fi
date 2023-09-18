 SUBROUTINE mdsur3(sname)
 !
 !     Modifies an existing surface Data_Base extracted from element data
 !
 IMPLICIT NONE
 ! routine arguments
 CHARACTER(len=*),INTENT(IN):: sname            !surface name

 ! local variables
 LOGICAL :: found
 INTEGER (kind=4) :: nsegm,ncnod,ncnxx,ipair,newn,i,j
 INTEGER (kind=4), PARAMETER :: nnseg = 3
 INTEGER (kind=4), ALLOCATABLE :: lcseg(:,:),lcnod(:)
 INTEGER (kind=4), POINTER :: issdb(:,:)
 REAL (kind=8), POINTER :: rssdb(:,:)
 TYPE (pair3_db), POINTER :: pair
 TYPE (surf3_db), POINTER :: surf,psurf

 INTERFACE
   INCLUDE 'elemnt.h'
 END INTERFACE

 ! search for surface named SNAME in actual database
 CALL srch3_srf(shead, psurf, surf, sname, found)
 IF (.NOT.found) RETURN  !surface SNAME not used for contact

 CALL new_srf(surfa)         ! allocate auxiliar surface (in SURF_DB)

 CALL elemnt('SURFAC',name=sname,flag2=found)  ! get surface from
                                               ! element data base

 nsegm = surfa%nelem*2    ! maximum number of segms

 ALLOCATE(lcseg(nnseg,nsegm))                  ! temporary auxiliar space
 CALL store_segs(surfa%head,lcseg,nnseg,nsegm) !surface ==> lcseg
 CALL dalloc_srf(surfa)                        ! release surface

 !                      for surfaces acting as SLAVE
 IF( surf%iscod )THEN
   ncnxx = 2*nsegm+2                         ! maximum number of nodes
   ALLOCATE ( lcnod( ncnxx ) )               ! temporary auxiliar space
   ncnod = 0                                 ! initializes
   ! get new nodes and keep old nodes
   CALL getno3(nsegm,nnseg,lcseg,surf%ncnod,surf%lcnod,ncnod,lcnod,newn)

   ! Updates pairs associated to SLAVE surface
   pair => headp                               ! initializes
   DO ipair=1,npair                            ! for each pair
     IF( TRIM(pair%slave) == TRIM(sname) )THEN             ! if the right surface
       ALLOCATE ( issdb(nisdb,ncnod), rssdb(nrsdb,ncnod) )
       issdb(:,1:newn) = 0                     ! initializes database
       rssdb(:,1:newn) = 0d0                   ! of new nodes
       DO i=1,surf%ncnod                       ! for each new node
         j = surf%lcnod(i)                     ! position in new list
         IF( j /= 0 )THEN                      ! if node exists
           issdb(:,j) = pair%issdb(:,i)        ! transfer old data
           rssdb(:,j) = pair%rssdb(:,i)
         END IF
       END DO
       DEALLOCATE ( pair%issdb, pair%rssdb )   ! release memory
       pair%issdb => issdb                     ! transfer new arrays
       pair%rssdb => rssdb
       NULLIFY ( issdb, rssdb )                ! nullify pointers
       pair%ncnod = ncnod                      ! new number of nodes
       pair%prev = .FALSE.
       IF(pair%press) THEN
         DEALLOCATE( pair%presn )
         ALLOCATE( pair%presn(pair%ncnod) )
         pair%presn = 0d0                      ! initializes
       END IF
       IF( pair%wrink ) THEN
         DEALLOCATE( pair%mingp )
         ALLOCATE( pair%mingp(pair%ncnod) )
         pair%mingp = 1d9                      !initializes
       END IF
     END IF
     pair => pair%next                         ! next pair
   END DO
   surf%ncnod = ncnod                        ! new number of nodes
   DEALLOCATE( surf%lcnod )                  ! release memory
   ALLOCATE( surf%lcnod(ncnod) )             ! get memory for list of nodes
   surf%lcnod = lcnod(1:ncnod)               ! assign list of nodes
   DEALLOCATE ( lcnod )                      ! free auxiliar space
 ELSE
   surf%ncnod = 0                            ! nodes list not used
 END IF

 ! Deallocate current pointers if exist
 IF ( ASSOCIATED (surf%lcseg) ) DEALLOCATE ( surf%lcseg )
 IF ( ASSOCIATED (surf%nhseg) ) DEALLOCATE ( surf%nhseg )
 IF ( ASSOCIATED (surf%lcseb) ) DEALLOCATE ( surf%lcseb )
 IF ( ASSOCIATED (surf%nhseb) ) DEALLOCATE ( surf%nhseb )
 IF ( ASSOCIATED (surf%xc   ) ) DEALLOCATE ( surf%xc    )

 IF( surf%imcod .OR. surf%iwrit .OR. surf%press )THEN ! If connections required
   surf%nsegm = nsegm                          ! new number of elements
   ALLOCATE ( surf%lcseg(nnseg,nsegm) )        ! get memory for connectivities
   surf%lcseg = lcseg(1:nnseg,1:nsegm)         ! assign connectivities
   IF( surf%imcod ) THEN                       ! if acts as MASTER
     ALLOCATE ( surf%nhseg(nnseg,nsegm) )      ! get memory for side elmts.
     CALL mastd3 (nsegm,nnseg,surf%lcseg,surf%nhseg) ! get side elements
     ALLOCATE ( surf%xc(3,nsegm) )             ! get memory for centers
     IF(surf%bottom )THEN                      ! if bottom surface is used
       ALLOCATE ( surf%lcseb(nnseg,nsegm) )    ! get memory for connections
       surf%lcseb(1,:) = surf%lcseg(1,:)       ! assign in reversed order
       surf%lcseb(2,:) = surf%lcseg(3,:)
       surf%lcseb(3,:) = surf%lcseg(2,:)
       ALLOCATE ( surf%nhseb(nnseg,nsegm) )    ! get memory for side elem.
       surf%nhseb(1,:) = surf%nhseg(1,:)       ! assign in reversed order
       surf%nhseb(2,:) = surf%nhseg(3,:)
       surf%nhseb(3,:) = surf%nhseg(2,:)
     END IF
   END IF
   ! Updates pairs associated to SLAVE surface
   pair => headp                               ! initializes
   DO ipair=1,npair                            ! for each pair
     IF( TRIM(pair%master) == TRIM(sname) )THEN            ! if the right surface
       pair%issdb = 0                          ! nullify data base
       pair%rssdb = 0d0
       pair%prev = .FALSE.
       IF(pair%press) pair%presn = 0d0         ! initializes
     END IF
     pair => pair%next                         ! next pair
   END DO
 ELSE
   surf%nsegm = 0
 END IF

 DEALLOCATE ( lcseg )                        ! free auxiliar space

 RETURN
 END SUBROUTINE mdsur3

