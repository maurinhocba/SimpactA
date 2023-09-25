 SUBROUTINE chsur3(surf,found,flag,lcseg)
 !
 !  update surface from an element set
 !
 IMPLICIT NONE
 ! routine arguments
 TYPE (surf3_db), POINTER :: surf        !old surface
 INTEGER(kind=4), POINTER :: lcseg(:,:)  !new connectivities
 LOGICAL :: found, & !surface comes from  an element set
            flag     !surface has changed from previous strategy

 ! local variables
 INTEGER (kind=4), PARAMETER :: nnseg = 3
 INTEGER (kind=4) :: nsegm,ncnod,ncnxx,ipair,i
 INTEGER (kind=4), POINTER :: lcnod(:)
 TYPE (pair3_db), POINTER :: pair

 INTERFACE
   INCLUDE 'getnod.h'
 END INTERFACE


   nsegm = surf%nsegm    !number of segments

   !  nodes
   IF( ASSOCIATED(surf%lcnod) .AND. ASSOCIATED(lcseg) )THEN !previus nodes and present segments exist
     ! regenerate list of nodes from new set of connectivities
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
   IF ( found )THEN        !if segments are part of an element set
     IF( ASSOCIATED (surf%lcseg) ) THEN
       IF( flag )THEN      !if surface have changed
         DEALLOCATE ( surf%lcseg )   !release old connectivities
         IF( nsegm > 0 )ALLOCATE ( surf%lcseg(nnseg,nsegm))  !if new connectivities exist, get memory
       END IF
       IF( nsegm > 0 )THEN       !if number of segments is possitive
         DO i=1,nsegm
           surf%lcseg(:,i) = lcseg(:,i)  !pass new connectivities
         END DO
       END IF
     END IF
   END IF

   ! element connections
    IF ( ASSOCIATED (surf%nhseg) )THEN  !if side connectios are necessary
      DEALLOCATE ( surf%nhseg )            !release array
      IF( nsegm > 0 ) THEN                     !if connectivities are present
        ALLOCATE ( surf%nhseg(nnseg,nsegm) )            ! get memory for side elmts.
        CALL mastd3 (nsegm,nnseg,surf%lcseg,surf%nhseg) ! get side elements
      END IF
    END IF

    ! reversed connectivities
    IF ( ASSOCIATED (surf%nhseb) ) THEN   !if reverse connectivities are necessary
      DEALLOCATE ( surf%nhseb ); NULLIFY(surf%nhseb) !release old data
      IF( nsegm > 0 )THEN                 !if connectivities are present
        ALLOCATE ( surf%nhseb(nnseg,nsegm) )    ! get memory for side elem.
        DO i=1,nsegm
          surf%nhseb(1,i) = surf%nhseg(1,i)       ! assign in reversed order
          surf%nhseb(2,i) = surf%nhseg(3,i)
          surf%nhseb(3,i) = surf%nhseg(2,i)
        END DO
      END IF
    END IF

    ! reverse connections
    IF ( ASSOCIATED (surf%lcseb) ) THEN   !if side (reverse) connectios are necessary
      DEALLOCATE ( surf%lcseb )           !release old data
      IF( nsegm > 0 )THEN                    !if connectivites are present
        ALLOCATE ( surf%lcseb(nnseg,nsegm) )    ! get memory for connections
        DO i=1,nsegm
          surf%lcseb(1,i) = surf%lcseg(1,i)       ! assign in reversed order
          surf%lcseb(2,i) = surf%lcseg(3,i)
          surf%lcseb(3,i) = surf%lcseg(2,i)
        END DO
      END IF
    END IF

    ! OTHER ARRAYS

    !  element center coodinates
    IF ( ASSOCIATED (surf%xc) ) THEN
      DEALLOCATE ( surf%xc )
      IF( nsegm > 0 ) THEN
        ALLOCATE ( surf%xc(3,nsegm) )             ! get memory for centers
        ! recomputed later
      END IF
    END IF

    ! inverse relation
    IF ( ASSOCIATED (surf%nr ) ) THEN
      DEALLOCATE ( surf%nr  )
      IF( nsegm > 0 )THEN
        ALLOCATE( surf%nr(npoin))
        surf%nr = 0          !generate inverse relation
        DO i=1,surf%ncnod
          surf%nr(surf%lcnod(i)) = i
        END DO
      END IF
    END IF

    ! element curvatures
    IF ( ASSOCIATED (surf%cu ) )THEN
      DEALLOCATE ( surf%cu )
      IF( nsegm > 0 )THEN
        ! here surf%lcnod is undefined because is unnecessary for master surface
        ! so an error occurs when curvatures and normals are computed
        IF( .NOT.ASSOCIATED (surf%lcnod) )THEN
          ncnxx = 2*nsegm+2                         ! maximum number of nodes
          ALLOCATE ( lcnod( ncnxx ) )               ! temporary auxiliar space
          ncnod = 0                                 ! initializes
          CALL getnod(ncnod,ncnxx,nsegm,nnseg,lcnod,surf%lcseg)
          surf%ncnod = ncnod                        ! new number of nodes
          ALLOCATE( surf%lcnod(ncnod) )             ! get memory for list of nodes
          surf%lcnod = lcnod(1:ncnod)               ! assign list of nodes
          DEALLOCATE ( lcnod )            ! free auxiliar space
        END IF
        ALLOCATE(surf%cu(3,nsegm))
        CALL curve3(coora,surf,npoin)
      END IF
    END IF
    ! element normals
    IF ( ASSOCIATED (surf%tc   ) )THEN
      DEALLOCATE ( surf%tc    )
      IF( nsegm > 0 ) THEN
        ALLOCATE(surf%tc(3,nsegm))
        DO i=1,nsegm
          surf%tc (:,i) = (/0d0,0d0,1d0/)
        END DO
      END IF
    END IF

   IF( .NOT.surf%imcod .AND. .NOT.surf%iwrit .AND. .NOT.surf%press ) surf%nsegm = 0  !keep as flag


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
           IF( pair%wrink ) THEN
             IF(ASSOCIATED(pair%mingp))DEALLOCATE( pair%mingp )
             ALLOCATE( pair%mingp(pair%ncnod) )
             pair%mingp = 1d9                      !initializes
           END IF
         END IF
       END IF
       pair => pair%next                         ! next pair
     END DO
   END IF

   IF( surf%nsegm == 0 .AND.  ASSOCIATED(surf%lcseg))DEALLOCATE(surf%lcseg)

 RETURN
 END SUBROUTINE chsur3
