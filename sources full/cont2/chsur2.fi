 SUBROUTINE chsur2(surf,found,flag,lcseg)
 !
 !  update surface from an element set
 !
 IMPLICIT NONE
 ! routine arguments
 TYPE (surf2_db), POINTER :: surf
 INTEGER(kind=4), POINTER :: lcseg(:,:)
 LOGICAL :: found,flag

 ! local variables
 INTEGER (kind=4), PARAMETER :: nnseg = 2
 INTEGER (kind=4) :: nsegm,ncnod,ncnxx,ipair,i
 INTEGER (kind=4), POINTER :: lcnod(:)
 TYPE (pair2_db), POINTER :: pair

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

   ! element connections
    IF ( ASSOCIATED (surf%nhseg) )THEN
      DEALLOCATE ( surf%nhseg )
      IF( nsegm > 0 ) THEN
        ALLOCATE ( surf%nhseg(nnseg,nsegm) )      ! get memory for side elmts.
        CALL mastdb (nsegm,nnseg,surf%lcseg,surf%nhseg) ! get side elements
      END IF
    END IF

    ! reversed connectivities
    IF ( ASSOCIATED (surf%nhseb) ) THEN
      DEALLOCATE ( surf%nhseb ); NULLIFY(surf%nhseb)
      IF( nsegm > 0 )THEN
        ALLOCATE ( surf%nhseb(nnseg,nsegm) )    ! get memory for side elem.
        surf%nhseb(1,:) = surf%nhseg(2,:)       ! assign in reversed order
        surf%nhseb(2,:) = surf%nhseg(1,:)
      END IF
    END IF

    ! reverse connections
    IF ( ASSOCIATED (surf%lcseb) ) THEN
      DEALLOCATE ( surf%lcseb )
      IF( nsegm > 0 )THEN
        ALLOCATE ( surf%lcseb(nnseg,nsegm) )    ! get memory for connections
        surf%lcseb(1,:) = surf%lcseg(2,:)       ! assign in reversed order
        surf%lcseb(2,:) = surf%lcseg(1,:)
      END IF
    END IF

    ! other arrays

    !  element center coodinates
    IF ( ASSOCIATED (surf%xc) ) THEN
      DEALLOCATE ( surf%xc )
      IF( nsegm > 0 ) THEN
        ALLOCATE ( surf%xc(2,nsegm) )             ! get memory for centers
      END IF
    END IF

    !! inverse relation
    !IF ( ASSOCIATED (surf%nr ) ) THEN
    !  DEALLOCATE ( surf%nr  )
    !  IF( nsegm > 0 )THEN
    !    ALLOCATE( surf%nr(npoin))
    !    surf%nr = 0          !generate inverse relation
    !    DO i=1,surf%ncnod
    !      surf%nr(surf%lcnod(i)) = i
    !    END DO
    !  END IF
    !END IF

    IF ( ASSOCIATED (surf%cu ) )THEN
      DEALLOCATE ( surf%cu )
      IF( nsegm > 0 )THEN
        IF( .NOT.ASSOCIATED (surf%lcnod) )THEN
          ncnxx = 2*nsegm+2                         ! maximum number of nodes
          ALLOCATE ( lcnod( ncnxx ) )               ! temporary auxiliar space
          ncnod = 0                                 ! initializes
          CALL getnod(ncnod,ncnxx,nsegm,nnseg,lcnod,lcseg)
          surf%ncnod = ncnod                        ! new number of nodes
          ALLOCATE( surf%lcnod(ncnod) )             ! get memory for list of nodes
          surf%lcnod = lcnod(1:ncnod)               ! assign list of nodes
          DEALLOCATE ( lcnod )            ! free auxiliar space
        END IF
        ALLOCATE(surf%cu(nsegm))
        CALL curve2(coora,surf,npoin)
      END IF
    END IF

    IF ( ASSOCIATED (surf%tc   ) )THEN
      DEALLOCATE ( surf%tc    )
      IF( nsegm > 0 ) THEN
        ALLOCATE(surf%tc(2,nsegm))
        DO i=1,nsegm
          surf%tc (:,i) = (/0d0,1d0/)
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
 END SUBROUTINE chsur2
