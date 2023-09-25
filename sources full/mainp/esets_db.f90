 MODULE esets_db

   !  control information for element sets

   USE param_db,ONLY: midn, mset, mnam
   IMPLICIT NONE
   SAVE

   INTEGER (kind=4), PARAMETER :: msets = 25 ! number of element types
   INTEGER (kind=4) :: eset(msets), &        !number of sets of each type
                       esets,       &        !number of type of elements used
                       nel(msets),  &        !number of elements for each set type
                       nelms(msets)          !number of elements sets for this type
   INTEGER (kind=4) :: nrenu=0,     &  !re-numbering flag
                       maxnn,       &  !maximum number of nodes per element
                       melem,       &  !total number of element
                       gelem           !auxiliar counter

   INTEGER (kind=4), POINTER :: gnods(:,:), & !(maxnn,melem) auxiliar array
                                bestn(:)      !(npoin)       auxiliar array

   LOGICAL :: rot_free

   CHARACTER(len=midn) :: ELTY(msets) = (/             &
             'SPOT  ','TRUSS ','NONE  ','SSOLAG','SPRISM',  &
             'SHELQ ','SHELT ','BEAM  ','SHREV ','RIGID ',  &
             'BEAM2 ','NONE  ','NBST  ','LBST  ','RBST  ',  &
             'PRISM ','QUADL ','SOLAG ','TTR2D ','TR2D  ',  &
             'HEAT2 ','HEAT3 ','NONE  ','TLBST ','BSQ   ' /)
   INTEGER(kind=4), PARAMETER :: mnode(msets) = (/ &
                     2, 2, 0, 8, 6, 4, 6, 3, 3, 8, &
                     2, 0, 3, 3, 3, 6, 4, 8, 3, 3, &
                     0, 0, 0, 3, 4  /)


   ! derived type for a set name list
   TYPE set_label
     CHARACTER (len=mnam) :: sname !set name
     INTEGER (kind=4) :: stype   !set type
     TYPE (set_label), POINTER :: next
   END TYPE set_label

   CHARACTER (len=mset), PARAMETER :: snames(8)= (/ 'node set',          & !(1)
                                                    'element set',       & !(2)
                                                    'velocity set',      & !(3)
                                                    'contact pair',      & !(4)
                                                    'contact surface',   & !(5)
                                                    'drawbead line',     & !(6)
                                                    'load set',          & !(7)
                                                    'heat set' /)          !(8)

   INTEGER (kind=4) :: n_names=0               ! number of set names
   TYPE (set_label), POINTER :: snhead,sntail  !first and last names

 CONTAINS

   SUBROUTINE elsets ( )

     !  checking element types used

     IMPLICIT NONE

     INTEGER (kind=4) :: iset
     !--------------------------------------------------

     esets = 0     !initializes
     eset = 0
     melem= 0
     maxnn = 0

     DO iset=1,msets    !for each element type
         IF (nelms(iset) > 0) THEN   !if there are elements for this elm. type
           esets = esets + 1         !increase the number of elm types used
           eset(esets) = iset        !keep for this set the element type
           melem= melem + nel(iset)
           maxnn = MAX(maxnn,mnode(iset))
         END IF
     END DO
     rot_free = (nelms(11) + SUM(nelms(13:15)) + SUM(nelms(24:25))) > 0     !true if rotation free elements exist
     RETURN
   END SUBROUTINE elsets

   SUBROUTINE add_name (name,stype)
    ! insert a name at the end of the list
    IMPLICIT NONE

     !Dummy arguments
     CHARACTER (len=*), INTENT(IN) :: name
     INTEGER (kind=4), INTENT(IN) :: stype

     ! local
     LOGICAL :: found
     TYPE (set_label), POINTER :: new
     INTEGER (kind=4) ::  i

     ! check if name was already used
     i = stype
     CALL search_name(name,found,i)
     IF( found )THEN
       WRITE( *,"(' label ',A,' already used for a ',A)") TRIM(name),TRIM(snames(stype))
       WRITE(55,"(' label ',A,' already used for a ',A)",ERR=9999) TRIM(name),TRIM(snames(stype))
       CALL runend('same name used for different objects')
     END IF

     ALLOCATE (new)                     !get memory and assign
     new%sname = name
     new%stype = stype
     IF ( n_names == 0 ) THEN           !if list is empty
       snhead => new                    !initializes list
       sntail => new                    !point to new name
       n_names = 1
     ELSE
       sntail%next => new               !point to new name
       sntail => new                    !new last name
       n_names = n_names + 1
     END IF

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE add_name

   SUBROUTINE delete_name (name,stype,iostat)
     ! deletes a name fron the list
     ! posic is redirected on the succeeding node, or nullified if end_ns

     !Dummy arguments
     CHARACTER(len=*), INTENT(IN) :: name
     INTEGER (kind=4), INTENT(IN) :: stype
     INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

     ! local
     TYPE (set_label), POINTER :: addr, addr_prev
     INTEGER (kind=4) :: i
     LOGICAL :: flag

     flag = .FALSE.
     addr => snhead                  ! point to first
     NULLIFY(addr_prev)              ! point to nothing
     DO i=1,n_names
       IF( TRIM(addr%sname) == TRIM(name) .AND. addr%stype == stype  )THEN  !name found
         flag = .TRUE.
         IF( i == 1 )THEN            !deleted name is the first (head)
           snhead => addr%next
         ELSE
           addr_prev%next => addr%next
           IF( i == n_names )sntail => addr_prev
         END IF
         IF(PRESENT(iostat))iostat = 1
         EXIT
       ELSE
         addr_prev => addr
         addr => addr_prev%next
       END IF
     END DO
     IF( flag )THEN
       IF(PRESENT(iostat))iostat =  1  ! O.K.
       n_names = n_names - 1
       DEALLOCATE(addr)
     ELSE
       IF(PRESENT(iostat))iostat = -1  ! error
     END IF

   END SUBROUTINE delete_name

   SUBROUTINE search_name (name, found, stype)
     ! procedure searching a

     !Dummy arguments
     CHARACTER (len=*), INTENT(IN)  :: name
     INTEGER (kind=4), INTENT(IN OUT) :: stype
     LOGICAL, INTENT(OUT)         :: found

     ! local
     INTEGER (kind=4) :: i
     TYPE (set_label), POINTER :: addr

     found = .FALSE.          !initializes
     addr => snhead

     DO i=1, n_names
       IF ( TRIM(name) == TRIM(addr%sname) ) THEN   !label Found?
         IF( stype == 0 )THEN                !for any set type
           found = .TRUE.                    !O.K.
           stype = addr%stype                !note that different types
           EXIT                              !may have the same name
         ELSE IF( stype == addr%stype )THEN  !specific type
           found = .TRUE.                    !O.K.
           EXIT                              !
         END IF
       END IF
       addr => addr%next
     END DO

   END SUBROUTINE search_name

   SUBROUTINE dump_esets()

     !   dumps esets database

     IMPLICIT NONE

     ! local
     INTEGER (kind=4) :: i
     TYPE (set_label), POINTER :: addr

     WRITE (50,ERR=9999) n_names
     addr => snhead
     DO i=1,n_names
       WRITE(50,ERR=9999) addr%sname, addr%stype
       addr => addr%next
     END DO
     WRITE (50,ERR=9999) eset(1:msets), esets, nelms(1:msets), nel(1:msets)

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE dump_esets

   SUBROUTINE rest_esets()

     !   restores esets database

     IMPLICIT NONE
     ! local
     INTEGER (kind=4) :: i
     TYPE (set_label), POINTER :: addr,nadd

     READ (51) n_names

     ALLOCATE (addr)
     snhead => addr
     DO i=1,n_names
       IF(i < n_names)THEN
         ALLOCATE (nadd)
         addr%next => nadd
       ELSE
         sntail => addr
       END IF
       READ(51) addr%sname, addr%stype
       addr => addr%next
     END DO
     READ (51) eset(1:msets), esets, nelms(1:msets), nel(1:msets)
     rot_free = (nelms(11) + SUM(nelms(13:15)) + SUM(nelms(24:25))) > 0     !true if rotation free elements exist

     RETURN

   END SUBROUTINE rest_esets

 END MODULE esets_db
