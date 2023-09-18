MODULE damp_db
  ! damping database + damping routines
  ! generated and read at RDDAMP
  ! regenerated for new strategies at RDDAMP and MRNDCR
  ! used at EXPLIT
  ! saved and re-read at DUMPIN RESTAR
  USE param_db,ONLY: mnam
  IMPLICIT NONE
  SAVE             !

  PRIVATE          !all private although explicitly stated

  PUBLIC ::   &   !variables
            dtype,      & !damping type, 1: standard viscous, 2:non-viscous
            ndamp,    & ! number of damping data sets
            damp,       & ! resultant damping
                 ! routines
            rddamp,     & ! reads damping DATA
            gres_damp,  & ! generate resultant damping
            updt_dval,  & ! updates automatic damping values (STATIC)
            dump_damp,  & ! dumps damping data for restart
            rest_damp,  & ! restores damping data at restart
            updlon_damp   ! updates damping data after reading/deleted nodes

  INTEGER (kind=4) :: dtype=1, &!damping type, 1:standard viscous, 2:non-viscous
                      ndamp=0   !number of damping data sets

  ! Derived type for a list of nodes with the same damping
  TYPE damp_nod
    INTEGER (kind=4) :: node            !node label
    TYPE (damp_nod), POINTER :: next    !point to next
  END TYPE damp_nod

  ! Derived type for a list of nodal sets with the same damping
  TYPE damp_nds
    CHARACTER (len=mnam) :: name        !nodal set name
    TYPE (damp_nds), POINTER :: next    !point to next
  END TYPE damp_nds

  ! Derived type for a set of damping data
  TYPE damp_set
    CHARACTER (len=mnam) :: name        ! damping set name
    REAL (kind=8) :: xitat(2)           ! damping coefficients (alpha, ? )
    INTEGER (kind=4) :: numpn           ! number of points in the set (excluding set points)
    INTEGER (kind=4) :: numst           ! number of node sets in the set
    TYPE (damp_nod), POINTER :: headn, tailn ! head and tail of the list of nodes
    TYPE (damp_nds), POINTER :: heads, tails ! head and tail of the list of sets
    !
    TYPE (damp_set), POINTER :: next    !pointer to next damping set
  END TYPE damp_set

  TYPE (damp_set), POINTER :: headd, taild  !first and last damping sets

  REAL (kind=8), POINTER :: damp(:)   !(neq) resultant damping coefficient
                                      !allocated if ndamp > 0
CONTAINS
                 !PUBLIC
  ! rddamp              reads damping DATA
  ! gres_damp           generate resultant damping (superposition of all the sets)
  ! dump_damp           dumps damping database for restart
  ! rest_damp           restores damping data at restart
  ! new_damp            updates damping data after remeshing
                 !PRIVATE
  ! rd_damps            reads damping data for one set
  ! alloc_damps         allocates and initializes a damping data set
  ! ini_dampn           initialize a list of nodes
  ! add_dampn           add node to a list of nodes
  ! ini_damps           initialize a list of damping data sets
  ! add_damps           adds data set to the end of the list of sets
  ! dalloc_damps        deallocates a damping data set
  ! del_all_damps       deletes all existing damping data sets releasing memory
  ! srch_dn             searches for an node labeled "node"
  ! del_dn              deletes a node pointed with posic

  SUBROUTINE rddamp (iwrit)
    !
    ! reads damping DATA
    !
    USE c_input,ONLY: listen,exists,ludat,lures,backs,getint
    IMPLICIT NONE
    ! dummy argument
    INTEGER(kind=4),INTENT(IN) :: iwrit   !echo flag
    ! local variables
    TYPE (damp_set), POINTER :: ds

    CALL listen('RDDAMP')                 !read first card
    IF (.NOT.exists('DAMPIN')) THEN       !key-word DAMPING not present
      backs = .TRUE.                        !one line back
      IF ( ndamp > 0 ) WRITE (lures,"('Damping from the previous strategy used.')",ERR=9999)
    ELSE                                  !key-word DAMPING exists
      dtype = getint('DTYPE ',1, ' Damping type (1:stand,2:nonvisco.')
      IF ( ndamp > 0 ) CALL del_all_damps    !delete previous damping definition
      CALL ini_damps (headd, taild)          !initializes pointers
      ndamp = 0                              !initializes
      DO                                     ! loop over damping data sets
        CALL listen('RDDAMP')                ! read a card
        IF (exists('ENDDAM')) EXIT           ! key-word END_DAMPING exit loop
        ndamp = ndamp + 1                    ! increase number of damping sets
        CALL alloc_damps (ds)                ! allocate a new set
        CALL rd_damps (ds, dtype, iwrit)     ! read a new set
        CALL add_damps (ds, headd, taild)    ! add a new set to the damping database
      END DO                                 ! end damping data
    END IF

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE rddamp

  SUBROUTINE rd_damps (ds, dtype, iwrit)
    !
    ! reads damping data for one set
    !
    USE param_db,ONLY: mnam
    USE nsets_db
    USE c_input,ONLY: get_name,listen,getrea,getint,exists,ludat,lures,backs,rdfrin,intrd,nintg,maxni
    !
    IMPLICIT NONE
    ! dummy arguments
    INTEGER(kind=4), INTENT(IN) :: dtype,  & ! damping type
                                   iwrit     ! echo flag
    TYPE (damp_set), POINTER :: ds           ! damping set

    ! local variables
    LOGICAL :: founds
    CHARACTER(len=mnam) :: sname
    INTEGER(kind=4) :: ii,ki,kf,kg,inode
    REAL (kind=8) :: a,b,freq
    TYPE (damp_nod), POINTER :: dn
    TYPE (damp_nds), POINTER :: nds
    TYPE (nset), POINTER :: ns

    ds%name = get_name('DAMPSE',founds,' DAMPING DATA SET: ',stype='NSET')!Damping set name

    CALL listen('RDDAMP')                                 !read a card
    IF (dtype == 1) THEN
      ! standard damping
      a     = getrea('AMPLI ',0d0, '!Amplitude of damped vibrations....')
      freq  = getrea('FREQ  ',0d0, '!Period ...........................')

      IF(a > 0  .AND. freq /= 0)THEN
        ds%xitat(1) = LOG(100/a)/ABS(freq)   !alpha value
        ds%xitat(2) = 1d0/freq               !inverse of period (frequency)
      ELSE
        ! write error or warning
        WRITE(lures,"(/,' DAMPING DATA ERROR, NOT considered: ',A)",ERR=9999) TRIM(ds%name)
        WRITE(55,"(/,' DAMPING DATA ERROR, NOT considered: ',A)",ERR=9999) TRIM(ds%name)
        ds%xitat(1) = 0d0
      END IF
    ELSE !IF (dtype == 2) THEN
      ! non-viscous damping
      a = getrea('ALPHA ',0d0,'!Damping constant for transl.dof...')
      IF (a < 0d0 .OR. a > 1d0) &
        CALL runend('Rd_damps: Damping constant should be > 0 and < 1')
      ds%xitat(1) = a
      b = getrea('BETA  ',a,' Damping constant for rotat. dof...')
      IF (b < 0d0 .OR. b > 1d0) &
        CALL runend('Rd_damps: Damping constant should be > 0 and < 1')
      ds%xitat(2) = b
    END IF
                                      !read associated nodes
    CALL listen('RDDAMP')             !read a card
    DO
      IF(exists('ENDDAM')) EXIT       !key-word END_DAMPING => exit loop

      IF (exists('FROM  ')) THEN      !List of nodes as an explicit loop
        ki = getint('FROM  ',0, '!From node  ......................')
        kf = getint('TO    ',0, '!To node .........................')
        kg = getint('STEP  ',1, ' Step (every) ....................')
        DO inode=ki,kf,kg             !loop
          ALLOCATE (dn)               !get memory for a new node
          dn%node = inode             !store node label
          ds%numpn = ds%numpn + 1     !increase number of nodes in the list
          CALL add_dampn (dn, ds%headn, ds%tailn) !add node to the list
        END DO
        IF(iwrit > 0)WRITE(lures,"(8i10)",ERR=9999) (inode,inode=ki,kf,kg) !echo

      ELSE IF(exists('SET   ',ii))THEN !list of nodes in a set
        sname = get_name(posin=ii,stype='NSET') !set name
        CALL nsdb_search (sname, founds, ns)    !see if nodes set exists
        IF (founds) THEN                        !if set found
          ALLOCATE (nds)               !get memory for a new nodal set
          nds%name = sname
          CALL add_dampns (nds, ds%heads, ds%tails) !add node to the list
          ds%numst = ds%numst + 1     !increase number of nodes in the list
          IF(iwrit > 0)WRITE(lures,"(' SET:',a)",ERR=9999) sname
        ELSE         !set not found => ERROR
          WRITE (lures, '(" Set ",A,"  not found")',ERR=9999) TRIM(sname)
          CALL runend('Intime: Set not found ')
        END IF

      ELSE      ! read a list of nodes

        backs = .TRUE.                       !one line back
        CALL rdfrin('RDDAMP',intrd,nintg,maxni)  !read a list of labels
        DO kg=1,nintg                            !loop for the nodes in the list
          inode = intrd(kg)            !node label
          ALLOCATE (dn)                !get memory
          dn%node = inode              !store node label
          ds%numpn = ds%numpn + 1      !increase number of node in the set
          CALL add_dampn (dn, ds%headn, ds%tailn) !add node to the list
        END DO
        IF(iwrit > 0)WRITE(lures,"(8i10)",ERR=9999) (intrd(ii),ii=1,nintg) !ECHO

      END IF

      CALL listen('RDDAMP')  !read a new card
    END DO

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE rd_damps

  SUBROUTINE gres_damp (ndime, ndofn, neq, ifpre)
    !
    ! Generate RESultant DAMPing (superposition of all the sets)
    !
    USE nsets_db
    IMPLICIT NONE
    ! dummy arguments
    INTEGER (kind=4), INTENT(IN) :: ndime, ndofn, neq, ifpre(:,:)
    ! local variables
    INTEGER (kind=4) :: n
    TYPE (damp_set), POINTER :: ds
    TYPE (damp_nod), POINTER :: dn
    TYPE (damp_nds), POINTER :: nds
    LOGICAL :: found
    TYPE(nset), POINTER :: ns

    IF (ASSOCIATED(damp)) DEALLOCATE( damp)  !release previous memory
    ALLOCATE( damp(neq) )                    !allocate new memory
    damp = 0d0                               !initializes

    ds => headd                              !point to first damping set
    ! loop over damping data sets
    DO
      IF (.NOT.ASSOCIATED (ds) ) EXIT        !last set reached, exit loop
      ! go through the nodes from the set
      dn => ds%headn                         !point to first node in the set
      DO
        IF (.NOT.ASSOCIATED (dn) ) EXIT      !last node reached, exit loop
        CALL add_node_to_array(dn%node,ndime,ndofn,ifpre,ds%xitat)
        dn => dn%next                        !point to next node in the list
      END DO
      ! go through the sets
      nds => ds%heads                       !point to first nodal set
      DO
        IF (.NOT.ASSOCIATED (nds) ) EXIT    !last nodal set reached, exit loop
        CALL nsdb_search (nds%name, found, ns)    !point to nodal set
        CALL ns_head (ns)                     !point to first node in the set
        DO                                    !loop over all nodes in the set
          IF ( end_ns(ns) ) EXIT              !if last node found Exit loop
          n = get_label(ns)                   !get node label
          CALL add_node_to_array(n,ndime,ndofn,ifpre,ds%xitat)
          CALL ns_next (ns)                !point to next node in the list
        END DO
        nds => nds%next                        !point to next node in the list
      END DO
      ds => ds%next                          !point to next set in the list
    END DO
    RETURN
  END SUBROUTINE gres_damp

  SUBROUTINE add_node_to_array(node,ndime,ndofn,ifpre,xitat)
    IMPLICIT NONE
    ! arguments
    INTEGER(kind=4) :: node, ndime, ndofn, ifpre(:,:)
    REAL(kind=8) :: xitat(2)
    ! local
    INTEGER (kind=4) :: n, idof, ieq

    ! external
    INTEGER (kind=4) :: chnode

    n = chnode (node)                 !internal node number
    DO idof = 1,ndime                      !for each transl. DOF
      ieq = ifpre(idof,n)                  !associated equation
      ! note that next line accumulates damping effects (FF)
      IF (ieq > 0) damp(ieq) = damp(ieq) + xitat(1)  !if active equation
    END DO

    DO idof = ndime+1,ndofn                !for each rotat. DOF
      ieq = ifpre(idof,n)                  !associated equation
      ! note that next line accumulates damping effects (FF)
      IF (ieq > 0) THEN                    !if active equation
        IF (dtype == 1) THEN  ! standard viscous damping
          damp(ieq) = damp(ieq) + xitat(1)
        ELSE !IF (dtype==2): non-viscous damping can be different for rot.
          damp(ieq) = damp(ieq) + xitat(2)
        END IF
      END IF
    END DO
  END SUBROUTINE add_node_to_array


  SUBROUTINE dump_damp (neq)

    ! dumps damping database for restart
    IMPLICIT NONE

    ! arguments
    INTEGER (kind=4) :: neq

    ! local
    INTEGER (kind=4) :: i
    TYPE (damp_set), POINTER :: ds
    TYPE (damp_nod), POINTER :: dn
    TYPE (damp_nds), POINTER :: nds

    WRITE (50,ERR=9999) dtype,ndamp
    IF (ndamp == 0) RETURN
    ! damping vector could be regenerated at restart (FF)
    WRITE (50,ERR=9999) (damp(i), i=1,neq)  !damping vector

    ds => headd                  !point to first damping set
    ! loop over damping data sets
    DO
      IF (.NOT.ASSOCIATED (ds) ) EXIT         !last set processed, Exit loop
      WRITE (50,ERR=9999) ds%name, ds%xitat, ds%numpn, ds%numst  !name, damping factor, No of nodes
      ! go through the nodes from the set
      dn => ds%headn      !point to first node in the set
      DO
        IF (.NOT.ASSOCIATED (dn) ) EXIT       !last node processed, exit loop
        WRITE (50,ERR=9999) dn%node                    !node label
        dn => dn%next                         !point to next node
      END DO
      nds => ds%heads      !point to first node in the set
      DO
        IF (.NOT.ASSOCIATED (nds) ) EXIT       !last node processed, exit loop
        WRITE (50,ERR=9999) nds%name                    !node label
        nds => nds%next                         !point to next node
      END DO
      ds => ds%next                         !point to next set
    END DO

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_damp

  SUBROUTINE rest_damp (neq)

    ! restores damping data at restart

    IMPLICIT NONE
    ! arguments
    INTEGER (kind=4) :: neq

    ! local
    INTEGER (kind=4) :: i, j
    TYPE (damp_set), POINTER :: ds
    TYPE (damp_nod), POINTER :: dn
    TYPE (damp_nds), POINTER :: nds

    READ (51) dtype,ndamp
    IF (ndamp == 0) RETURN

    ALLOCATE( damp(neq) )           !allocate new memory
    READ (51) (damp(i), i=1,neq)    !damping vector

    CALL ini_damps (headd, taild)   !initializes list of sets

    ! loop over damping data sets
    DO i = 1,ndamp                  !for each damping set
      CALL alloc_damps (ds)                  !get memory for set
      READ (51) ds%name, ds%xitat, ds%numpn, ds%numst  !name, damping factor, No of nodes
      ! go through the nodes from the set
      DO j = 1, ds%numpn                     !for each node in the list
        ALLOCATE (dn)                        !get memory for node
        READ (51) dn%node                    !node label
        CALL add_dampn (dn, ds%headn, ds%tailn)   !add node to the list
      END DO
      DO j = 1, ds%numst                     !for each set in the list
        ALLOCATE (nds)                      !get memory for set
        READ (51) nds%name                  !set name
        CALL add_dampns (nds, ds%heads, ds%tails)   !add set to the list
      END DO
      CALL add_damps (ds, headd, taild)      !add set to the list
    END DO
    RETURN
  END SUBROUTINE rest_damp

  SUBROUTINE alloc_damps (damps)
    !allocates and initializes a damping data set
    TYPE (damp_set), POINTER :: damps

    ALLOCATE (damps)        !get memory
    damps%name  = ''        !initializes name
    damps%numpn = 0         !initializes number of nodes
    damps%numst = 0         !initializes number of nodal sets
    damps%xitat = 0d0       !initializes damping parameter
    CALL ini_dampn (damps%headn, damps%tailn, damps%heads, damps%tails)  !initializes pointers
    RETURN

  END SUBROUTINE alloc_damps

  SUBROUTINE ini_dampn (head, tail, heads, tails)
    !initialize a list of nodes

    !Dummy arguments
    TYPE (damp_nod), POINTER :: head, tail
    TYPE (damp_nds), POINTER :: heads, tails

    NULLIFY (head, tail)
    NULLIFY (heads, tails)

  END SUBROUTINE ini_dampn

  SUBROUTINE add_dampn (new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (damp_nod), POINTER :: new, head, tail

    IF (.NOT. ASSOCIATED (head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    ENDIF
  END SUBROUTINE add_dampn

  SUBROUTINE add_dampns (new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (damp_nds), POINTER :: new, head, tail

    IF (.NOT. ASSOCIATED (head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    ENDIF
  END SUBROUTINE add_dampns

  SUBROUTINE ini_damps (head, tail)
    !initialize a list of damping data sets

    !Dummy arguments
    TYPE (damp_set), POINTER :: head, tail

    NULLIFY (head, tail)

  END SUBROUTINE ini_damps

  SUBROUTINE add_damps (new, head, tail)
    !adds data set to the end of the list of sets
    !Dummy arguments
    TYPE (damp_set), POINTER :: new, head, tail

    IF (.NOT. ASSOCIATED (head)) THEN    !Check if a list is empty
      head => new             !list is empty, start it
      tail => new
      NULLIFY (tail%next)

    ELSE                      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    END IF
    RETURN
  END SUBROUTINE add_damps

  SUBROUTINE dalloc_damps (damps)
    ! deallocates a damping data set
    ! dummy arguments
    TYPE (damp_set), POINTER :: damps
    !local
    TYPE (damp_nod), POINTER :: n, naux
    TYPE (damp_nds), POINTER :: nds, nsaux

    n => damps%headn       !point to first node in the set
    DO
      IF (.NOT.ASSOCIATED (n) ) EXIT   !if last node deallocated , exit loop
      naux => n%next                   !keep next position
      DEALLOCATE (n)                   !release memory
      n => naux                        !point to next node
    END DO
    nds => damps%heads       !point to first set
    DO
      IF (.NOT.ASSOCIATED (nds) ) EXIT   !if last set deallocated , exit loop
      nsaux => nds%next                   !keep next position
      DEALLOCATE (nds)                    !release memory
      nds => nsaux                        !point to next set
    END DO

    DEALLOCATE (damps)                 !release set memory
    RETURN
  END SUBROUTINE dalloc_damps

  SUBROUTINE del_all_damps
    ! deletes all existing damping data sets releasing memory
    !local
    TYPE (damp_set), POINTER :: ds, dsaux

    ds => headd         !point to first set
    DO
      IF (.NOT.ASSOCIATED (ds) ) EXIT   !if last set deleted, exit loop
      dsaux => ds%next                  !keep pointer to next set
      CALL dalloc_damps (ds)            !delete set
      ds => dsaux                       !point to next set
    END DO
    RETURN
  END SUBROUTINE del_all_damps

  SUBROUTINE srch_dn (head, anter, posic, node, found)
    ! searches for an node labeled "node"
    !Dummy arguments
    LOGICAL :: found
    INTEGER (kind=4) :: node
    TYPE (damp_nod), POINTER :: head, anter, posic

    found = .FALSE.
    NULLIFY (posic,anter)
    !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      posic => head
      DO
        IF(posic%node == node) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next) ) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    ENDIF
    IF (.NOT.found) NULLIFY (posic,anter)
  END SUBROUTINE srch_dn


  SUBROUTINE del_dn (head, tail, anter, posic)
    ! deletes a node pointed with posic
    TYPE (damp_nod), POINTER :: head, tail, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    ! if posic == tail
    IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
    DEALLOCATE (posic)
    NULLIFY (anter)
  END SUBROUTINE del_dn

  SUBROUTINE updlon_damp ( )

    ! updates list of nodes

    USE nsets_db
    IMPLICIT NONE

    ! local
    INTEGER (kind=4) ::  j, lab,chnode,n
    TYPE (damp_set), POINTER :: ds,prev       !damping set
    TYPE (damp_nod), POINTER :: dn,da,anter   !damping nodes
    TYPE (damp_nds), POINTER :: nds,aset,bset !damping node set
    TYPE (nset), POINTER :: ns                !nodal sets
    LOGICAL :: found

    IF (ndamp == 0) RETURN                    !if no damping sets

    ds => headd   !point to first damping set
    NULLIFY( prev ) !no previous set
    ! loop over damping data sets
    DO
      IF( .NOT.ASSOCIATED(ds) )EXIT
      ! go through the nodes from the set
      dn => ds%headn    !point to first node in the set
      NULLIFY(anter)    !no previous node
      n = 0             !initializes counter
      DO j = 1, ds%numpn                   !for each individual node in the list
        lab = chnode( dn%node )            !internal node
        IF( lab /= 0 )THEN                   ! if node yet exists
          anter => dn                        ! keep previous node
          da => anter                        ! previous node
          dn => anter%next                   ! point to next  position
          n = n+1                            ! increase number of actual node
        ELSE                               !if node deleted
          IF( n == 0 )THEN                   ! if no node in the list
            ds%headn => dn%next              ! change first node to present
            DEALLOCATE(dn)                   ! deallocate memory
            dn => ds%headn                   ! point again to first node
          ELSE
            anter%next => dn%next            !update pointer
            DEALLOCATE(dn)                   !release memory
            dn => anter%next                 !update pointer
          END IF
        END IF
      END DO
      ds%numpn = n                   !keep number of nodes

      ! go through the nodal sets in the set
      nds => ds%heads                !first set of nodes in damping set
      NULLIFY(bset)                          !nullify pointer to previous set
      n = 0                                  !initializes number of node sets
      DO j = 1, ds%numst                     !for each set in the list
        CALL nsdb_search (nds%name, found, ns)    !see if node set exists
        IF( found  )THEN                     !set exist
          bset => nds                        !??? keep pointer to previous set
          aset => bset                       !???
          nds => aset%next                   !???
          n = n+1                            !increase number of sets
        ELSE
          IF( n == 0 )THEN
            ds%heads => nds%next
            DEALLOCATE(nds)
            nds => ds%heads
          ELSE
            aset%next => nds%next
            DEALLOCATE(nds)
            nds => aset%next
          END IF
        END IF
      END DO
      ds%numst = n    !keep number of node sets in present damping set
      IF( ds%numst == 0 .AND. ds%numpn == 0 )THEN !if no nodes and sets delete
        !delete set
        IF( ASSOCIATED(prev) )THEN
          prev%next => ds%next
        ELSE
          headd => ds%next
        END IF
        CALL  dalloc_damps (ds)
        ndamp = ndamp - 1
        IF( ASSOCIATED( prev) )THEN
          ds => prev%next   !point to first damping set
        ELSE
          ds => headd
        END IF
      ELSE
        prev => ds
        ds => ds%next   !point to first damping set
        taild => prev
      END IF

    END DO
    RETURN
  END SUBROUTINE updlon_damp

  SUBROUTINE updt_dval (ndime,pmin,pmax,flag,fdamp,resid,resid_o,incdis,smass,itera)
    !
    ! Update damping factors automatically
    ! This routine is inefficient if called continuosly
    !
    USE nsets_db
    IMPLICIT NONE
    ! dummy arguments
    INTEGER (kind=4), INTENT(IN) :: ndime,pmin,pmax
    LOGICAL, INTENT(IN OUT) :: flag
    INTEGER (kind=4), OPTIONAL, INTENT(IN) :: itera
    REAL (kind=8), INTENT(IN), OPTIONAL :: resid(:,:),resid_o(:,:),incdis(:,:),smass(:,:),fdamp
    ! local variables
    INTEGER (kind=4) :: n,chnode
    REAL (kind=8) :: we,ki,alpha,alpha_o,ratio
    TYPE (damp_set), POINTER :: ds
    TYPE (damp_nod), POINTER :: dn
    TYPE (damp_nds), POINTER :: nds
    LOGICAL :: found
    TYPE(nset), POINTER :: ns

    ! automatic damping computations
    flag = .FALSE.                           !no modifications
    ds => headd                              !point to first damping set
    ! loop over damping data sets
    DO
      IF (.NOT.ASSOCIATED (ds) ) EXIT        !last set reached, exit loop
      ! go through the nodes from the set
      IF( ds%xitat(2) > 0d0 )THEN
        alpha_o = ds%xitat(2)  !old value (inverse of numbers of steps)
        IF( PRESENT(itera) )THEN             !use number of iterations (kinetic damping)
          alpha = 1d0/REAL(4*itera)
        ELSE
          we = 0d0          !initializes incremental work
          ki = 0d0          !initializes incremental kinetic energy
          dn => ds%headn                         !point to first node in the set
          DO
            IF (.NOT.ASSOCIATED (dn) ) EXIT      !last node reached, exit loop
            n = chnode(dn%node)
            we = we + DOT_PRODUCT(resid(1:ndime,n) - resid_o(:,n),incdis(:,n))
            ki = ki + smass(1,n)*DOT_PRODUCT(incdis(:,n),incdis(:,n))
            dn => dn%next                        !point to next node in the list
          END DO
          ! go through the sets
          nds => ds%heads                       !point to first nodal set
          DO
            IF (.NOT.ASSOCIATED (nds) ) EXIT    !last nodal set reached, exit loop
            CALL nsdb_search (nds%name, found, ns)    !point to nodal set
            CALL ns_head (ns)                     !point to first node in the set
            DO                                    !loop over all nodes in the set
              IF ( end_ns(ns) ) EXIT              !if last node found Exit loop
              n = get_label(ns)                   !get node label
              n = chnode(n)
              we = we + DOT_PRODUCT(resid(1:ndime,n) - resid_o(:,n),incdis(:,n))
              ki = ki + smass(1,n)*DOT_PRODUCT(incdis(:,n),incdis(:,n))
              CALL ns_next (ns)                !point to next node in the list
            END DO
            nds => nds%next                        !point to next node in the list
          END DO
          alpha = SQRT(ABS(we)/ki)   !new period inverse
        END IF
        IF( 1d0/alpha > pmax ) alpha = 1d0/pmax      !check not exceeding max
        IF( 1d0/alpha < pmin ) alpha = 1d0/pmin      !check not exceeding min
        ratio = alpha/alpha_o      !ratio
        IF( PRESENT(fdamp) )THEN
          ratio = MIN(ratio,fdamp)       !compare with upper limit
          ratio = MAX(ratio,1d0/fdamp)   !compare with lower limit
          WRITE(58,"(2i10,f7.2)")INT(1d0/alpha_o),INT(1d0/alpha),ratio
        END IF
        ds%xitat(1) = ds%xitat(1)*ratio  !recompute alpha
        ds%xitat(2) = ratio*alpha_o      !recompute alpha
        flag = .TRUE.
      END IF
      ds => ds%next                          !point to next set in the list
    END DO
    RETURN
  END SUBROUTINE updt_dval

END MODULE damp_db
