 MODULE sms_db

   ! module to handle fibre information for Selected Mass Scaling

   USE param_db,ONLY: mnam

   IMPLICIT NONE
   LOGICAL :: selective_mass_scaling = .FALSE.  !may be in ctrl_db
   INTEGER(kind=4) :: sms_ns = 0,  & !number of element sets
                      sms_nf = 0,  & !number of fibres
                      sms_nnf = 2    !maximum number of nodes in the fibre

   CHARACTER(len=mnam), ALLOCATABLE :: sms_name(:) ! set names

   REAL(kind=8), ALLOCATABLE :: sms_alp(:),  & !(sms_ns) alpha factor for each set name
                                sms_thl(:),  & !(sms_ns) max layer thickness
                                sms_mf(:,:), & !(0:ndime,sms_nf) mass and factors
                                sms_faux(:,:)  !(ndime,sms_nf)   nodal forces and accelerations
   INTEGER(kind=4), ALLOCATABLE :: sms_fibre(:,:)  !(sms_nnf,sms_nf) nodes in each fibre

   ! types used to generate

   TYPE f_nod         !list of nodes
     INTEGER(kind=4) :: n
     TYPE (f_nod), POINTER :: next    !point to next
   END TYPE f_nod

   TYPE fiber         !fiber
     INTEGER(kind=4) :: nn   !number of nodes in the fiber
     REAL(kind=8) :: alp     !alpha factor
     TYPE (f_nod), POINTER :: head, tail    !first and last nodes
     TYPE (fiber), POINTER :: next          !pointer to next fiber
   END TYPE fiber

   TYPE(fiber), POINTER :: f_head,f_tail

 CONTAINS

 SUBROUTINE read_sms_inf( )
 ! read selective mass scaling     sms_nformation
 ! in MAINPG after damping isms_nformation

    USE c_input,ONLY: listen,exists,ludat,lures,backs,getint,getrea,get_name
    USE meshmo_db,ONLY: l_old
    USE npo_db, ONLY : coora,label
    USE ctrl_db, ONLY : npoin,ndime
    IMPLICIT NONE
    !Dummy argument
 
    !local variables
    INTEGER(kind=4), PARAMETER :: smax = 100,  & !maximum number of element sets
                                  mp = 9         !maximum number of nodes along a fibre + 1 (7 layers)
    CHARACTER(len=mnam) :: sname(smax) ! set names
    REAL(kind=8) :: alpha(smax),  & !(sms_ns) alpha factor for each set name
                    thlay(smax),  & !(sms_ns) max layer thickness
                    beta,y(ndime),lg(2),smin,largest
    LOGICAL :: found
    INTEGER(kind=4) :: i,j,k,l,n,ii,jj,ll,nnode,nelem,ip,iel,nar,max_nn,counter
    INTEGER(kind=4), ALLOCATABLE :: nod(:)
    REAL(kind=8), ALLOCATABLE :: x(:,:),alp(:)
    ! for this approach element connectivities have been ordered to minimize checks
    INTEGER (kind=4) :: arista(2,6,5) = (/ 1,2, 0,0, 0,0, 0,0, 1,2, 2,3, & !triangle
                                           2,3, 4,1, 0,0, 0,0, 1,2, 3,4, & !quad4
                                           1,2, 0,0, 0,0, 0,0, 1,3, 1,4, & !tetrahedra  (not programed yet)
                                           1,4, 2,5, 3,6, 0,0, 1,2, 1,3, & !prism
                                           1,5, 2,6, 3,7, 4,8, 1,2, 1,4  /) !hexahedra
    INTEGER (kind=4) :: naris(5) = (/ 1,2,1,3,4 /)   !number of aristas a each element type
    TYPE(f_nod), POINTER :: fn,n_new,fi
    TYPE(fiber), POINTER :: fib,f_new,fa

    INTERFACE
      INCLUDE 'elemnt.h'
    END INTERFACE

    CALL listen('SMSINF')                 !read first card

    IF (.NOT.exists('SELECT')) THEN       !key-word SELECT  not present

      backs = .TRUE.                       !one line back
      IF ( sms_ns > 0 ) &
        WRITE (lures,"('SMS Sets from the previous strategy used.')",ERR=9999)

    ELSE                                  !key-word SELECT exists
      ! data can not be added nor partially deleted, must be input completely
      sms_ns = 0                            !initializes
      IF( ALLOCATED(sms_name)) DEALLOCATE(sms_name,sms_alp,sms_thl,sms_mf,sms_faux,sms_fibre)
      max_nn = getint('MAXNN',  8,' Maximum number of nodes in the thickness ..')
      DO                                     ! loop over data sets
        CALL listen('SMSINF')                ! read a card
        IF (exists('ENDSEL')) EXIT           ! key-word END_SELECT exit loop
        sms_ns = sms_ns + 1                  ! update number of sets
        IF( sms_ns > smax) CALL runend('SMS_INF: number of sets exceeded')
        sname(sms_ns) = get_name('ELSNAM',found,stype='ESET')       !set name
        IF( .NOT.found ) CALL runend('SMS_INF: Element Set not found ')
        beta = getrea('BETA', 0d0,' Scaling factor for mass scaling...')
        IF( beta ==  0d0 )THEN
          alpha(sms_ns) = getrea('ALPHA',1d0,'!Inverse Scaling factor for mass scaling...')
        ELSE
          alpha(sms_ns) = 1d0/beta
        END IF
        thlay(sms_ns) = getrea('THLAY',0d0,' Maximum thickness layer for mass scaling...')
        IF( thlay(sms_ns) == 0d0 )THEN
          beta = getrea('THICN',0d0,'!lamina  thickness of the set ..............')
          i    = getint('NLAYR',0  ,'!Number of layer in the thickness ..........')
          thlay(sms_ns) = beta/REAL(i)*1.5d0
        END IF
      END DO                                 ! end SMS data
      IF( sms_ns > 0 )THEN                   ! if sets exist
        ALLOCATE( sms_name(sms_ns),sms_alp(sms_ns),sms_thl(sms_ns))  !get memory for basic data
        sms_name = sname(1:sms_ns)         ! keep names
        sms_alp  = alpha(1:sms_ns)         ! keep alpha values
        sms_thl  = thlay(1:sms_ns)         ! keep maximum distances
        alp = 1d0
        sms_nf = 0                         ! initializes number of fibres
        DO i=1,sms_ns                      ! for each set
          counter = 0                      ! initializes counter for this set
          largest = 0d0                    ! initializes largest arista for this set
          CALL elemnt('SLNODS',name=sms_name(i),lnod=l_old,flag2=found,flag1=.TRUE.) ! ==>  l_old(nnod,nelm)
          IF(.NOT.found) CALL runend('SMS_INF: Elm set not found')    !check existence
          IF(.NOT.ASSOCIATED(l_old))CALL runend('SMS_INF: elm_type not admissible')
          nnode = SIZE(l_old,1)                            !number of nodes per element
          ! fix position in array ARISTA
          IF( ndime == 2 )THEN                             !2D problems
            IF( nnode == 3) ip = 1                           !3-node triangle
            IF( nnode == 4) ip = 2                           !4-node quad
          ELSE                                             !3D problems
            IF( nnode == 4) ip = 3                           !4-node tetrahedra (NOT IN USE)
            IF( nnode == 6) ip = 4                           !6-node prism
            IF( nnode == 8) ip = 5                           !8-node hexahedra
          END IF
          nar = naris(ip)                                  !number of aristas per element
          nelem = SIZE(l_old,2)                            !number of elements
          ALLOCATE(x(ndime,nnode))                         !memory for element coordinates
          DO iel=1,nelem                                   !for each element
            x = coora(:,l_old(:,iel))                        !element coordinates
            DO j=1,2                                         ! compute the in-plane sides
              ii = arista(1,j+4,ip)                              !first node
              jj = arista(2,j+4,ip)                              !second node
              y(:) = x(:,jj) - x(:,ii)                           !vector
              lg(j) = DOT_PRODUCT(y,y)                           !distance(squared)
            END DO
            smin = MIN(lg(1),lg(2))                            !minimum in-plane length
            DO j=1,nar                                       !for each arista
              ii = arista(1,j,ip)                              !first node
              jj = arista(2,j,ip)                              !second node
              y(:) = x(:,jj) - x(:,ii)                         !vector
              beta = DOT_PRODUCT(y,y)                          !distance squared
              IF( ip /= 3 .AND. beta/smin > sms_alp(i) )THEN   ! not for tetrahedra
                WRITE(55,"('Warning element squared aspect ratio: ',e12.4,    &
                & ' larger than provided alpha (beta^-1) value: ',e12.4, /    &
                & ' alpha (beta^-1) changed to Minimum value' )")beta/smin,sms_alp(i)
                sms_alp(i) = beta/smin
              END IF
              beta = SQRT(beta)                                !distance
              IF( beta > sms_thl(i) )THEN
                IF( ip == 3 ) CYCLE      !too long (tetrahedra only) CYCLE
                counter = counter + 1                          !increase counter
                IF( beta > largest ) largest = beta            !save value
              END IF
              found = .FALSE.                                  !initializes
              ii = l_old(ii,iel)                               !global node
              jj = l_old(jj,iel)                               !global node
              fib => f_head                                    !point to first fiber
              outer : DO k=1,sms_nf                            !search in list of fibres
                fn => fib%head                                 !point to first node in the fiber
                DO l=1,fib%nn                                  !for each node in the fiber
                  n = fn%n                                       !fiber node
                  IF( n == ii .OR. n == jj)THEN                   !if one of the nodes
                    found = .TRUE.                                  !modify flag
                    ll = l                                          !keep position
                    fn => fn%next                                   !go to next node
                    EXIT outer                                      !exit loop
                  END IF
                  fn => fn%next                                   !point to next node
                END DO
                fib => fib%next                                !point to next fibre
              END DO outer
              IF( found )THEN                                  !fibre exists
                found = .FALSE.                                   !initializes flag
                IF( n == ii )THEN                                !if first node was found
                  DO l=ll+1,fib%nn                                 !look for second node
                    n = fn%n                                         !node
                    IF( n == jj)THEN
                      found  = .TRUE.
                      EXIT                                 !if second node already in the list
                    END IF
                    fn => fn%next                                   !go to next node
                  END DO
                  IF( .NOT.found )THEN                               !second node not in the list
                    ALLOCATE(n_new)
                    n_new%n = jj                                   !add second node
                    CALL add_node(fib,n_new)
                    CYCLE
                  END IF
                ELSE                                             !if second node was found
                  DO l=ll+1,fib%nn                                 !look for second node
                    n = fn%n                                         !node
                    IF( n == ii)THEN
                      found  = .TRUE.
                      EXIT                                 !if second node already in the list
                    END IF
                    fn => fn%next                                   !go to next node
                  END DO
                  IF( .NOT.found )THEN                               !second node not in the list
                    ALLOCATE(n_new)
                    n_new%n = ii                                   !add second node
                    CALL add_node(fib,n_new)
                    CYCLE
                  END IF
                END IF
                fib%alp = MIN(fib%alp,sms_alp(i))
                IF( fib%nn > max_nn )THEN !check
                  WRITE(lures,"(' Selective Mass Scaling ERROR',/, &
                       ' The maximum number of nodes to be colapsed MAXNN have been exceeded',i4,/,&
                       ' Increase if necessary or decrement THLAY, List follows')")max_nn
                  ALLOCATE( nod(max_nn+1) )
                  fn => fib%head
                  DO l=1,fib%nn
                    nod(l) = label(fn%n)
                    fn => fn%next
                  END DO
                  WRITE(lures,"( 10i7)")nod
                  CALL runen3(' SMS error, maximum number of nodes along a line exceeded')
                END IF
              ELSE                                             !fibre does not exist
                ALLOCATE(f_new)                            !allocate new fibre
                f_new%nn = 2                               !first two nodes
                ALLOCATE(n_new)                            !space for first node
                n_new%n = ii                               !pass information
                f_new%head => n_new                        !point to first
                NULLIFY(n_new)                             !
                ALLOCATE(n_new)                            !space for second node
                n_new%n = jj                               !pass information
                f_new%head%next => n_new                   !link with first node
                f_new%tail => n_new                        !point to last
                f_new%alp = sms_alp(i)
                CALL add_fibre(f_head,f_tail,f_new,sms_nf)
              END IF
            END DO
          END DO
          IF( counter > 0 )THEN
            WRITE(55,"('Warning ',i6,' aristas founded with length larger than provided (',e12.4,')', &
                        ' maximum layer thickness found: ',e12.4)")counter,sms_thl(i),largest
          END IF
          DEALLOCATE(x,l_old)                           !release memory
        END DO
        ! check and collapse fibres sharing a node
        fib => f_head      !point to first fiber
        DO                     !loop over each fibre
          fn => fib%head           !point to first node of the fiber
          k=1
          DO                       !loop for every node in the fiber
            found = .FALSE.          !initializes to no coincident nodes
            n = fn%n                 !fiber nodes
            fa    => fib             !point to actual fibre
            f_new => fib%next        !point to next fibre
            IF( .NOT.ASSOCIATED(f_new) )EXIT
            loop : DO                       !over each other fiber (below)
              n_new => f_new%head           !point to first node of other fibre
              DO l=1,f_new%nn               !por each node in other fibre
                ii = n_new%n                   !node
                IF( ii == n )THEN              !collapse fiber f_new
                  found = .TRUE.               !
                  EXIT loop                    !exit loop
                END IF
                n_new => n_new%next         !point to next node in the fiber
              END DO
              fa => f_new                ! keep previous fiber
              f_new => f_new%next        !point to next fibre
              IF( .NOT.ASSOCIATED(f_new) )EXIT loop !all fibers processed
            END DO loop
            IF( found )THEN                !if coincident nodes collapse second fiber into first
              n_new => f_new%head           !point to first node of fibre to delete
              DO i=1,f_new%nn               !por each node in other fibre
                n = n_new%n                 !internal node
                fi => fib%head              !point to first node of the fiber to be kept
                found = .FALSE.             !initializes search of node n
                DO j=1,fib%nn               !por each node in other fibre
                  jj = fi%n                 !internal node
                  IF( jj == n )THEN           !check
                    found = .TRUE.            !node exist
                    EXIT                      !exit search
                  END IF
                  fi => fi%next
                END DO
                IF( .NOT. found ) THEN     !search finished and the node is not on the fibre
                  ALLOCATE(fi)             !get memory
                  fi%n = n                 !pass data
                  CALL add_node(fib,fi)    !add to fiber
                END IF
                n_new => n_new%next
              END DO
              fa%next => f_new%next        !circunvect fiber
              DEALLOCATE (f_new)
              !f_new => fa%next             !point to next fiber
              sms_nf = sms_nf - 1
              CYCLE
            END IF
            k = k + 1             !counter on number of checked nodes
            IF( k > fib%nn )EXIT  !if number of nodes exceeded exit
            fn => fn%next         !else go to next nodes
            IF( .NOT.ASSOCIATED(fn)) EXIT  !it should not happen
          END DO
          fib => fib%next
          IF( .NOT.ASSOCIATED(fib) )EXIT
        END DO
        !find maximum number of nodes per fibre
        sms_nnf = 2                                            !initializes
        fib => f_head
        DO j=1,sms_nf                                          !for each fibre
          IF( fib%nn > sms_nnf ) sms_nnf = fib%nn
          fib => fib%next
        END DO
        WRITE(lures,"(' Number of fibres:', i5,/,' Max No of nodes:',i5)")sms_nf,sms_nnf
        IF( sms_nf > 0 )THEN
          ALLOCATE( sms_fibre(sms_nnf,sms_nf), sms_mf(0:ndime,sms_nf), sms_faux(ndime,sms_nf))
          sms_fibre = 0 !initializes
          fib => f_head          !point to first
          !for each fibre
          DO j=1,sms_nf
            sms_mf(0,j) = fib%alp  !keep alpha value
            fn => fib%head           !point to first node
            DO k=1,fib%nn            !for each node
              sms_fibre(k,j) = fn%n     !pass information
              fn => fn%next
            END DO
            fib => fib%next          !next fibre
          END DO
          selective_mass_scaling = .TRUE.
          CALL sms_erase()
        END IF
      ELSE
        selective_mass_scaling = .FALSE.
        IF( ALLOCATED(sms_name)) DEALLOCATE(sms_name,sms_alp,sms_thl,sms_mf,sms_faux,sms_fibre)
      END IF
    END IF

  RETURN
  9999 CALL runen2('')

 END SUBROUTINE read_sms_inf

 SUBROUTINE fibre_mass ( emass )
    ! compute fibre_mass
    ! called from lumass, slumas and explit
 USE kinc_db, ONLY : nn,nn1
 USE ctrl_db, ONLY : ndime,ndofn,npoin
 USE npo_db, ONLY : ymass,ifpre
 IMPLICIT NONE
 INTEGER(kind=4) :: i,j,k,l,ieq
 REAL(kind=8) :: emass(ndofn,npoin)

    IF( sms_nf == 0 )RETURN !no fibres
    sms_mf(1:ndime,:) = 0d0    !initializes
    DO i=1,sms_nf        !for each fibre
      DO j=1,sms_nnf        !for each node in the fibre
        k = sms_fibre(j,i)    !node number
        IF(k == 0)EXIT        !all nodes processed EXIT fiber
        DO l=1,ndime                 !for each space dimension
          ieq = ifpre(l,k)           !global DOF
          SELECT CASE (ieq)
          CASE (1: )                 !active DOF
            sms_mf(l,i) =  sms_mf(l,i) + ymass(ieq)   !DOF equivalent mass
          CASE(-nn:-1 )              !if a slave DOF
            sms_mf(l,i) =  sms_mf(l,i) + emass(l,k)  !add nodal mass
          CASE( :-nn1)              !if a prescribed DOF
            sms_mf(l,i) = 1e20        !large mass
          END SELECT
        END DO
      END DO
      sms_mf(:,i) = sms_mf(:,i) / (1d0-sms_mf(0,i))   !modify
    END DO
    RETURN
 END SUBROUTINE fibre_mass

 SUBROUTINE dump_sms_data (ndime)
   !called from DUMPIN
   IMPLICIT NONE
   INTEGER(kind=4), INTENT(IN) :: ndime
   ! local variables
   INTEGER(kind=4) :: i,j

   WRITE(50,ERR=9999) selective_mass_scaling,sms_ns,sms_nf,sms_nnf

   IF( .NOT. selective_mass_scaling ) RETURN

   WRITE(50,ERR=9999) sms_name, sms_alp,  sms_thl
   WRITE(50,ERR=9999) ((sms_fibre(i,j),i=1,sms_nnf),j=1,sms_nf)
   WRITE(50,ERR=9999) (sms_mf(0:ndime,j),j=1,sms_nf)

   RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump_sms_data

 SUBROUTINE rest_sms_data (ndime)
   !called from RESTAR
   IMPLICIT NONE
   ! dummy arguments
   INTEGER(kind=4), INTENT(IN) :: ndime
   ! local variables
   INTEGER(kind=4) :: i,j

   READ(51) selective_mass_scaling,sms_ns,sms_nf,sms_nnf

   IF( .NOT. selective_mass_scaling ) RETURN

   ALLOCATE (sms_name(sms_ns), sms_alp(sms_ns),  sms_thl(sms_ns), &
             sms_mf(0:ndime,sms_nf), sms_faux (ndime,sms_nf),     &
             sms_fibre(sms_nnf,sms_nf))

   READ(51) sms_name, sms_alp,  sms_thl
   READ(51) ((sms_fibre(i,j),i=1,sms_nnf),j=1,sms_nf)
   READ(51) (sms_mf(0:ndime,j),j=1,sms_nf)

   RETURN
 END SUBROUTINE rest_sms_data

  SUBROUTINE add_node(fib, node)
    ! insert a node at the end of the list

    !Dummy arguments
    TYPE (fiber), POINTER :: fib
    TYPE (f_nod), POINTER :: node


    fib%tail%next => node           !generate link
    fib%tail => node                !keep position of last
    fib%nn = fib%nn + 1             !increase node counter
    NULLIFY(fib%tail%next)          !point present to new

  END SUBROUTINE add_node

  SUBROUTINE add_fibre(head,tail,new, nf)
    ! insert a node at the end of the list

    !Dummy arguments
    TYPE (fiber), POINTER :: head,tail, new
    INTEGER (kind=4), INTENT(IN OUT) :: nf

    IF( ASSOCIATED (head) )THEN
      tail%next => new           !generate link
      tail => new                !keep position of last

    ELSE
      head => new
      tail => new
    END IF
    NULLIFY(tail%next)
    nf = nf + 1                   !increase number of fibres

  END SUBROUTINE add_fibre

  SUBROUTINE sms_erase
  ! release auxiliar memory used
  IMPLICIT NONE
    TYPE(f_nod), POINTER :: fn    !auxiliar pointer to node
    TYPE(fiber), POINTER :: fib   !auxiliar pointer to fiber
    INTEGER (kind=4) :: i,j       !auxiliar counters

    DO i=1,sms_nf                   !for
      fib => f_head                 !point to first fiber
      DO j=1,fib%nn                    !for every node
        fn => fib%head                 !point to node
        fib%head => fn%next            !update
        DEALLOCATE(fn)                 !release node memory
      END DO
      f_head => fib%next            !update
      DEALLOCATE(fib)               !release fiber memory
    END DO

  END SUBROUTINE sms_erase

  SUBROUTINE updlon_sms ()
  USE npo_db, ONLY : oldlb
  IMPLICIT NONE
  INTEGER(kind=4) :: i,j,k,lab,nn,nf,new(sms_nnf),chnode

    IF( sms_nf == 0 )RETURN !no fibres
    nf = 0           !initializes number of fibres
    DO i=1,sms_nf        !for each fibre
      nn = 0               !initializes nodes in this fibre
      new = 0               !intializes
      DO j=1,sms_nnf        !for each node in the fibre
        k = sms_fibre(j,i)  !node
        IF(k == 0)EXIT      !if node does not exist
        lab  = oldlb(k)     !original label
        lab  = chnode(lab)  !new internal number
        IF( lab > 0 )THEN      !if node exists yet
          nn = nn + 1          !updates number of nodes on the fiber
          new(nn) = lab        !pass new internal position
        END IF
      END DO
      IF( nn > 0 )THEN      !if there are nodes on the fiber
        nf = nf+1              !increase number of fibers
        sms_fibre(:,nf) = new  !update internal data base
      ELSE
        sms_fibre(:,i) = 0  !null values of this old fiber
      END IF
    END DO
    IF( nf == 0 )THEN       !if there are not active fibers
      DEALLOCATE(sms_alp,sms_thl,sms_mf,sms_faux,sms_fibre)  !release memory
      selective_mass_scaling = .FALSE.       !deactive sms
    ELSE
      sms_nf = nf           !keep number of fibers but do not modify array sizes
    END IF
    RETURN
  END SUBROUTINE updlon_sms

 END MODULE sms_db
