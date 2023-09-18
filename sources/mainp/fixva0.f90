SUBROUTINE fixva0 (iwrit,ndime,nrotd,actio,rot_free,nvelr)

  !***  READ the fixed values.

  USE param_db,ONLY: mnam
  USE ifx_db
  USE c_input
  USE nsets_db
  USE rve_db

  IMPLICIT NONE
  LOGICAL :: rot_free
  CHARACTER(len=*),INTENT(INOUT):: actio
  INTEGER (kind=4) :: iwrit,ndime,nrotd,nvelr

  ! Local
  LOGICAL set, found
  CHARACTER(len=mnam) :: sname
  TYPE (ifx_nod), POINTER :: ifx, anter, posic
  INTEGER :: i,j,k,n,ki,kf,icod,nnods,nod,nrve
  INTEGER (kind=4), ALLOCATABLE :: nods(:)
  REAL (kind=8) :: v
  TYPE (rve_set), POINTER :: rves
  TYPE (rve_nod), POINTER :: rven
  TYPE (nset), POINTER :: ns

  !------------------------------------------------------------

  IF (iwrit == 1) WRITE(lures,"(/,' B O U N D A R Y  C O N D I T I O N S',/)",ERR=9999)

  !Initialize empty list
  nd1 = nrotd
  IF (rot_free) nd1 = nd1 + 1
  IF (nifx == 0 .OR. .NOT.exists('ADD   ')) THEN
    CALL ini_ifx(ihead,itail)
    nifx=0                       !initializes for new boundary conditions
  ELSE IF (exists('ADD   ')) THEN
    IF (iwrit == 1) WRITE(lures,"(/, &
      & ' Kept Boundary Conditions from the previous strategy.',/)",ERR=9999)
  END IF

  IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'

  IF (iwrit == 1) THEN
    IF(ndime == 2) WRITE(lures,"(/,'     Node   XYA')",ERR=9999)
    IF(ndime == 3) WRITE(lures,"(/,'     Node   XYZABG')",ERR=9999)
  END IF

  !Loop to read data and add them to the list
  DO
    CALL listen('FIXVA0')          !read a line
    IF (exists('ENDBOU')) EXIT     !end of data =>exit

    IF (exists('SET   ',n))THEN                    !key-word SET found
      sname = get_name(posin=n,stype='NSET')       !set name
      CALL nsdb_search (sname, found, ns)          !search in list of sets
      IF (found) THEN                      !set exists, position is NS
        nnods = get_length (ns)            !number of nodes in the set
        ALLOCATE (nods(nnods))             !get memory
        CALL ns_head (ns)                  !go to top of list
        DO i =1, nnods                     !for each node in the list
          nods(i) = get_label(ns)          !get node label
          CALL ns_next (ns)                !go to next node
        END DO
      ELSE                                 !else SET name not found => error
        WRITE (lures, '(" Set ",A,"  not found")',ERR=9999) TRIM(sname)
        CALL runend('Intime: Set not found ')
      END IF
      ki = 1                               !first node
      kf = nnods                           !last node
      set = .TRUE.                         !set flag to TRUE
    ELSE                             !read fixities for a single node
      ki  = INT(param(1))                  !node label
      kf = ki                              !same as last
      set = .FALSE.                        !set flat to FALSE
    END IF
                                     !store fixities in list
    DO n = ki,kf                           !for each node in the list
      ALLOCATE (ifx)                       !get memory
      IF (set) THEN                        !for a whole set
        nod = nods(n)                      !node label
      ELSE
        nod = n                            !only node
      END IF
      ifx%ifix(1) = nod                    !store node label
      ifx%ifix(2:nd1+1)= INT (param(2:nd1+1))   !store fixities
      IF (iwrit == 1) WRITE(lures, '(i9,3x,7i1)',ERR=9999) ifx%ifix(1:nd1+1)  !echo

      !    search in all the prescribed values to change conditions

      rves => headv      !point to the first set of prescribed values

      DO i = 1,nvelr     !for each prescribed velocities set
        ! the convention: the last read overrides the previous data
        ! if prescribed velocity has been read previously it will
        ! be cancelled at fixed or released degrees of freedom
        nrve = rves%nrv       !number of prescribed values in the set
        rven => rves%head     !point to first node in the set
        DO j=1,nrve                        !for each node in the set
          IF (nod == rven%node) THEN       !compare labels
            DO k = 1,nrotd                    !for each DOF
              v = rven%v(k)                   !prescribed velocity in that DOF
              IF (v /= 0d0)THEN               !if original velocity /= 0
                icod= INT (param(k+1))        !fixity code
                rven%v(k) = 0d0               !initializes to 0
                IF (icod == 0) THEN           !if DOF is free now
                  WRITE (lures,"(' Warning. Node:',i8,', DOF:',i2,' released.'/&
                     & ' Previously prescribed velocity cancelled.')",ERR=9999) nod,k
                ELSE !icod == 1 (prescribed value)
                  WRITE (lures,"(' Warning. Node:',i8,', DOF:',i2,' fixed.'/ &
                     & ' Previously prescribed velocity cancelled.')",ERR=9999) nod,k
                END IF
              END IF
            END DO
            EXIT      !node label was found => exit this list of nodes
          END IF
          rven => rven%next   !point to next node in the list
        END DO        !list of nodes of a set of prescribed values
        rves => rves%next     !point to next list of nodes
      END DO   ! list of sets

      ! before adding check if the node has not been fixed previously
      CALL srch_ifx (ihead, anter, posic, ifx%ifix(1), found)
      IF (found) THEN
        WRITE (lures, "(' Warning. At node',i8,' previously prescribed ', &
                      & 'boundary conditions overwritten.')",ERR=9999)  ifx%ifix(1)
        CALL del_ifx (ihead, itail, anter, posic)
      ELSE
        nifx=nifx+1       !increase counter of fixed nodes
      END IF
      CALL add_ifx( ifx, ihead, itail )    !add node to the list
    END DO   !list of nodes (if set)
    IF (set) DEALLOCATE (nods)  !deallocate temporary array
  END DO  !fixed nodes, end of data of boundaries found

RETURN
 9999 CALL runen2('')
END SUBROUTINE fixva0
