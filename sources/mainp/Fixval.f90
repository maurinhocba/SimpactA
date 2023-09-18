SUBROUTINE fixval(iwrit,ndime,ndofn,ifpre,nvfix,rot_free,iffix,label)

  !***  APPLY fixities
  USE ctrl_db, ONLY : nrotd
  USE c_input
  USE ifx_db
  USE kinc_db, ONLY : nn,nn1
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: rot_free
  INTEGER (kind=4), INTENT(IN) :: iwrit,ndime,ndofn,label(:)
  INTEGER (kind=4), INTENT(IN OUT) :: ifpre(:,:),iffix(:)
  INTEGER (kind=4), INTENT(OUT) :: nvfix

  INTEGER (kind=4) :: ifix(7),n,i,ipoin,chnode,np,kfix(7)
  TYPE (ifx_nod), POINTER :: ifx

  nvfix = 0  !initializes number of prescribed values

  IF(iwrit == 1) THEN
    IF(ndime == 2) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node   XYA')",ERR=9999)
    IF(ndime == 3) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node   XYZABG')",ERR=9999)
  END IF

  ifx => ihead                !point to first node in the list
  DO n=1,nifx                !for each node in the list

    np  = ifx%ifix(1)        !node label
    ipoin = chnode(np)       !node internal number
    IF( ipoin == 0 )THEN  !verify node exist
      ifx => ifx%next                     !point to next node
      CYCLE
    END IF
    ifix(1:nd1) = ifx%ifix(2:nd1+1)     !restriction codes
    kfix = ifix                         !keep original codes
    DO i=1,nrotd                   !for each DOF
      SELECT CASE (ifpre(i,ipoin))   !according to previous codes
      CASE (0)                         ! Active DOF
        IF(ifix(i) == 1) THEN            !if a restriction is included
          nvfix = nvfix+1                !increase number of fixed values
          ifpre(i,ipoin) = -(nn+nvfix)   !assign a position
        END IF
      CASE (1)                         ! Not an active DOF
        IF(ifix(i) == 0) THEN            !if a release code is included
          IF( i <= ndime .OR. ndofn <= 6 )&
          WRITE(*,"(' WARNING, To release Node ',i5,' DOF ',i2, &
                      & ' use Code 2 instead')",ERR=9999) label(ipoin),i
        ELSE IF(ifix(i) == 2) THEN            !if a release code is included
          IF (iwrit == 1) WRITE(lures,"(' WARNING, Node ',i5,' DOF ',i2, &
                      & ' Has been released')",ERR=9999) label(ipoin),i
          ifpre(i,ipoin) = 0             !release DOF
        END IF
      CASE (-nn:-1)                    !for slave DOF
        kfix(i) = 2                      !Modify KFIX for reference
      CASE (:-nn1)                     !for a prescribed DOF
        IF(ifix(i) == 0) THEN            !if a release code included, ERROR
          WRITE(lures,"(' ERROR, Node ',i5,' DOF ',i2, &
                      & ' was previously constrained')",ERR=9999) label(ipoin),i
          CALL runend('FIXVAL: Inconsistent input data    ')
        ELSE IF(ifix(i) == 1) THEN       !constrained twice, WARNING
          WRITE(lures,"(' WARNING, Node ',i5,' DOF ',i2, &
                      & ' was previously constrained')",ERR=9999) label(ipoin),i
        END IF
      END SELECT
    END DO
    IF(rot_free)kfix(nd1) = ifix(nd1)       !BST element type
    IF(iwrit == 1) WRITE(lures,"(i7,3x,7i1)",ERR=9999) label(ipoin),kfix(1:nd1)
    IF(rot_free)iffix(ipoin) = ifix(nd1)    !keep BST type code
    IF( ndofn == 8 )THEN      !for Psi functions
      IF ( ifpre(7,ipoin) == 0 )THEN     !if the DOF exist
        IF(ifix(5) == 1) THEN          !if a restriction is included
          nvfix = nvfix+1                !increase number of fixed values
          ifpre(7,ipoin) = -(nn+nvfix)   !assign a position
        END IF
      END IF
      IF ( ifpre(8,ipoin) == 0 )THEN     !if the DOF exist
        IF(ifix(4) == 1) THEN          !if a restriction is included
          nvfix = nvfix+1                !increase number of fixed values
          ifpre(8,ipoin) = -(nn+nvfix)   !assign a position
        END IF
      END IF
    END IF
    ifx => ifx%next                     !point to next node
  END DO

RETURN
 9999 CALL runen2('')
END SUBROUTINE fixval
