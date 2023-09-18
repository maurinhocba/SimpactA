SUBROUTINE updl10 (set, oldlb)

  !***  updates internal node numbering

  IMPLICIT NONE
  INTEGER (kind=4) :: oldlb(:)
  TYPE(ele10_set) :: set

  INTEGER (kind=4) ielem,lab,n,chnode,i

  ! change label of connectivities
  IF( set%ntype == 0 )THEN   !particles

    ielem = 0
    DO n=1,set%nnode
      lab  = oldlb(set%lnods(n,ielem))
      lab  = chnode(lab)
      IF( lab > 0 )THEN
        ielem = ielem + 1
        set%lnods(ielem,1) = lab
      END IF
    END DO
    set%nnode = ielem

  ELSE                       !meshes

    ielem = 0
    outer : DO i=1,set%nelem
      DO n=1,set%nnode
        lab = oldlb(set%lnods(n,i))
        IF( lab == 0 ) EXIT
        lab = chnode(lab)
        IF( lab > 0 )THEN
          set%lnods(n,i) = lab
        ELSE
          set%lnods(:,i) = 0
          CYCLE outer
        END IF
      END DO
      ielem = ielem + 1
    END DO outer
    IF( ielem < set%nelem )THEN
      n = 0
      DO i=1,set%nelem
        IF(set%lnods(1,i) /= 0)THEN
          n = n + 1
          set%lnods(:,n) = set%lnods(:,i)
        END IF
      END DO
      set%nelem = ielem
    END IF
  END IF

  ! change label of rigid body node
  IF(set%rbnod > 0) THEN
    lab = chnode(set%rbnod)
    IF( lab == 0 ) set%rbnod = 0     !master node deleted
    set%nmast = lab
  ELSE
    set%nmast = 0
  END IF
  RETURN

END SUBROUTINE updl10
