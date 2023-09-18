SUBROUTINE del_e_mrea(mrea,dim,elo,elf)
!-----------------------------------------------------------------------
! Resize a vector of integers maintaining old data
! Used in QLOOP for solid remeshing
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: dim,         & !Dimension to modify
                               elo            !First element to delete in array
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: elf   !Last element to delete in array
  REAL(kind=8),POINTER:: mrea(:,:)  !Array to modify
  !--- Local variables
  INTEGER(kind=4):: ndim, ndim1, ndim2, nno, nnf, ndl, ndl1, ndl2, iaux
  REAL(kind=8),ALLOCATABLE:: work(:,:)

  IF (.NOT.ASSOCIATED(mrea)) RETURN  !Array must be allocated

  IF ((dim /= 1).AND.(dim /= 2)) CALL runen3('Internal error in subroutine "del_e_mrea".')  !Internal error: wrong dimension

  ndim1 = SIZE(mrea,DIM=1)    !Determine number of files in array
  ndim2 = SIZE(mrea,DIM=2)    !Determine number of columns in array
  IF (dim == 1) THEN
    ndim = ndim1
  ELSE IF (dim == 2) THEN
    ndim = ndim2
  END IF

  nno = MAX0(1,elo)        !Check first element to delete or select default value
  IF (PRESENT(elf)) THEN   !Check last element to delete or select default value
    nnf = MIN0(ndim,elf)
  ELSE
    nnf = nno   !Only delete one element defined by the first element
  END IF

  IF (nno > nnf) THEN
    iaux = nno  ;  nno = nnf  ;  nnf = iaux   !Sort limit values
  END IF

  ndl = nnf - nno + 1   !Number of elements to delete
  IF (dim == 1) THEN
    ndl1 = ndl
    ndl2 = 0
  ELSE IF (dim == 2) THEN
    ndl1 = 0
    ndl2 = ndl
  END IF

  ALLOCATE(work(ndim1,ndim2))   !Allocate auxiliar array

  work(1:ndim1,1:ndim2) = mrea(1:ndim1,1:ndim2)   !Store old data in the auxiliar array

  DEALLOCATE(mrea)                        !ARRAY resize
  IF (ndl == ndim) RETURN                 !   <--- All elements has been deleted (no more calc. needed)
  ALLOCATE(mrea(ndim1-ndl1,ndim2-ndl2))   !   <--- Some elements has been deleted

  IF (dim == 1) THEN
    IF (nno > 1) mrea(1:nno-1,1:ndim2)=work(1:nno-1,1:ndim2)                 !Restore old data in
    IF (nnf < ndim1) mrea(nno:ndim1-ndl1,1:ndim2)=work(nnf+1:ndim1,1:ndim2)  !the array
  ELSE IF (dim == 2) THEN
    IF (nno > 1) mrea(1:ndim1,1:nno-1)=work(1:ndim1,1:nno-1)                 !Restore old data in
    IF (nnf < ndim2) mrea(1:ndim1,nno:ndim2-ndl2)=work(1:ndim1,nnf+1:ndim2)  !the array
  END IF

  DEALLOCATE(work)       !Deallocate auxiliar vector

RETURN
END SUBROUTINE del_e_mrea
