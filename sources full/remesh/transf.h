SUBROUTINE transf(numpn, coorn, coora, l_old, v_old, v_new, elemp )

  ! transfers data from the old mesh to the new one (from old to new nodes)

  !USE lispa0
  !USE ctrl_db, ONLY : ndime, ndofn
  IMPLICIT NONE

  ! dummy arguments
  INTEGER (kind=4) :: numpn,     &    ! number of new points
                      l_old(:,:)      ! connectivities of old elements
  INTEGER (kind=4), POINTER, OPTIONAL :: elemp(:) !element where projects
  REAL (kind=8)    :: coorn(:,:),  &  ! coordinates of new points
                      coora(:,:),  &  ! coordinates of old points
                      v_old(:,:)      ! nodal velocities (old mesh)
  REAL (kind=8), POINTER :: &
                      v_new(:,:)      ! nodal velocities (new mesh)

END SUBROUTINE transf
