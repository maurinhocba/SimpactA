SUBROUTINE transf1(ne_new,nnod,ng,coorn,coora,l_old,l_new,sg_new,strnd)

  ! transfers data from the old mesh to the new one (from old nodes to new Gauss points)

  IMPLICIT NONE

  ! dummy arguments
  INTEGER (kind=4), INTENT(IN) :: ne_new, &       ! number of elements in the new mesh
                                  nnod,   &       ! number of nodes per elemen (new)
                                  ng              ! number of Gauss points
  INTEGER (kind=4), POINTER    :: l_old(:,:),  &  ! connectivities of old elements
                                  l_new(:,:)      ! connectivities of new elements
  REAL (kind=8), INTENT(IN)    :: coorn(:,:),  &  ! coordinates of new points
                                  coora(:,:)      ! coordinates of old points

  REAL (kind=8), POINTER :: &
                      sg_new(:,:),  &  ! internal variables at gauss points (new mesh)
                      strnd(:,:)       ! internal variables at nodes (old mesh)

END SUBROUTINE transf1
