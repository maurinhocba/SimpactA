      SUBROUTINE compare_faces(nnf,nf,nf2,flag)
      !Compare two faces, it checks if containts the same nodes
      IMPLICIT NONE
      INTEGER (kind=4), INTENT(IN) :: nnf,     &  !number of nodes in each face
                                      nf(nnf), &  !nodes on face one
                                      nf2(nnf)    !nodes on face two
      LOGICAL, INTENT(OUT) ::  flag               !.TRUE. if containts the same nodes

      INTEGER (kind=4) :: i,j

      flag = .FALSE.  !initializes
      outer : DO i=1,nnf        !for each node in face one
        DO j=1,nnf                !loop on each node in face two
          IF( nf2(j) == nf(i) ) CYCLE outer  !if node exist loop to next node
        END DO
        RETURN                    !if node not found exit routine
      END DO outer
      flag = .TRUE.             !all nodes checked
      RETURN
      END SUBROUTINE compare_faces
