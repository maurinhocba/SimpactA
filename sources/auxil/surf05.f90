      SUBROUTINE surf05 ( lnods, nelem, nnode)
!******************************************************************
!
!*** get the boundary definition from the element set
!    segments are oriented (outward normal)
!
!******************************************************************
      USE surf_db
      IMPLICIT NONE
      !dummy arguments
      INTEGER (kind=4), INTENT(IN) :: nelem,  &  !number of elements
      &                               nnode      !number of nodes per element
      INTEGER (kind=4), INTENT(IN) :: lnods(:,:) !(nnode,nelem) connectivities
      !TYPE (cont_srf), POINTER :: surfa   !INTENT(OUT) surface data

      ! local variables
      INTEGER (kind=4) :: in,ie,je,jn,nf(4),nnf,nf2(4),nnf2,ne,i,nface,nns
      INTEGER, ALLOCATABLE :: nn(:,:), cn(:)
      INTEGER (kind=4), ALLOCATABLE :: lface(:,:)
      TYPE (srf_seg), POINTER :: seg
      LOGICAL :: flag

      IF(nnode == 4 ) THEN
        nface = 4
        nns = 3
        ALLOCATE ( nn(nns,nface), cn(nns) )
        nn = RESHAPE( (/ 1,3,2, 1,2,4, 2,3,4, 1,3,4 /), (/ 3,4 /) )
        cn = (/ 2,3,1 /)
      ELSE
        nface = 6
        nns = 4
        ALLOCATE ( nn(nns,nface), cn(nns) )
        nn = RESHAPE( (/ 1,4,3,2, 1,2,6,5, 2,3,7,6,  &
                         3,4,8,7, 4,1,5,8, 5,6,7,8 /), (/ 4,6 /) )
        cn = (/ 2,3,4,1 /)
      END IF

      ALLOCATE (lface(nface,nelem))  !this arrays stores -1 if inner face
                                     !  > 1 if part of the boundary
      lface = 0                      ! initializes to unknown
      !GENERATE SURFACE data_base Initialize empty list
      IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
      NULLIFY(surfa%head,surfa%tail)

      !   SEARCH FOR THE 'TWIN' FACES
      ne = 0                     !initializes number of boundary faces
      DO ie=1,nelem              !for each element
        jsearch : DO in=1,nface        !for each face in the element
          IF( lface(in,ie) /=0 )CYCLE  !if processed face then cycle
          nf=lnods(nn(:,in),ie)        !nodes of the face
          nnf = 1
          DO i=1,nns                   !for each node of the face
            IF( nf(nnf) /= nf(cn(i)) )THEN  !compare with next node
              nnf = nnf+1                 !increase number of nodes of the face
            ELSE                          !if nodes are equal
              nf(cn(nnf)) = nf(cn(cn(i)))  !modify next node
            END IF
          END DO
          nnf = nnf-1
          IF( nnf <= 2 ) THEN      !if face has two nodes or less ==> not a face
            lface(in,ie) = -1      !set as if where an inner face
            CYCLE                  !next face
          END IF
          !   search for the 'twin' face defined by nf(1:nnf)
          DO je=ie+1,nelem  !for each remaining element
            DO jn=1,nface                 !for each face
              IF( lface(jn,je) /= 0 )CYCLE  !if processed face then cycle
              nf2=lnods(nn(:,jn),je)        !nodes of the face
              nnf2 = 1
              DO i=1,nns                    !for each node of the face
                IF( nf2(nnf2) /= nf2(cn(i)) )THEN  !compare with next node
                  nnf2 = nnf2+1            !increase number of nodes of the face
                ELSE                               !if nodes are equal
                  nf2(cn(nnf2)) = nf2(cn(cn(i)))  !modify next node
                END IF
              END DO
              nnf2 = nnf2-1
              IF( nnf2 <= 2 ) THEN      !if face has two nodes or less ==> not a face
                lface(jn,je) = -1      !set as if where an inner face
                CYCLE                  !next face
              END IF
              IF( nnf /= nnf2 )CYCLE   !faces has different number of nodes
              CALL compare_faces(nnf,nf,nf2,flag)
              IF( flag )THEN
                lface(in,ie) = -1         !set to Inner face
                lface(jn,je) = -1         !set to Inner face
                CYCLE jsearch             !TWIN face found, exit search
              END IF
            END DO
          END DO
          ne = ne + 1         !increase number of boundary faces
          lface(in,ie) = ne   !assign number to boundary face
          ! store segment in list
          CALL new_seg(seg)                   !allocate surface segment
          seg%nodes(1:nnf) = nf(1:nnf)        !
          !WRITE(lures, '(10x,4i8)',ERR=9999) seg%nodes(1:nnf)
          CALL add_seg(seg, surfa%head, surfa%tail)  !add segment to the d_b
        END DO jsearch
      END DO

      surfa%nelem = ne                        !keep number of segments
      DEALLOCATE (lface)                      !release memory

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE surf05
