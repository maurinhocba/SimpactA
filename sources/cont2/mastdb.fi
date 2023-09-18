      SUBROUTINE mastdb(nsegm,nnseg,lcseg,nhseg,confor,x)
!
!   Generates data base for surfaces that will act as master in some pairs.
!   This data includes the neighbours segments
!
!.... input
!....   nsegm = number of segments of the master surface considered
!....   nnseg = number of nodes per segment of the surface considered
!....   lcseg(inseg,icseg) = global node number for the local element node
!....                        [inseg] of the segment [icseg]
!....   confor : .TRUE.   conforming surface is considered
!....            .FALSE.  non-conforming surface is considered
!....   x     = nodal coordinates
!.... output
!....   nhseg(nnseg,nsegm) = neighbours segments to each side

      IMPLICIT NONE
!     Dummy arguments
      INTEGER (kind=4), INTENT(IN) :: nsegm,nnseg,lcseg(:,:)
      INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
      LOGICAL, INTENT(IN), OPTIONAL :: confor
      REAL (kind=8), INTENT(IN), OPTIONAL :: x(:,:)
!     local variables
      INTEGER (kind=4) jm,i,j,jn,icseg,inseg,ns
      INTEGER (kind=4), PARAMETER :: nextn(2) = (/2,1/)
      INTEGER (kind=4), ALLOCATABLE :: ni(:,:)

!     Generate standard Data Base for conforming elements

      nhseg = 0   !initializes connected segments
      DO icseg = 1,nsegm                  !for each segment
        outer : DO inseg = 1,nnseg                !for each node in the segment
          IF( nhseg(inseg,icseg) /= 0)CYCLE       !if already considered cycle
          jm = lcseg(inseg,icseg)                 !node
          ! search the TWIN node
          DO i=icseg+1,nsegm                      !for the rest of the segments
            DO j=1,nnseg                          !for each node
              IF(nhseg(j,i) /= 0)CYCLE            !if already considered cycle
              jn = lcseg(j,i)                     !node
              IF( jn == jm )THEN                  !if same node
                nhseg(inseg,icseg) = i            !this is the node
                nhseg(j,i) = icseg                !
                CYCLE outer                       !Go to next node
              END IF
            END DO
          END DO
        END DO outer
      END DO

      IF( .NOT.PRESENT(confor) .OR. .NOT.PRESENT(x) )RETURN
      IF( confor )RETURN

!....  computed 'best' neighbour segment for boundary sides
!....  in Non-conforming meshes

      ns = COUNT(nhseg(1:2,1:nsegm) == 0) !number of non-paired nodes
      IF( ns == 0 )RETURN
      ALLOCATE ( ni(4,ns) )  !Auxiliar array
      ni = 0                 !initializes array

      ! search for boundary nodes (sides), make a list
      ns = 0                   !initializes position in array
      DO icseg=1,nsegm         !for each segment
        DO inseg=1,nnseg          !for each node
          IF(nhseg(inseg,icseg) /= 0 )CYCLE   !if paired side cycle
          ns = ns+1                           !increase counter
          ni(1,ns) = icseg                    !assign segment number
          ni(2,ns) = inseg                    !assign side number
          ni(3,ns) = lcseg(inseg,icseg)       !global node
          ni(4,ns) = lcseg(nextn(inseg),icseg)       !opposite global node
        END DO
      END DO

      !     compare boundaries and assign alternate nodes

      CALL mastd2 (ns,ni,x,nhseg)
      DEALLOCATE ( ni )

      RETURN
      END SUBROUTINE mastdb
