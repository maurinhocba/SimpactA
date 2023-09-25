 SUBROUTINE mastd3(nsegm,nnseg,lcseg,nhseg,confor,x)
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
 REAL (kind=8), POINTER, OPTIONAL :: x(:,:)
!     local variables
 INTEGER (kind=4) jm,km,i,j,jn,kn,icseg,inseg,ns
 INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/)
 INTEGER (kind=4), ALLOCATABLE :: ni(:,:)

!     Generate standard Data Base for conforming elements

 nhseg = 0   !initializes connected segments
 DO icseg = 1,nsegm                  !for each segment
   outer : DO inseg = 1,nnseg                !for each node in the segment
     IF( nhseg(inseg,icseg) /= 0)CYCLE       !if already considered cycle
     jm = lcseg(nextn(inseg),icseg)          !first node of the side
     km = lcseg(nextn(nextn(inseg)),icseg)   !second node of the side
     ! search the TWIN side
     DO i=icseg+1,nsegm                      !for the rest of the segments
       DO j=1,nnseg                          !for each side
         IF(nhseg(j,i) /= 0)CYCLE            !if already considered cycle
         jn = lcseg(nextn(j),i)              !first node of the side
         IF( jn == km )THEN                  !if second = first
           kn = lcseg(nextn(nextn(j)),i)     !second node of the side
           IF( kn == jm )THEN                !if first = second
             nhseg(inseg,icseg) = i          !this is the side
             nhseg(j,i) = icseg              !
             CYCLE outer                     !Go to next side
           END IF
         END IF
       END DO
     END DO
   END DO outer
 END DO

 IF( .NOT.PRESENT(confor) .OR. .NOT.PRESENT(x) )RETURN

 IF( confor ) THEN
   !....  Check non-paired sides
   !....  for conforming meshes

   ! search for boundary nodes (sides), make a list
   DO icseg=1,nsegm         !for each segment
     loop : DO inseg=1,nnseg          !for each side
       IF(nhseg(inseg,icseg) /= 0 )CYCLE   !if paired side cycle
       i = nextn(inseg)                    !first local node of segment
       j = nextn(i)                        !second local node of segmnet
       jm= lcseg(i,icseg)           !first global node
       km= lcseg(j,icseg)           !second global node
       DO i=icseg+1,nsegm         !for each segment
         DO j=1,nnseg          !for each side
           IF(nhseg(j,i) /= 0 )CYCLE   !if paired side cycle
           jn = nextn(j)                    !first local node of segment
           kn = nextn(jn)                   !second local node of segmnet
           jn= lcseg(jn,i)                  !first global node
           kn= lcseg(kn,i)                  !second global node
           IF( jm == jn .AND. km == kn )THEN
             WRITE(55,"('REVERSED segm:',i6,' side:',i1,' =  segm:',i6,' side:',i1)")icseg,inseg,i,j
             EXIT loop
           END IF
         END DO
       END DO
     END DO loop
   END DO

 ELSE
  !....  computed 'best' neighbour segment for boundary sides
  !....  in Non-conforming meshes
   ns = COUNT(nhseg(1:3,1:nsegm) == 0) !number of non-paired sides
   ALLOCATE ( ni(4,ns) )  !Auxiliar array
   ni = 0                 !initializes array

   ! search for boundary nodes (sides), make a list
   ns = 0                   !initializes position in array
   DO icseg=1,nsegm         !for each segment
     DO inseg=1,nnseg          !for each side
       IF(nhseg(inseg,icseg) /= 0 )CYCLE   !if paired side cycle
       ns = ns+1                           !increase counter
       ni(1,ns) = icseg                    !assign segment number
       ni(2,ns) = inseg                    !assign side number
       i = nextn(inseg)                    !first local node of segment
       j = nextn(i)                        !second local node of segmnet
       ni(3,ns) = lcseg(i,icseg)           !first global node
       ni(4,ns) = lcseg(j,icseg)           !second global node
     END DO
   END DO

   CALL mastda (ns,ni,x,nhseg)
   DEALLOCATE ( ni )

 END IF
 RETURN
 END SUBROUTINE mastd3
