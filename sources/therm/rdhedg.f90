 SUBROUTINE rdhedg(nedge,iwrit,heads,tails)

 !     Reads edge distributed flows (solid 2D elements)

 USE c_input
 USE heat_db
 IMPLICIT NONE
 INTEGER (kind=4) :: nedge,iwrit
 TYPE (edg_nod),POINTER :: heads, tails
 !local variables
 INTEGER (kind=4), PARAMETER :: nnmax = 3 !maximum nodes in segment
 INTEGER (kind=4) :: i,kount,nn,noprs(nnmax),nnode
 REAL (kind=8) :: p(nnmax)

 TYPE (edg_nod),POINTER :: edg

 LOGICAL :: flag = .FALSE.

 !*** distributed edge heat flows
 CALL listen ('RDHEDG')

 nedge = 0
 IF ( .NOT.exists('EDGE  ')) THEN
   backs = .TRUE.

 ELSE
   IF(iwrit == 1) WRITE(lures,"(/6x,'List of Edges and Applied Heat Flows'/)",ERR=9999)
   !*** loop over each face with a heat flow
   CALL ini_edg (heads, tails)

   !*** evaluate nodes in segment and default system for flows definition
   IF( exists('NNODE ',i))THEN
     nnode = INT(param(i)) !nodes in segment
   ELSE
     nnode = 2  !default: 2-node defined segment
   END IF

   !*** loop over each loaded edge
   p = 0d0    !initializes
   noprs = 0
   DO
     CALL listen ('RDHEDG')
     IF (exists('ENDEDG')) EXIT

     !*** READ DATA locating the heat edge and applied flow
     nedge = nedge + 1

     IF (exists('N1    ',i)) THEN
       noprs(1) = INT(param(i))    !read first node
       IF (exists('N2    ',i)) THEN
         noprs(2) = INT(param(i))   !read second node
         nn = 2
       ELSE
         CALL runend('RDHEDG: Specifiy Second node N2')
       END IF
       IF (exists('N3    ',i)) THEN                 ! two node segment are allowed
         noprs(3) = INT(param(i))  !read third node
         nn = 3
       END IF
     ELSE
       nn = nnode  !use default node definition
       noprs(1:nn) = INT(param(i+1:i+nn)) !read all nodes in segment
       i = i+nn
     END IF

     IF( nwopa > i )THEN  !more arguments in the line
       kount = nwopa - i
       CALL vecasi(kount,param(i+1),p)
     ELSE                 !read flows in the next line
       CALL rdfrre('RDHEDG',p,kount,nn,flag)
     END IF

     IF (iwrit == 1) THEN
       WRITE(lures,"(' EDGE',i3,' NN=',i2,' Nodes',6i6)",ERR=9999)  &
                    nedge,nn,noprs(1:nn)
       WRITE(lures,"(' Nodal p ',3g12.3)",ERR=9999) p(1:nn)
     END IF

     ! add new segment to list
     ALLOCATE (edg)
     edg%nn2   = nn
     edg%flow2 = p
     edg%lnod2 = noprs
     CALL add_edg (edg, heads, tails)

    END DO

  ENDIF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rdhedg
