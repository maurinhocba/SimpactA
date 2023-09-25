 SUBROUTINE rdhsur(nsurf,iwrit,heads,tails)

 !     Reads surface distributed flows

 USE c_input
 USE heat_db
 IMPLICIT NONE
 INTEGER (kind=4) :: nsurf,iwrit
 TYPE (srf_nod),POINTER :: heads, tails

 INTEGER (kind=4), PARAMETER :: nnmax = 6
 INTEGER (kind=4) :: i,n,nn,nvert,ivert,kount,noprs(nnmax),sp
 REAL (kind=8) :: p(nnmax)
 CHARACTER (len=6), PARAMETER :: pos(2) = (/'BOTTOM','TOP   '/)

 TYPE (srf_nod),POINTER :: srf

 LOGICAL :: flag = .FALSE.

 !*** distributed surface heat flows

 CALL listen ('RDHSUR')

 nsurf = 0

 IF( .NOT.exists('SURF  ')) THEN
   backs = .TRUE. !
 ELSE

   IF(iwrit == 1) WRITE(lures,"(/6x,'List of Faces read in this strategy'/)",ERR=9999)

   !*** loop over each face with a heat flow

   CALL ini_srf (heads, tails)
   p = 0d0    !initializes
   noprs = 0
   DO
     CALL listen('RDHSUR')
     IF (exists('ENDSUR')) EXIT

     nsurf = nsurf + 1
     IF(exists('NN    ',i))THEN
       nn = INT(param(i))
       IF( nn /= 3 )CALL runend('NN must be three now')
       DO n=1,nn
         IF( words(i+n) == '      ')THEN
           noprs(n) = INT(param(i+n))
         ELSE
           CALL runend('RDHSUR: ERROR reading face nodes   ')
         END IF
       END DO
       nvert = nn
       IF( nn > 4) nvert = nn/2
     ELSE
       CALL runend('DSURF0: Specifiy Number of nodes NN')
     END IF

     sp = 1        !top surface is the default
     IF(exists('TOP   ')) sp = 1
     IF(exists('BOTTOM')) sp = 0

     CALL rdfrre('RDHSUR',p,kount,nvert,flag)
     IF( kount == 1)THEN
         p(2:nn) = p(1)
     ELSE IF ( nn > nvert )THEN
       DO ivert = nvert+1,nn-1
         p(ivert) = (p(ivert-nvert)+p(ivert-nvert+1))/2d0
       END DO
       p(nn) = (p(1)+p(nvert))/2d0
     END IF

     !         prints actual DATA
     IF(iwrit == 1) THEN
       WRITE(lures,"(a6,' FACE',i3,' NN=',i2,' Nodes',6i6)",ERR=9999)  &
                    pos(sp+1),nsurf,nn,noprs(1:NN)
       WRITE(lures,"(' Nodal p ',6g12.3)",ERR=9999) p(1:nn)
     END IF

     ALLOCATE (srf)

     srf%pos   = sp == 1
     srf%nn3   = nn
     srf%flow  = p
     srf%lnod3 = noprs

     CALL add_srf (srf, heads, tails)

   END DO

 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rdhsur
