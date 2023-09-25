 SUBROUTINE outdy8 (b1,b2,nelem,nreqs,narch,iwrit,ngaus, &
                    ngrqs,ttime,heade)
 !********************************************************************
 !
 !**   output routine for beam elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem,nreqs,narch,iwrit,ngaus
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL    (kind=8) :: ttime
 LOGICAL  :: b1,b2
 TYPE (ele08), POINTER :: heade

 INTEGER (kind=4), PARAMETER :: nstre = 6

 INTEGER (kind=4) g,i,ielem
 TYPE (ele08), POINTER :: e

 ! ** writes stress history at requested gauss points at every toutd time

 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime
   DO g=1,nreqs
     i = MOD(ngrqs(g)-1,ngaus)+1
     ielem = (ngrqs(g)-i)/ngaus+1
     e => heade
     DO
       IF( e%numel == ielem) EXIT
       e => e%next
     END DO
     WRITE(narch,ERR=9999) e%stres(:,i),e%epdef(i),e%sedef(i)
   END DO
   !CALL flushf(narch)
 END IF

 ! *** writes stresses on output file

 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900,ERR=9999)
   e => heade
   DO ielem=1,nelem

     ! writes complete stress state on tape 16

     DO g=1,ngaus
       WRITE(16,ERR=9999) (e%stres(i,g),i=1,nstre),e%epdef(g),e%sedef(g)
     END DO
     IF(iwrit == 1) THEN

       ! PRINT stresses

       WRITE(lures,910,ERR=9999) e%numel
       DO g = 1,ngaus
         WRITE(lures,920,ERR=9999) g,(e%stres(i,g),i=1,nstre),e%epdef(g),e%sedef(g)
       END DO
     END IF
     e => e%next
   END DO
 END IF

 RETURN
 900 FORMAT(5x,'stresses for beam elements'/)
 910 FORMAT(/,6x,'element no.',i5,/,6x,16('-')/,                   &
   1x,'g.p.',5x,'n',10x,'q2',10x,'q3',10x,'mt',10x,'m2',10x,'m3',  &
   10x,'epd',7x,'a/a0' )
 920 FORMAT(i3,1x,7e12.4,f8.5)
 9999 CALL runen2('')
 END SUBROUTINE outdy8
