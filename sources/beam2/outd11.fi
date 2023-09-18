 SUBROUTINE outd11 (b1,b2,nelem,nstre,nreqs,narch,ntype,ngaus,iwrit, &
                    ngrqs,ttime,heade,stint,sh)
 !********************************************************************
 !
 !**   output routine for beam/shell elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) nelem,nstre,nreqs,narch,ntype,ngaus,iwrit
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL    (kind=8) ttime,sh(2,2)
 REAL    (kind=8), POINTER :: stint(:,:)
 LOGICAL   b1,b2
 TYPE(ele11), POINTER :: heade

 INTEGER (kind=4) i,g,ielem,nlayr,ns
 REAL    (kind=8) epdef(2),stres(nstre+1),pp(2,2)
 TYPE(ele11), POINTER :: e

 ! ** writes stress history at requested gauss points at every toutd time

 ns = nstre/2                  !number of components
 i = 1                               !initializes position
 DO g=1,ngaus
   pp(1,g) = i                       !N_1 position for this gauss point
   pp(2,g) = pp(1,g)+nstre           !M_1 position for this gauss point
   i = pp(2,g) + 1                   !add shear position
 END DO

 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime
   DO i=1,nreqs
     g = MOD(ngrqs(i)-1,ngaus) + 1  !g= MOD(ngrqs(i),10)
     ielem = (ngrqs(i)-1)/ngaus + 1  !ielem = ngrqs(i)/10
     !search
     e => heade
     DO
       IF( e%numel == ielem )EXIT
       e => e%next
       IF( .NOT.ASSOCIATED (e) )STOP !error
     END DO
     IF(ASSOCIATED( e%gausv )) THEN
       nlayr = SIZE(e%gausv,DIM=2)
       epdef(1) = e%gausv(ns+1,1,g)
       epdef(2) = e%gausv(ns+1,nlayr,g)
     ELSE
       epdef = 0d0
     END IF
     WRITE(narch,ERR=9999) stint(pp(1,g):pp(2,g),ielem),epdef,e%lambd(2)
   END DO
   !CALL flushf(narch)
 END IF

 ! *** writes stresses on output file

 IF(b2) THEN
   IF(iwrit == 1) THEN
     WRITE(lures,"(5x,'STRESSES FOR BEAM/SHELL 2-D ELEMENTS'/)",ERR=9999)

     IF(ntype == 1) THEN
       WRITE(lures,"(/,' Element No.',6x,'N',11x,'M',10x,'Q',9x,'EPS-1',7x,     &
                   &   'EPS-N',6x,'h/h0')",ERR=9999)
     ELSE IF(ntype == 2 .OR. ntype == 3 )THEN
       WRITE(lures,"(/,' Element No.',6x,'N1',10x,'N2',10x,'M1',10x,    &
                   &   'M2',10x,'Q',9x,'EPS-1',7x,'EPS-2',6x,'h/h0' )",ERR=9999)
     END IF
   END IF

   e => heade
   DO ielem=1,nelem

     DO g=1,ngaus
       IF(ASSOCIATED( e%gausv )) THEN
         nlayr = SIZE(e%gausv,DIM=2)
         epdef(1) = e%gausv(ns+1,1,g)
         epdef(2) = e%gausv(ns+1,nlayr,g)
       ELSE
         epdef = 0d0
       END IF

       ! writes complete stress state on tape 16
       IF( ngaus == 2 )THEN
         stres = sh(1,g)*stint(pp(1,1):pp(2,1),ielem)+sh(2,g)*stint(pp(1,2):pp(2,2),ielem)
         WRITE(16,ERR=9999) stres,epdef,e%lambd(2)
       ELSE
         WRITE(16,ERR=9999) stint(pp(1,g):pp(2,g),ielem),epdef,e%lambd(2)
       END IF
       IF(iwrit == 1) &
         WRITE(lures,"(i10,2x,8e12.4)",ERR=9999) ielem,stint(pp(1,g):pp(2,g),ielem),epdef,e%lambd(2)
     END DO
     e => e%next
   END DO
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outd11
