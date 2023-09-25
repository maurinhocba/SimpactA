 SUBROUTINE outdy9 (b1,b2,nelem,ngaus,nreqs,narch,iwrit,ntype, &
                    ngrqs,ttime,heade)
 !********************************************************************
 !
 !**   output routine for beam/shell elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) nelem,ngaus,nreqs,narch,iwrit,ntype
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL    (kind=8) ttime
 LOGICAL   b1,b2
 TYPE(ele09), POINTER :: heade

 INTEGER (kind=4) g,i,ielem
 REAL    (kind=8) epdef(2,ngaus)
 TYPE(ele09), POINTER :: e

 ! ** writes stress history at requested gauss points at every toutd time

 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime
   DO i=1,nreqs
     g = MOD(ngrqs(i)-1,ngaus)+1
     ielem = (ngrqs(i)-1)/ngaus+1
     !search
     e => heade
     DO
       IF( e%numel == ielem )EXIT
       e => e%next
       IF( .NOT.ASSOCIATED (e) )STOP !error
     END DO
     IF(ASSOCIATED( e%ehist )) THEN
       epdef(1,1) = e%ehist(2,g)
       epdef(2,1) = SQRT(e%ehist(3,g))
     ELSE
       epdef(:,1) = 0d0
     END IF
     WRITE(narch,ERR=9999) e%stres(:,g),epdef(:,1),e%sedef(g)
   END DO
   !CALL flushf(narch)
 END IF

 ! *** writes stresses on output file

 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,"(5x,'STRESSES FOR BEAM/SHELL 2-D ELEMENTS'/)",ERR=9999)
   e => heade
   DO ielem=1,nelem
     IF(ASSOCIATED( e%ehist )) THEN
       epdef = e%ehist(1:2,:)
     ELSE
       epdef = 0d0
     END IF

     ! writes complete stress state on tape 16

     DO g=1,ngaus
       epdef(2,g) = SQRT(ABS(epdef(2,g)))
       WRITE(16,ERR=9999) e%stres(:,g),epdef(:,g),e%sedef(g)
     END DO

     IF(iwrit == 1) THEN

       ! PRINT stresses

       IF(ntype == 1) THEN
         WRITE(lures,"(/,6x,'Element No.',i5,/,6x,16('-')/,1x,'G.P.',5x, &
         &  'N',11x,'M',11x,'Q',10x,'EPD',7x,'Eq-Str',6x,'h/h0'          &
         &       )",ERR=9999) ielem
       ELSE
         WRITE(lures,"(/,6x,'Element No.',i5,/,6x,16('-')/,1x,'g.p.',5x,  &
         &  'N1',10x,'N2',10x,'M1',10x,'M2',10x,'Q',11x,'EPD',7x,'Eq-str',&
         &  6x,'h/h0' )",ERR=9999) ielem
       END IF
       DO g=1,ngaus
         WRITE(lures,"(i3,1x,8e12.4)",ERR=9999) g,e%stres(:,g),epdef(1:2,g),e%sedef(g)
       END DO
     END IF
     e => e%next
   END DO
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outdy9
