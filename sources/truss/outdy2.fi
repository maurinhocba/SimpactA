 SUBROUTINE outdy2 (b1,b2,nelem,nreqs,narch,iwrit,heade,     &
                    ngrqs,lambd)
 !********************************************************************
 !
 !**   output routine for 2/3-d truss elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) nelem,nreqs,narch,iwrit
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL (kind=8) lambd
 LOGICAL b1,b2
 TYPE(ele02), POINTER :: heade
 TYPE(section), POINTER :: sec

 INTEGER (kind=4) i,ielem,imat,omat
 REAL (kind=8) stres,epdef,sigma
 LOGICAL :: newm,yield
 TYPE(ele02), POINTER :: e

 ! ** writes stress history at requested gauss points at every toutd time

 IF(b1 .AND. (nreqs > 0)) THEN
   omat = -1
   WRITE(narch,ERR=9999) lambd
   DO i = 1,nreqs
     ielem = ngrqs(i)
     e => heade
     DO
       IF( e%numel == ielem )EXIT
       e => e%next
     END DO
     imat  = e%matno
     newm = imat /= omat
     IF ( newm )THEN
       sec => psecs(imat)%p
       omat = imat
       yield = sec%mtbas%matdef(3) > 1
     END IF
     sigma = e%gausv(3)
     stres = e%gausv(4)
     IF( yield ) THEN
       epdef = e%gausv(5)
     ELSE
       epdef = 0d0
     END IF
     WRITE(narch,ERR=9999) stres,sigma,epdef
   END DO
    !CALL flushf(narch)
 END IF

 ! *** writes stresses on output file

 IF(b2) THEN
   IF(iwrit == 1)WRITE(lures,"(5x,'stresses for truss elements'/)",ERR=9999)
   omat = -1
   e => heade
   DO ielem = 1,nelem
     imat  = e%matno
     newm = imat /= omat
     IF ( newm )THEN
       sec => psecs(imat)%p
       omat = imat
       yield = sec%mtbas%matdef(3) > 1
     END IF
     sigma = e%gausv(3)
     stres = e%gausv(4)
     IF( yield ) THEN
       epdef = e%gausv(5)
     ELSE
       epdef = 0d0
     END IF
     !         writes complete stress state on tape 16

     WRITE(16,ERR=9999) stres,sigma,epdef
     IF(iwrit == 1) THEN

       !         PRINT stresses

       WRITE(lures,"(5x,'element No.',i5,'    N =',e13.5,'   S =', &
&                 e13.5,'    EPD=',e12.4)",ERR=9999) ielem,stres,sigma,epdef
     END IF
     e => e%next
   END DO
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outdy2
