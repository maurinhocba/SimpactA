 SUBROUTINE outdy1 (b1,b2,nnode,nelem,nreqs,narch,iwrit,heade,ndime,     &
                    ngrqs,lambd)
 !********************************************************************
 !
 !**   output routine for 2/3-d SPOT elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) nnode,nelem,nreqs,narch,iwrit,ndime
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL (kind=8) lambd
 LOGICAL b1,b2
 TYPE(ele01), POINTER :: heade

 INTEGER (kind=4) i,ielem
 REAL (kind=8) momen,force,epdef
 LOGICAL :: same
 TYPE(ele01), POINTER :: e

 ! ** writes momens history at requested gauss points at every toutd time

 epdef = 0d0
 IF(b1 .AND. (nreqs > 0)) THEN
   !omat = -1
   WRITE(narch,ERR=9999) lambd
   DO i = 1,nreqs
     ielem = ngrqs(i)
     e => heade
     DO
       IF( e%numel == ielem )EXIT
       e => e%next
     END DO
     !imat  = e%matno
     !newm = imat /= omat
     !IF ( newm )THEN
     !  CALL sect_search (imat,found,sec)
     !  omat = imat
     !  yield = sec%mtbas%matdef(3) > 1
     !END IF
     IF( nnode == 2) THEN
       same = ( e%lnods(2) == 0 )
       force = e%gausv(1)            !force
       IF( e%rot )THEN
         IF( same )THEN
           momen = e%gausv(2+ndime)
         ELSE
           momen = e%gausv(2)
         END IF
       ELSE
         momen = 0d0
       END IF
       !IF( yield ) THEN
       !  epdef = e%gausv(?)
       !ELSE
       !  epdef = 0d0
       !END IF
       WRITE(narch,ERR=9999) force,momen,epdef
     ELSE
       WRITE(narch,ERR=9999) e%gausv(7:12)
     END IF
   END DO
    !CALL flushf(narch)
 END IF

 ! *** writes forces on output file

 IF(b2) THEN
   IF(iwrit == 1)WRITE(lures,"(5x,'forces for SPOT elements'/)",ERR=9999)
   !omat = -1
   e => heade
   DO ielem = 1,nelem
     !IF ( newm )THEN
     !  CALL sect_search (imat,found,sec)
     !  omat = imat
     !  yield = sec%mtbas%matdef(3) > 1
     !END IF
     IF( nnode == 2)THEN
       same = ( e%lnods(2) == 0 )
       force = e%gausv(1)            !force
       IF( e%rot )THEN
         IF( same )THEN
           momen = e%gausv(2+ndime)
         ELSE
           momen = e%gausv(2)
         END IF
       ELSE
         momen = 0d0
       END IF
       !IF( yield ) THEN
       !  epdef = e%gausv(?)
       !ELSE
       !  epdef = 0d0
       !END IF
       !         writes complete momens state on tape 16

       WRITE(16,ERR=9999) force,momen,epdef
       IF(iwrit == 1) THEN

         !         PRINT moments

         WRITE(lures,"(5x,'element No.',i5,'    N =',e13.5,'   M =', &
&                   e13.5,'    EPD=',e12.4)",ERR=9999) ielem,force,momen,epdef
       END IF
     ELSE
       momen = SQRT(e%gausv(18)**2+e%gausv(19)**2)
       WRITE(16,ERR=9999) e%gausv(17),momen,0d0
       IF(iwrit == 1) WRITE(lures,"(5x,'element No.',i5,        &
                   & '  N =',e13.5,' T1 =',e13.5,' T2 =',e13.5, &
                   & ' R1 =',e13.5,' R2 =',e13.5,' R3 =',e13.5)",ERR=9999) ielem,e%gausv(17:22)

     END IF
     e => e%next
   END DO
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outdy1
