 SUBROUTINE outd17(b1,b2,iwrit,heade,nreqs,narch,ngrqs,ttime,ngaus)
 !--------------------------------------------------------------------
 !
 !       output routine for 2-D solid element (TLF)
 !
 !--------------------------------------------------------------------
 USE lispa0
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: iwrit,    &!flag for output into ASCII file
                                  nreqs,    &!number of elements for output
                                  ngaus,    &!number of gauss points
                                  narch      !output file number
 INTEGER (kind=4), POINTER :: ngrqs(:)       !element numbers for output
 REAL (kind=8), INTENT(IN) :: ttime
 LOGICAL, INTENT(IN) :: b1,   & !flag for output at selected points
                        b2      !flag for global output
 TYPE (ele17), POINTER  :: heade  !pointer to first element

 INTEGER (kind=4) isec,k,l,m,ielem,osec,mtype
 REAL (kind=8) stres(4),strpl,sxx,syy,sxy,szz,vm,preys
 TYPE (ele17), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: elast

 !     output at selected elements for HISTORY

 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime                  !write process time
   DO l=1,nreqs                        !for each element
     k = MOD(ngrqs(l)-1,ngaus) + 1
     ielem = (ngrqs(l)-1) / ngaus + 1
     e => heade                   !point to first element
     DO                           !search for element with label IELEM
       IF( e%numel == ielem )EXIT !if found, exit search
       e => e%next                !point to next element
     END DO
     stres = e%stint(:,k)     !integrated Kirchhoff stresses
     strpl = 0d0
     ! strpl = e%gausv(5,k)     !plastic strain of first component
     vm    = 0d0              !von Mises value
     preys = 0d0              !actual yield stress
     WRITE(narch,ERR=9999) (stres(m),m=1,4),strpl,vm,preys
   END DO
 END IF
 !
 ! *** writes complete stress state on file 16 for post-process
 !
 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900,ERR=9999)        !if echo to ASCII file
   osec = -1                              !initializes
   e => heade                             !point to first element
   DO                                     !loop over all elements
     IF (.NOT.ASSOCIATED (e) ) EXIT       !if last element processed, exit
     isec = e%matno                      !associated material
     IF( isec /= osec )THEN
       sec => psecs(isec)%p
       mtype = sec%mtbas%mtype             !type of base material
       elast = sec%mtbas%matdef(3) == 1
       osec = isec
     END IF
     !-------------------------------------------------------------
     DO k=1,ngaus
       stres = e%stint(:,k)     !integrated Kirchhoff stresses on material axis
       strpl = 0d0                   !equivalent plastic strain
       vm    = 0d0                   !von Mises equivalent stress
       IF( elast )THEN
         strpl = 0d0            !equivalent plastic strain
       ELSE
         strpl = e%gausv(5,k)   !equivalent plastic strain
       END IF
       vm = (stres(1)+stres(2)+stres(4))/3d0
       sxx = stres(1)-vm         !
       syy = stres(2)-vm         !
       sxy = stres(3)            !Stress XY
       szz = stres(4)-vm         !
       vm = SQRT(1.5d0*(sxx**2+syy**2+szz**2+2d0*sxy**2))  !Von Mises stress
       WRITE(16,ERR=9999) (stres(l),l=1,4),strpl,vm
       IF(iwrit == 1)WRITE(lures,920,ERR=9999) e%numel,k,stres(1:4),strpl,vm
     END DO

     e => e%next                   !point to next element
   END DO
 END IF

 RETURN
 900 FORMAT(5X,'STRESSES FOR 2-D SOLID ELEMENTS in MATERIAL AXES'//   &
 &       ' ELEM G   S-XX',9X,'S-YY',9X,'S-XY',9X,'S-ZZ',9X,'EPST',9X,'VMST')
 920 FORMAT(i5,i2,6e13.4)
 9999 CALL runen2('')
 END SUBROUTINE outd17
