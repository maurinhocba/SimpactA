 SUBROUTINE outd20(b1,b2,iwrit,heade,nreqs,narch,ngrqs,ttime)
 !--------------------------------------------------------------------
 !
 !       output routine for 2-D solid triangular element (TLF)
 !
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: iwrit,   &!flag for output into ASCII file
                                  nreqs,   &!number of elements for output
                                  narch     !output file number
 INTEGER (kind=4), POINTER ::     ngrqs(:)  !element numbers for output
 REAL (kind=8), INTENT(IN) :: ttime         !total time
 LOGICAL, INTENT(IN) :: b1,  &   !flag for output at selected points
                        b2       !flag for global output
 TYPE (ele20), POINTER  :: heade  !pointer to first element

 INTEGER (kind=4) l,ielem,m,pos
 REAL (kind=8) stres(4),sxx,syy,sxy,szz,vm,preys,strpl
 TYPE (ele20), POINTER :: e

 !     output at selected elements for HISTORY

 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime         !write process time
   DO l=1,nreqs                        !for each point
     ielem  = ngrqs(l)         !element number
     e => heade                        !point to first element
     DO                           !search for element with label IELEM
       IF( e%numel == ielem )EXIT !if found, exit search
       e => e%next                !point to next element
     END DO
     pos = 0
     IF(ASSOCIATED(e%gausv)) pos = SIZE( e%gausv )
     stres = e%stint     !integrated Kirchhoff stresses on material axis

     IF( pos == 0 )THEN
       strpl = 0d0            !equivalent plastic strain
     ELSE
       strpl = e%gausv(pos)        !equivalent plastic strain
     END IF
     vm = (stres(1)+stres(2)+stres(4))/3d0
     sxx = stres(1)-vm         !
     syy = stres(2)-vm         !
     sxy = stres(3)            !Stress XY
     szz = stres(4)-vm         !
     vm = SQRT(1.5d0*(sxx**2+syy**2+szz**2+2d0*sxy**2))  !Von Mises stress
     WRITE(narch,ERR=9999) (stres(m),m=1,4),strpl,vm,preys
   END DO
 END IF
 !
 ! *** writes complete stress state on file 16 for post-process
 !
 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900,ERR=9999)        !if echo to ASCII file
   e => heade                             !point to first element
   DO                                     !loop over all elements
     IF (.NOT.ASSOCIATED (e) ) EXIT       !if last element processed, exit
     pos = 0
     IF(ASSOCIATED(e%gausv)) pos = SIZE( e%gausv )

     !-------------------------------------------------------------
     !
     !      WRITE ASCII    Nxx, Nyy, Nxy, Mxx, Myy, Mxy, t/t0
     !      WRITE for GiD  Nxx, Nyy, Nxy, Mxx, Myy, Mxy,
     !                     t/t0, ep1, epN, VMises1, VMisesN
     !
     !-------------------------------------------------------------

     stres = e%stint     !integrated Kirchhoff stresses on material axis

     IF( pos == 0 )THEN
       strpl = 0d0            !equivalent plastic strain
     ELSE
       strpl = e%gausv(pos)        !equivalent plastic strain
     END IF
     vm = (stres(1)+stres(2)+stres(4))/3d0
     sxx = stres(1)-vm         !
     syy = stres(2)-vm         !
     sxy = stres(3)            !Stress XY
     szz = stres(4)-vm         !
     vm = SQRT(1.5d0*(sxx**2+syy**2+szz**2+2d0*sxy**2))  !Von Mises stress
     stres(1:4) = stres(1:4)

     WRITE(16,ERR=9999) (stres(l),l=1,4),strpl,vm
     IF(iwrit == 1)WRITE(lures,"(i7,6e12.4)",ERR=9999) e%numel,stres,strpl,vm

     e => e%next  !point to next element
   END DO
 END IF

 RETURN
 900 FORMAT(5X,'STRESSES FOR 2-D TRIANGULAR  ELEMENTS'//  &
&       '   ELEM   S-XX',9X,'S-YY',9X,'S-XY',   &
&              9X,'S-ZZ',9X,'EPST',9X,'VMS')
 9999 CALL runen2('')
 END SUBROUTINE outd20
