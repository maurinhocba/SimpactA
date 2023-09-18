 SUBROUTINE outdy6(b1,b2,nreqs,heade,narch,iwrit,ngrqs,ttime,nstre)
 !********************************************************************
 !
 !**   output routine for shell elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: iwrit,   &!flag for output into ASCII file
                                  nstre,   &!number of stress values
                                  nreqs,   &!number of elements for output
                                  narch     !output file number
 INTEGER (kind=4), POINTER :: ngrqs(:)      !element numbers for output
 REAL (kind=8), INTENT(IN) :: ttime         !total time
 LOGICAL, INTENT(IN) :: b1,  &   !flag for output at selected points
                        b2       !flag for global output
 TYPE(ele06), POINTER  :: heade  !pointer to first element

 TYPE(ele06), POINTER :: e

 INTEGER (kind=4) g,iaux,ielem
 REAL    (kind=8) ehist(2,ngaus)
 LOGICAL :: extra

 ! *** writes stress history at requested gauss points at every toutd time

 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime
   DO g=1,nreqs
     iaux  = ngrqs(g)-1
     ielem = iaux/ngaus+1
     iaux  = MOD(iaux,ngaus)
     e => heade                   !point to first element
     DO                           !search for element with label IELEM
       IF( e%numel == ielem )EXIT !if found, exit search
       e => e%next                !point to next element
     END DO
           !integrated stresses (forces and moments)
     IF( ASSOCIATED (e%ehist) ) THEN
       ehist(1,1) = e%ehist(1,iaux)
       ehist(2,1) = SQRT(e%ehist(2,iaux))
     ELSE
       ehist(1:2,1) = 0d0
     END IF
     WRITE(narch,ERR=9999) e%strsg(1:8,iaux),ehist(1:2,1),e%ambda(1,iaux)
   END DO
   !  CALL flushf(narch)
 END IF

 ! *** writes stresses on output file

 IF(b2) THEN
   IF(iwrit == 1)WRITE(lures,"(5x,'STRESSES FOR SHELL ELEMENTS'/)",ERR=9999)
   extra = nstre == 14
   e => heade                   !point to first element
   DO
     IF( .NOT.ASSOCIATED(e) )EXIT
     ielem = e%numel
     IF( ASSOCIATED (e%ehist) ) THEN
       ehist(1:2,:) = e%ehist(1:2,:)
       DO g=1,ngaus
         ehist(2,g) = SQRT(ehist(2,g))
       END DO
     ELSE
       ehist = 0d0
     END IF

     !         writes complete stress state on tape 16

     DO g=1,ngaus
       IF( extra )THEN
         WRITE(16,ERR=9999) e%strsg(1:8,g),ehist(1:2,g),e%ambda(2,g),e%strsg(9:14,g)
       ELSE
         WRITE(16,ERR=9999) e%strsg(1:8,g),ehist(1:2,g),e%ambda(2,g)
       END IF
     END DO
     IF(iwrit == 1) THEN

       !           PRINT membrane stresses

       WRITE(lures,"(/,6x,'Element No.',i5,/,6x,16('-')/,1x,'g.p.', &
                    &  6x,'N-xx',9x,'N-yy',9x,'N-xy',8x,'Eepd',6x,  &
                    & 'Equi.Stres',3x,'Lambda')",ERR=9999) ielem
       DO g=1,ngaus
         WRITE(lures,"(i3,3x,3e13.5,1x,2e12.5,f8.5)",ERR=9999) g,   &
                       e%strsg(1:3,g),ehist(1:2,g),e%ambda(2,g)
       END DO

       !           PRINT flexional stresses

       WRITE(lures,"(' G.p.',5x,'M-xx',11x,'M-yy',11x,'M-xy',11x, &
                   & 'Q-x',12x,'Q-y')",ERR=9999)
       DO g=1,ngaus
         WRITE(lures,"(i3,1x,3e15.6,2e15.6)",ERR=9999) g,e%strsg(4:8,g)
       END DO
       !           PRINT aditional stresses
       IF(extra)THEN
         WRITE(lures,"(' G.p.',5x,'M-xx',9x,'M-yy',9x,'M-xy',9x,'M-yx',9x, &
                     & 'Q-x',10x,'Q-y')",ERR=9999)
         DO g=1,ngaus
           WRITE(lures,"(i3,1x,6e13.5)",ERR=9999) g,e%strsg(9:14,g)
         END DO
       END IF
     END IF
     e => e%next          !point to next element
   END DO
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outdy6
