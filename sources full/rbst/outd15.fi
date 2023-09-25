 SUBROUTINE outd15(b1,b2,logst,iwrit,heade,nreqs,narch,ngrqs,ttime,stint,rhead)
 !--------------------------------------------------------------------
 !
 !       output routine for shell element BBST (TLF)
 !
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), POINTER :: ngrqs(:)      !element numbers for output
 INTEGER (kind=4), INTENT (IN) :: iwrit,   &!flag for output into ASCII file
                                  nreqs,   &!number of elements for output
                                  narch     !output file number
 REAL (kind=8), INTENT(IN) :: ttime         !total time
 REAL (kind=8), POINTER :: stint(:,:) !forces, moments and shear
 LOGICAL, INTENT(IN) :: logst, & ! large or small strain formulation
                        b1,  &   !flag for output at selected points
                        b2       !flag for global output
 TYPE (ele15), POINTER  :: heade  !pointer to first element
 TYPE (reinf), POINTER  :: rhead  !pointer to first element

 INTEGER (kind=4) nlayr,isec,l,m,ielem,ierr,osec,mtype,nvarn,nvare
 REAL (kind=8) stres(8),strpl(2),sigma(6),g(3),lb(2),r1,r2, &
               vm(2),thick,thic0
 TYPE (ele15), POINTER :: e
 TYPE (reinf), POINTER :: rf
 TYPE (section), POINTER  :: sec  !pointer to first element
 REAL (kind=8), ALLOCATABLE :: auxil(:)
 LOGICAL :: found,plast

 !INTERFACE
 !   INCLUDE 'ud_sh2_out.h'
 !END INTERFACE


 !     output at selected elements for HISTORY
 IF(b1 .AND. (nreqs > 0)) THEN
   WRITE(narch,ERR=9999) ttime         !write process time
   DO l=1,nreqs                        !for each element
     ielem  = ngrqs(l)                 !element number
     e => heade                        !point to first element
     DO                           !search for element with label IELEM
       IF( e%numel == ielem )EXIT !if found, exit search
       e => e%next                !point to next element
     END DO
     sigma = stint(1:6,ielem)         !integrated stresses (forces and moments)
     stres(7:8) = stint(7:8,ielem)    !integrated stresses shears
     g     = e%stra1(1:3)         !1st Fundamental form at mid-surface
    ! ******* consistent modifications follow
     IF( logst )THEN
       CALL lgst14(g,r1,r2,lb,'OUTD13',ierr)   !compute eigenvalues and eigenvectors
      ! Computes 2ndP-K stress on the natural Frame
       stres(1) = sigma(1)*r1*r1+sigma(2)*r2*r2+2d0*sigma(3)*r1*r2
       stres(2) = sigma(1)*r2*r2+sigma(2)*r1*r1-2d0*sigma(3)*r1*r2
       stres(3) =(sigma(2)-sigma(1))*r1*r2+sigma(3)*(r1*r1-r2*r2)
       stres(4) = sigma(4)*r1*r1+sigma(5)*r2*r2+2d0*sigma(6)*r1*r2
       stres(5) = sigma(4)*r2*r2+sigma(5)*r1*r1-2d0*sigma(6)*r1*r2
       stres(6) =(sigma(5)-sigma(4))*r1*r2+sigma(6)*(r1*r1-r2*r2)
      ! Computes Hencky stress on the natural Frame
       sigma(1) = stres(1)*lb(1)**2
       sigma(2) = stres(2)*lb(2)**2
       sigma(4) = stres(4)*lb(1)**2
       sigma(5) = stres(5)*lb(2)**2
       IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN
         sigma(3) = stres(3)*(lb(1)**2-lb(2)**2)/LOG(lb(1)/lb(2))/2d0
         sigma(6) = stres(6)*(lb(1)**2-lb(2)**2)/LOG(lb(1)/lb(2))/2d0
       ELSE
         sigma(3) = stres(3)*lb(1)*lb(2)
         sigma(6) = stres(6)*lb(1)*lb(2)
       END IF
      ! Computes Hencky on the Lagrangian Frame
       stres(1) = sigma(1)*r1*r1+sigma(2)*r2*r2-2d0*sigma(3)*r1*r2
       stres(2) = sigma(1)*r2*r2+sigma(2)*r1*r1+2d0*sigma(3)*r1*r2
       stres(3) =(sigma(1)-sigma(2))*r1*r2+sigma(3)*(r1*r1-r2*r2)
       stres(4) = sigma(4)*r1*r1+sigma(5)*r2*r2-2d0*sigma(6)*r1*r2
       stres(5) = sigma(4)*r2*r2+sigma(5)*r1*r1+2d0*sigma(6)*r1*r2
       stres(6) =(sigma(4)-sigma(5))*r1*r2+sigma(6)*(r1*r1-r2*r2)
     ELSE
       stres(1:6) = sigma
     END IF
     WRITE(narch,ERR=9999) stres(1:8)     !print to file
!     WRITE(narch,ERR=9999) stres(1:3),e%gausv(1:3,1)  !stresses and plastic strains
   END DO
 END IF
 !
 ! *** writes complete stress state on file 16 for post-process
 !
 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900,ERR=9999)        !if echo to ASCII file
   osec = -1                              !initializes
   e => heade                             !point to first element
   ielem = 0
   DO                                     !loop over all elements
     IF (.NOT.ASSOCIATED (e) ) EXIT       !if last element processed, exit
     ielem = ielem + 1
     isec = e%matno                      !associated material
     IF( isec /= osec )THEN
       CALL sect_search(isec,found,sec)
       nlayr = sec%iprop(1)                !number of layers
       nvare = sec%iprop(2)                !number of material components
       thic0 = sec%rprop(1)                !original thickness
       IF( sec%secty == 12 )THEN
         mtype = sec%mtbas%mtype             !type of base material
         plast = sec%mtbas%matdef(3) > 1     !plasticity
         IF( mtype == 30 ) THEN
           nvarn = sec%secdef(5)
           IF( nvarn > 0 )THEN
             IF( ALLOCATED(auxil) ) DEALLOCATE(auxil)
             ALLOCATE(auxil(nvarn))
           ELSE
             IF( ALLOCATED(auxil) ) DEALLOCATE(auxil)
             ALLOCATE(auxil(nvare))
           END IF
         END IF
       END IF
       osec = isec
     END IF

     thick = thic0*e%lb        !current thickness

     !-------------------------------------------------------------
     !
     !      WRITE ASCII    Nxx, Nyy, Nxy, Mxx, Myy, Mxy, t/t0
     !      WRITE for GiD  Nxx, Nyy, Nxy, Mxx, Myy, Mxy,
     !                     t/t0, ep1, epN, VMises1, VMisesN
     !
     !-------------------------------------------------------------
     g     = e%stra1(1:3)      !first fundamental form at mid-surface
     sigma = stint(1:6,ielem)      !integrated 2nd Piola Kirchoff stresses
     stres(7:8) = stint(7:8,ielem) !integrated stresses shears
     vm    = stint(9:10,ielem)     !equivalent stresses
    ! ******* consistent modifications follow
     IF( logst .or. plast )THEN
       CALL lgst14(g,r1,r2,lb,'OUTD13',ierr)   !compute eigenvalues and eigenvectors
      ! Computes 2ndP-K stress on the natural Frame
       stres(1) = sigma(1)*r1*r1+sigma(2)*r2*r2+2d0*sigma(3)*r1*r2
       stres(2) = sigma(1)*r2*r2+sigma(2)*r1*r1-2d0*sigma(3)*r1*r2
       stres(3) =(sigma(2)-sigma(1))*r1*r2+sigma(3)*(r1*r1-r2*r2)
       stres(4) = sigma(4)*r1*r1+sigma(5)*r2*r2+2d0*sigma(6)*r1*r2
       stres(5) = sigma(4)*r2*r2+sigma(5)*r1*r1-2d0*sigma(6)*r1*r2
       stres(6) =(sigma(5)-sigma(4))*r1*r2+sigma(6)*(r1*r1-r2*r2)
      ! Computes Hencky stress on the natural Frame
       sigma(1) = stres(1)*lb(1)**2
       sigma(2) = stres(2)*lb(2)**2
       sigma(4) = stres(4)*lb(1)**2
       sigma(5) = stres(5)*lb(2)**2
       IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN
         sigma(3) = stres(3)*(lb(1)**2-lb(2)**2)/LOG(lb(1)/lb(2))/2d0
         sigma(6) = stres(6)*(lb(1)**2-lb(2)**2)/LOG(lb(1)/lb(2))/2d0
       ELSE
         sigma(3) = stres(3)*lb(1)*lb(2)
         sigma(6) = stres(6)*lb(1)*lb(2)
       END IF
      ! Computes Hencky on the Lagrangian Frame
       stres(1) = sigma(1)*r1*r1+sigma(2)*r2*r2-2d0*sigma(3)*r1*r2
       stres(2) = sigma(1)*r2*r2+sigma(2)*r1*r1+2d0*sigma(3)*r1*r2
       stres(3) =(sigma(1)-sigma(2))*r1*r2+sigma(3)*(r1*r1-r2*r2)
       stres(4) = sigma(4)*r1*r1+sigma(5)*r2*r2-2d0*sigma(6)*r1*r2
       stres(5) = sigma(4)*r2*r2+sigma(5)*r1*r1+2d0*sigma(6)*r1*r2
       stres(6) =(sigma(4)-sigma(5))*r1*r2+sigma(6)*(r1*r1-r2*r2)
     ELSE
       stres(1:6) = sigma
     END IF

     !IF( mtype /= 30 )THEN
       strpl = 0d0                   !equivalent plastic strain
       IF( ASSOCIATED (e%gausv)) THEN
         m = 1                         !bottom (first layer)
         DO l=1,2
           strpl(l) = e%gausv(4,m)   !equivalent plastic strain
           m = nlayr                   !top (last layer)
         END DO
       END IF
       WRITE(16,ERR=9999) (stres(l),l=1,8),e%lb,strpl(1:2),vm(1:2)
     !ELSE
     !  CALL ud_sh2_out(strpl,vm,auxil,nvarn,nlayr,sec%mtbas%props, &
     !                  stres,e%gausv,thick)
     !  WRITE(16,ERR=9999) (stres(l),l=1,6),e%lb,strpl(1:2),vm(1:2),auxil(1:nvarn)
     !END IF

    IF(iwrit == 1)WRITE(lures,"(i7,9e12.4)",ERR=9999) e%numel,stres,e%lb,strpl(1:2)

     e => e%next                   !point to next element
   END DO
   ! reinforcements
   IF(iwrit == 1 .AND. ASSOCIATED(rhead) ) WRITE(lures,901,ERR=9999)        !if echo to ASCII file
   rf => rhead                             !point to first element
   DO                                     !loop over all elements
     IF (.NOT.ASSOCIATED (rf) ) EXIT       !if last element processed, exit
     !      WRITE for GiD  N, M, T
     WRITE(16,ERR=9999) rf%stint,0d0,0d0,0d0  !3 null values added
     !      WRITE ASCII    N, M, T
     IF(iwrit == 1)WRITE(lures,"(i7,9e12.4)",ERR=9999) rf%numel,rf%stint

     rf => rf%next                   !point to next element
   END DO
 END IF
 900 FORMAT(5X,'STRESSES FOR SHELL ELEMENTS'//  &
&       '   ELEM   N-XX',9X,'N-YY',9X,'N-XY',   &
&              9X,'M-XX',9X,'M-YY',9X,'M-XY',9X,'Q-X ',9X,'Q-Y ',9X,'Th-rat')

 901 FORMAT(5X,'STRESSES FOR REINFORCED ELEMENTS'//  &
     &         '   ELEM     N',10X,'M',12X,'T')

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outd15
