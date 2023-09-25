 SUBROUTINE outd04(b1,b2,iwrit,heade,nreqs,narch,ngrqs,ttime,nodvo)
 !--------------------------------------------------------------------
 !
 !       output routine for 3-D solid element (TLF)
 !
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: iwrit,    &!flag for output into ASCII file
                                  nreqs,    &!number of elements for output
                                  narch      !output file number
 INTEGER (kind=4), POINTER :: ngrqs(:)   !element numbers for output
 REAL (kind=8), INTENT(IN) :: ttime
 LOGICAL, INTENT(IN) :: b1,   & !flag for output at selected points
                        b2,   & !flag for global output
                        nodvo   !flag for nodal integration
 TYPE (ele04), POINTER  :: heade  !pointer to first element

 INTEGER (kind=4) isec,l,m,ielem,osec,nlayr,il,nvarn,i,n
 REAL (kind=8) stres(6),strpl,sxx,syy,sxy,szz,sxz,syz,vm,preys,damal
 TYPE (ele04), POINTER :: e
 TYPE (section), POINTER :: sec
 TYPE (mater), POINTER :: mat
 INTEGER (kind=4), ALLOCATABLE :: lvari(:), lmats(:)
 REAL (kind=8), ALLOCATABLE :: auxil(:)
 LOGICAL :: newm,yield

 !     output at selected elements for HISTORY

 IF(b1 .AND. (nreqs > 0)) THEN
   osec = -1
   WRITE(narch,ERR=9999) ttime                  !write process time
   DO l=1,nreqs                        !for each element
     ielem = ngrqs(l)
     e => heade                   !point to first element
     DO                           !search for element with label IELEM
       IF( e%numel == ielem )EXIT !if found, exit search
       e => e%next                !point to next element
     END DO
     isec = e%matno
     newm = isec /= osec
     IF( newm )THEN
       sec => psecs(isec)%p
       yield = sec%mtbas%matdef(3) > 1
       osec = isec
     END IF
     stres = e%stint(1:6)     !integrated Kirchhoff stresses
     IF( yield )THEN
       IF( ASSOCIATED(e%gausv) )THEN
         strpl = e%gausv( 7)   !plastic strain of first component
       ELSE
         strpl = 0d0           !plastic strain of first component
       END IF
       preys = 0d0 !actual yield stress ??
     ELSE
       strpl = 0d0
       preys = 0d0 !actual yield stress
     END IF
     vm = (stres(1)+stres(2)+stres(3))/3d0
     sxx = stres(1)-vm         !
     syy = stres(2)-vm         !
     szz = stres(3)-vm         !
     sxy = stres(4)            !Stress XY
     sxz = stres(5)            !Stress XZ
     syz = stres(6)            !Stress YZ
     vm = SQRT(1.5d0*(sxx**2+syy**2+szz**2+2d0*sxy**2+2d0*sxz**2+2d0*syz**2))  !Von Mises stress
     WRITE(narch)(stres(m),m=1,6),strpl,vm,preys
     !WRITE(narch)e%gausv(1:4),stres(1:4),strpl
     WRITE(narch,ERR=9999)(stres(m),m=1,6),strpl,vm,preys
     !WRITE(narch,ERR=9999)e%gausv(1:4),stres(1:4),strpl
   END DO
 END IF
 !
 ! *** writes complete stress state on file 16 for post-process
 !
 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900)        !if echo to ASCII file
   osec = -1                              !initializes
   e => heade                             !point to first element
   DO                                     !loop over all elements
     IF (.NOT.ASSOCIATED (e) ) EXIT       !if last element processed, exit
     isec = e%matno                      !associated material
     ielem = e%numel
     newm = isec /= osec
     IF( newm )THEN
       sec => psecs(isec)%p
       nlayr = sec%iprop(1)
       nvarn = sec%secdef(5)
       IF( ALLOCATED(auxil) ) DEALLOCATE(auxil)
       IF( nvarn > 0 )   ALLOCATE(auxil(nvarn))
       IF (nlayr /= 0) THEN
         !mtype =  0                   !type of base material
         IF( ALLOCATED(lvari) ) DEALLOCATE(lvari,lmats)
         ALLOCATE(lvari(nlayr),lmats(nlayr))
         lmats(1:nlayr)= sec%iprop(      5:  nlayr+4)    !pointers to layer material
         lvari(1:nlayr)= sec%iprop(nlayr+5:2*nlayr+4)    !pointers to layer internal variables
       ELSE
         yield = sec%mtbas%matdef(3) > 1
       END IF
       osec = isec
     END IF
     !-------------------------------------------------------------

     stres = e%stint(1:6)          !Kirchhoff stresses on material axis
     strpl = 0                     !equivalent plastic strain
     vm    = 0d0                   !von Mises equivalent stress
     IF( nvarn > 0 )   auxil(1:nvarn) = 0d0

     IF (nlayr == 0) THEN          ! Non laminated solid
       IF( nodvo )THEN
         strpl = 0d0 ! initializes equivalent plastic strain
         DO i=1,nnode
           n = e%lnods(i)
           IF( ASSOCIATED(ngauv(n)%p )) strpl = strpl + ngauv(n)%p(7)
         END DO
         strpl = strpl/4d0
       ELSE
         IF( ASSOCIATED( e%gausv ) )THEN
           IF (sec%mtbas%mtype == 7) THEN
             strpl = e%gausv(2)    ! Damage parameter
           ELSE IF(sec%mtbas%mtype == 92) then
             strpl = e%gausv(25)   ! Damage parameter
           ELSE
             IF( yield ) strpl = e%gausv(7)   !equivalent plastic strain
           END IF
         END IF
       END IF

     ELSE  ! laminated solid
       ! Gets damage internal variable from each layer
       ! and evaluates its average to post-process
       strpl = 0.0
       DO il = 1, nlayr
         IF( ASSOCIATED( e%gausv ) )THEN
             mat => pmats(lmats(il))%p
             SELECT CASE (mat%mtype)
             CASE (7)
               damal= e%gausv(lvari(il)+1)
             CASE (92)
               damal= e%gausv(lvari(il)+24)
             CASE DEFAULT
               damal= 0d0
               yield = mat%matdef(3) > 1
               IF (yield) strpl = strpl + e%gausv(7)   !equivalent plastic strain
             END SELECT
             IF( nvarn >= 2 )THEN
               auxil(1) = auxil(1) + damal/nlayr !layer contribution to damage index
               auxil(2) = MAX( auxil(2), damal ) !Maximum damage in laminate
             END IF
         END IF
       END DO
       strpl = strpl/nlayr !Average equivalent plastic strain
       IF(nvarn>0) strpl = auxil(1)    !Average damage (provisoriamente!!!)
     END IF

     vm = (stres(1)+stres(2)+stres(3))/3d0
     sxx = stres(1)-vm         !
     syy = stres(2)-vm         !
     szz = stres(3)-vm         !
     sxy = stres(4)            !Stress XY
     sxz = stres(5)            !Stress XZ
     syz = stres(6)            !Stress YZ
     vm = SQRT(1.5d0*(sxx**2+syy**2+szz**2+2d0*sxy**2+2d0*sxz**2+2d0*syz**2))  !Von Mises stress
     IF( nvarn > 0 )THEN
     !  auxil(1:nvarn) = 0d0 !              !INCLUDE HERE THE ADDITIONAL VARIABLES !!! FR!!!
       WRITE(16,ERR=9999) (stres(l),l=1,6),strpl,vm,(auxil(l),l=1,nvarn)
       IF(iwrit == 1)WRITE(lures,920,ERR=9999) e%numel,stres(1:6),strpl,vm,auxil(1:nvarn)
     ELSE
       WRITE(16,ERR=9999) (stres(l),l=1,6),strpl,vm
       IF(iwrit == 1)WRITE(lures,920,ERR=9999) e%numel,stres(1:6),strpl,vm
     END IF

     e => e%next                   !point to next element
   END DO

 END IF

 900 FORMAT(5X,'STRESSES FOR 3-D SOLID ELEMENTS in MATERIAL AXES'//   &
 &       ' ELEM   S-XX',9X,'S-YY',9X,'S-ZZ',9X,'S-XY',9X,'S-XZ',9X,'S-YZ', &
 &       9X,'EPST',9X,'VMST')
 920 FORMAT(i5,8e13.4)
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE outd04
