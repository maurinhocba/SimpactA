 SUBROUTINE outd18(b1,b2,iwrit,heade,nreqs,narch,ngrqs,ttime,ngaus)
 !--------------------------------------------------------------------
 !
 !       output routine for 3-D solid element (TLF)
 !
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: iwrit,    &!flag for output into ASCII file
                                  nreqs,    &!number of elements for output
                                  ngaus,    &!number of integration points
                                  narch      !output file number
 INTEGER (kind=4), POINTER :: ngrqs(:)   !element numbers for output
 REAL (kind=8), INTENT(IN) :: ttime
 LOGICAL, INTENT(IN) :: b1,   & !flag for output at selected points
                        b2      !flag for global output
 TYPE (ele18), POINTER  :: heade  !pointer to first element

 INTEGER (kind=4) isec,k,l,m,ielem,osec,nlayr,il,nvarn
 INTEGER (kind=4) nvarsm,nvarsf,ivm,ivf
 REAL (kind=8) stres(6),strpl,sxx,syy,sxy,szz,sxz,syz,vm,preys
 REAL (kind=8) dlay,dmat,dfib,thlay,angly,add
 TYPE (ele18), POINTER :: e
 TYPE (section), POINTER :: sec
 TYPE (mater), POINTER :: mat
 INTEGER (kind=4), ALLOCATABLE :: lvari(:), lmats(:)
 REAL (kind=8), ALLOCATABLE :: auxil(:)
 LOGICAL :: found,newm,yield

 !     output at selected elements for HISTORY

 IF(b1 .AND. (nreqs > 0)) THEN
   osec = -1
   WRITE(narch,ERR=9999) ttime                  !write process time
   DO l=1,nreqs                        !for each element
     k = MOD(ngrqs(l)-1,ngaus) + 1
     ielem = (ngrqs(l)-1) / ngaus + 1
     e => heade                   !point to first element
     DO                           !search for element with label IELEM
       IF( e%numel == ielem )EXIT !if found, exit search
       e => e%next                !point to next element
     END DO
     isec = e%matno
     newm = isec /= osec
     IF( newm )THEN
       CALL sect_search(isec,found,sec)
       yield = sec%mtbas%matdef(3) > 1
       osec = isec
     END IF
     stres = e%stint(:,k)     !integrated Kirchhoff stresses
     IF( yield )THEN
       IF( ASSOCIATED(e%gausv) )THEN
         strpl = e%gausv( 7,k) !plastic strain of first component
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
     WRITE(narch,ERR=9999)(stres(m),m=1,6),strpl,vm,preys
     !WRITE(narch,ERR=9999)e%gausv(1:4,k),stres(1:4),strpl
   END DO
 END IF
 !
 ! *** writes complete stress state on file 16 for post-process
 !
 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900,ERR=9999)  !if echo to ASCII file
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
       IF (nlayr > 0) THEN             !Layered solid
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
     DO k=1,ngaus
       stres = e%stint(:,k)     !integrated Kirchhoff stresses on material axis
       strpl = 0d0                    !equivalent plastic strain (or damage parameter)
       vm    = 0d0                    !von Mises equivalent stress
       IF(nvarn>0) auxil(1:nvarn)=0d0 !initialization of auxiliar variables

       IF( ASSOCIATED( e%gausv ) )THEN
         IF (nlayr == 0) THEN          ! Non laminated solid
             mat => sec%mtbas
             SELECT CASE (mat%mtype)
             CASE (7)
               dlay = e%gausv(2,k)                    !OK FR
             CASE (92)
               nvarsm = mat%matdef(6)  
               nvarsf = mat%matdef(7)  !decommented by FF
               ivm=24
               ivf=ivm+nvarsm
               dmat = e%gausv(ivm+2,k)                !OK FR
               IF( nvarsf > 1) THEN             !added by FF
                 dfib = e%gausv(ivf+2,k)              !OK FR
               ELSE            
                 dfib = 0d0
               END IF
               dlay = (mat%prope(6) * dmat) +    &
                      (mat%prope(7) * dfib)
               IF ((dmat /= 0.0d0) .OR. (dfib /= 0.0d0)) THEN
                 add = 1.0
               ELSE
                 add = 0.0
               END IF
             CASE DEFAULT
               dlay = 0d0
               yield = mat%matdef(3) > 1
               IF(yield) strpl = e%gausv(7,k)   !equivalent plastic strain
             END SELECT
             IF( nvarn >= 4 )THEN
               auxil(1) = dlay     !Total damage
               auxil(2) = dlay     !Maximum damage in laminate
!               auxil(6) = 1        !Layer id for max damage
!               auxil(7) = 0d0      !Layer orientation for max damage
               auxil(3) = add      !Qty of layers damaged in laminate
               auxil(4) = dmat     !Matrix damage
               auxil(5) = dfib     !Fiber damage
               strpl = auxil(1)    !Total damage (provisoriamente!!!)
             END IF

         ELSE  ! laminated solid
           ! Gets damage internal variable from each layer
           ! and evaluates its average to post-process
           DO il = 1, nlayr
             mat => pmats(lmats(il))%p
             thlay = sec%rprop(41+il)
             angly = sec%rprop(41+nlayr+il) !layer angle
             SELECT CASE (mat%mtype)
             CASE (7)
               ivm=lvari(il)-1
               dlay = e%gausv(ivm+2,k)                    !OK FR
             CASE (92)
               nvarsm = mat%matdef(6)  
               nvarsf = mat%matdef(7)                     !decommented by FF
               ivm=lvari(il)-1+24
               ivf=ivm+nvarsm
               dmat = e%gausv(ivm+2,k)                    !OK FR
               IF( nvarsf > 1 )THEN
                 dfib = e%gausv(ivf+2,k)                    !OK FR
               ELSE
                 dfib = 0d0
               END IF
               dlay = (mat%prope(6) * dmat) +    &
                      (mat%prope(7) * dfib)
               IF ((dmat /= 0.0d0) .OR. (dfib /= 0.0d0)) THEN
                 add = 1.0
               ELSE
                 add = 0.0
               ENDIF
             CASE DEFAULT
               ivm=lvari(il)-1
               dlay = 0d0
               yield = mat%matdef(3) > 1
               IF(yield) strpl = strpl + thlay*e%gausv(ivm+7,k)   !equivalent plastic strain
             END SELECT
             IF( nvarn >= 4 )THEN
               auxil(1) = auxil(1) + thlay*dlay !layer contribution to damage index 
               IF ( dlay > auxil(2) ) THEN
                 auxil(2) = dlay                !Maximum damage in laminate
!                 auxil(6) = il                  !Layer id for max damage
!                 auxil(7) = angly               !Layer orientation for max damage
               END IF
               auxil(3) = auxil(3) + add        !Qty of layers damaged in laminate
               auxil(4) = auxil(4) + thlay*dmat !Average of matrix damage
               auxil(5) = auxil(5) + thlay*dfib !Average of fiber damage
             END IF
           END DO !nlayr
           strpl = strpl/sec%rprop(1)       !Average equivalent plastic strain
           IF( nvarn >= 4 )THEN
             auxil(1) = auxil(1)/sec%rprop(1) ! Average total damage
             auxil(4) = auxil(4)/sec%rprop(1) ! Average matrix damage
             auxil(5) = auxil(5)/sec%rprop(1) ! Average fiber damage
             strpl = auxil(1)                 ! Average total damage (provisoriamente!!!)
           END IF
         END IF

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
         WRITE(16,ERR=9999) (stres(l),l=1,6),strpl,vm,(auxil(l),l=1,nvarn)
         IF(iwrit == 1)WRITE(lures,920,ERR=9999) e%numel,k,stres(1:6),strpl,vm,(auxil(l),l=1,nvarn)
       ELSE
         WRITE(16,ERR=9999) (stres(l),l=1,6),strpl,vm
         IF(iwrit == 1)WRITE(lures,920,ERR=9999) e%numel,k,stres(1:6),strpl,vm
       END IF
     END DO

     e => e%next                   !point to next element
   END DO
 END IF

 RETURN
 900 FORMAT(5X,'STRESSES FOR 3-D SOLID ELEMENTS in MATERIAL AXES'//   &
 &       ' ELEM G   S-XX',9X,'S-YY',9X,'S-ZZ',9X,'S-XY',9X,'S-XZ',9X,'S-YZ', &
 &       9X,'EPST',9X,'VMST')
 920 FORMAT(i5,i2,13e13.4)
 9999 CALL runen2('')
 END SUBROUTINE outd18
