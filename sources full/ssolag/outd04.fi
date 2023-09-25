 SUBROUTINE outd04(b1,b2,iwrit,heade,nreqs,narch,ngrqs,ttime,ngaus,isg,psg)
 !--------------------------------------------------------------------
 !
 !       output routine for 3-D solid element (TLF)
 !
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: iwrit,    &!flag for output into ASCII file
                                  nreqs,    &!number of elements for output
                                  ngaus,    &!number of integration points
                                  narch,isg(2,2)      !output file number
 INTEGER (kind=4), POINTER :: ngrqs(:)   !element numbers for output
 REAL (kind=8), INTENT(IN) :: ttime,psg(2,2)
 LOGICAL, INTENT(IN) :: b1,   & !flag for output at selected points
                        b2      !flag for global output
 TYPE (ele04), POINTER  :: heade  !pointer to first element

 INTEGER (kind=4) isec,k,l,m,ielem,osec,i,j
 REAL (kind=8) stres(6),strpl,sxx,syy,sxy,szz,sxz,syz,vm,preys,var(8,2)
 TYPE (ele04), POINTER :: e
 TYPE (section), POINTER :: sec
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
     IF( yield .AND. ASSOCIATED (e%gausv)  )THEN
       strpl = e%gausv( 7,k) !effective plastic strain of first component
     ELSE
       strpl = 0d0
     END IF
     preys = 0d0 !actual yield stress ??
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
 ! *** writes complete stress state on file 04 for post-process
 !
 IF(b2) THEN
   IF(iwrit == 1) WRITE(lures,900,ERR=9999)  !if echo to ASCII file
   osec = -1                              !initializes
   e => heade                             !point to first element
   DO                                     !loop over all elements
     IF (.NOT.ASSOCIATED (e) ) EXIT       !if last element processed, exit
     isec = e%matno                      !associated material
     newm = isec /= osec
     IF( newm )THEN
       CALL sect_search(isec,found,sec)
       yield = sec%mtbas%matdef(3) > 1
       osec = isec
     END IF
     !-------------------------------------------------------------
      IF(iwrit == 1)WRITE(lures,"('elem = ',i5,' alpha= ',e12.4)",ERR=9999) e%numel,e%alpha
     var = 0d0
     DO k=1,ngaus
       stres = e%stint(:,k)     !integrated Kirchhoff stresses on material axis
       strpl = 0d0                   !equivalent plastic strain
       vm    = 0d0                   !von Mises equivalent stress
       IF( yield .AND. ASSOCIATED(e%gausv) )  strpl = e%gausv(7,k)   !equivalent plastic strain
       vm = (stres(1)+stres(2)+stres(3))/3d0
       sxx = stres(1)-vm         !
       syy = stres(2)-vm         !
       szz = stres(3)-vm         !
       sxy = stres(4)            !Stress XY
       sxz = stres(5)            !Stress XZ
       syz = stres(6)            !Stress YZ
       vm = SQRT(1.5d0*(sxx**2+syy**2+szz**2+2d0*sxy**2+2d0*sxz**2+2d0*syz**2))  !Von Mises stress
       IF(iwrit == 1)WRITE(lures,920,ERR=9999) k,stres(1:6),strpl,vm
       DO i=1,2
         DO j=1,2
           IF( k == isg(j,i) )var(:,i) = var(:,i)+psg(j,i)*(/stres(1:6),strpl,vm/)
         END DO
       END DO
     END DO
     DO i=1,2
       WRITE(16,ERR=9999) var(:,i)
       IF(iwrit == 1)WRITE(lures,920,ERR=9999) i,var(:,i)
     END DO
     e => e%next                   !point to next element
   END DO
 END IF

 RETURN
 900 FORMAT(5X,'STRESSES FOR 3-D SOLID ELEMENTS in MATERIAL AXES'//   &
 &       '      G   S-XX',9X,'S-YY',9X,'S-ZZ',9X,'S-XY',9X,'S-XZ',9X,'S-YZ', &
 &       9X,'EPST',9X,'VMST')
 920 FORMAT(5x,i2,8e13.4)
 9999 CALL runen2('')
 END SUBROUTINE outd04
