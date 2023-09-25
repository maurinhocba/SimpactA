 SUBROUTINE cdump2(npoin)
 !
 !     writes contact data for restar
 !
 !USE cont2_db   !INTENT(IN)
 !USE c_input,ONLY: runen2
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: npoin

 !local variables
 TYPE (pair2_db), POINTER :: pair
 TYPE (surf2_db), POINTER :: surf
 INTEGER (kind=4) :: ipair,isurf,i,j
 LOGICAL :: area


 WRITE(50,ERR=9999) nsurf,npair,oldis,ctime,wear   !Control variables

 IF (wear) WRITE(50,ERR=9999) (wwear(i),i=1,npoin)
 !     Dumps pair data base

 pair => headp
 DO  ipair=1,npair
   WRITE(50,ERR=9999) pair%pname , pair%master, pair%slave , &
                      pair%imast , pair%islav , pair%indcon, pair%ncnod , &
                      pair%mtsur , pair%slsur , pair%mtdof , pair%sldof , &
                      pair%freq  , &
                      pair%npenal, pair%tpenal, pair%static, pair%kinet , &
                      pair%cutoff, pair%gapinc, pair%hcont , pair%econt , &
                      pair%hv    , pair%rec   , pair%start , &
                      pair%end   , &
                      pair%prev  , pair%press , pair%cpress, pair%auto  , &
                      pair%wrink

   WRITE(50,ERR=9999) ((pair%issdb(i,j),i=1,nisdb),j=1,pair%ncnod)
   WRITE(50,ERR=9999) ((pair%rssdb(i,j),i=1,nrsdb),j=1,pair%ncnod)
   IF (pair%press) WRITE(50,ERR=9999) (pair%presn(j),j=1,pair%ncnod)
   IF (pair%wrink) WRITE(50,ERR=9999) (pair%mingp(j),j=1,pair%ncnod)

   pair => pair%next
 END DO

 !     Dumps surface data base

 surf => shead
 DO isurf=1,nsurf
   area =  therm .AND. ASSOCIATED( surf%area)
   WRITE(50,ERR=9999) surf%sname, surf%cxc,    surf%bottom, surf%auto,   &
                      surf%ncnod, surf%curved, surf%nsegm,  surf%iwrit,  &
                      surf%imcod, surf%iscod,  surf%press,  surf%density, area
   IF( surf%iscod )WRITE(50,ERR=9999) (surf%lcnod(i),i=1,surf%ncnod)
   IF( surf%nsegm > 0 )WRITE(50,ERR=9999) ((surf%lcseg(i,j),i=1,2),j=1,surf%nsegm)
   IF( surf%imcod )THEN
     WRITE(50,ERR=9999) ((surf%nhseg(i,j),i=1,2),j=1,surf%nsegm)
     WRITE(50,ERR=9999) ((surf%xc(i,j),i=1,2),j=1,surf%nsegm)
   END IF
   IF( surf%bottom )THEN
     WRITE(50,ERR=9999) ((surf%lcseb(i,j),i=1,2),j=1,surf%nsegm)
     WRITE(50,ERR=9999) ((surf%nhseb(i,j),i=1,2),j=1,surf%nsegm)
   END IF
   IF((surf%curved .OR. surf%auto) .AND. surf%iscod ) &
     WRITE(50,ERR=9999) ((surf%tn(i,j),i=1,2),j=1,surf%ncnod)
!   IF( surf%curved .OR. surf%auto ) &
!     WRITE(50,ERR=9999) (surf%nr(j),j=1,npoin)
   IF( surf%curved  .AND. surf%imcod ) &
     WRITE(50,ERR=9999) (surf%cu(j),j=1,surf%nsegm)
   IF( surf%auto  ) &
     WRITE(50,ERR=9999) ((surf%tc(i,j),i=1,2),j=1,surf%ncnod)
   IF( area )  WRITE(50,ERR=9999) (surf%area(i),i=1,surf%ncnod)

   surf => surf%next
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE cdump2
