 SUBROUTINE cdump1( )
 !
 !     writes contact data for restar
 !
 !USE cont1_db   !INTENT(IN)
 !USE c_input,ONLY: runen2
 IMPLICIT NONE

 !local variables
 TYPE (pair1_db), POINTER :: pair
 TYPE (surf1_db), POINTER :: surf
 INTEGER (kind=4) :: ipair,isurf,i,j


 WRITE(50,ERR=9999) nsurf,npair,ctime     !Control variables

 !     Dumps pair data base

 pair => headp
 DO  ipair=1,npair
   WRITE(50,ERR=9999) pair%pname , pair%master, pair%slave ,              &
                      pair%imast , pair%islav , pair%indcon, pair%ncnod , &
                      pair%npenal, pair%tpenal, pair%static, pair%kinet,  &
                      pair%start,  pair%end,    pair%prev,   pair%press

   WRITE(50,ERR=9999) ((pair%issdb(i,j),i=1,nisdb),j=1,pair%ncnod)
   WRITE(50,ERR=9999) ((pair%rssdb(i,j),i=1,nrsdb),j=1,pair%ncnod)
   IF (pair%press) WRITE(50,ERR=9999) (pair%presn(j),j=1,pair%ncnod)

   pair => pair%next
 END DO

 !     Dumps surface data base

 surf => shead
 DO isurf=1,nsurf
   WRITE(50,ERR=9999) surf%sname, surf%ncnod, surf%nsegm, surf%iwrit,  &
                      surf%imcod, surf%iscod, surf%density, surf%diam
   IF( surf%iscod )WRITE(50,ERR=9999) (surf%lcnod(i),i=1,surf%ncnod)
   IF( surf%nsegm > 0 )WRITE(50,ERR=9999) ((surf%lcseg(i,j),i=1,2),j=1,surf%nsegm)

   surf => surf%next
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE cdump1
