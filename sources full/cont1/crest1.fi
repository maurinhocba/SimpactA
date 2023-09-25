 SUBROUTINE crest1()
 !
 !     Reads contact data from a restar file
 !
 IMPLICIT NONE
 !Local variables
 TYPE (pair1_db), POINTER :: pair
 TYPE (surf1_db), POINTER :: surf
 INTEGER (kind=4) :: ipair,isurf,i,j

 READ(51) nsurf,npair,ctime     !control variables

 !...  Read Pair Data Base

 CALL ini_cont1(headp,tailp)              !initializes pair list

 DO  ipair=1,npair
   CALL new_pair1 (pair)
   READ(51)  pair%pname , pair%master, pair%slave ,              &
             pair%imast , pair%islav , pair%indcon, pair%ncnod , &
             pair%npenal, pair%tpenal, pair%static, pair%kinet,  &
             pair%start,  pair%end,    pair%prev,   pair%press
   ALLOCATE ( pair%issdb(nisdb,pair%ncnod) )
   ALLOCATE ( pair%rssdb(nrsdb,pair%ncnod) )
   READ(51) ((pair%issdb(i,j),i=1,nisdb),j=1,pair%ncnod)
   READ(51) ((pair%rssdb(i,j),i=1,nrsdb),j=1,pair%ncnod)
   IF(pair%press) THEN
     ALLOCATE( pair%presn(pair%ncnod) )
     READ(51) (pair%presn(j),j=1,pair%ncnod)
   END IF
   CALL add_pair1 (pair, headp, tailp)
 END DO

 !...  Read Surface Database

 CALL ini_srf1 (shead,stail)              !initializes surface list

 DO isurf=1,nsurf
   CALL new_surf1 (surf)
   READ(51)  surf%sname, surf%ncnod, surf%nsegm, surf%iwrit,  &
             surf%imcod, surf%iscod, surf%density, surf%diam
   IF( surf%iscod )THEN
     ALLOCATE ( surf%lcnod(surf%ncnod) )
     READ(51) (surf%lcnod(i),i=1,surf%ncnod)
   END IF
   IF(surf%nsegm > 0)THEN
     ALLOCATE ( surf%lcseg(2,surf%nsegm) )
     READ(51) ((surf%lcseg(i,j),i=1,2),j=1,surf%nsegm)
   END IF
   CALL add1_srf (surf, shead, stail)

 END DO

 RETURN

 END SUBROUTINE crest1
