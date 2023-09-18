 SUBROUTINE crest2(npoin)
 !
 !     Reads contact data from a restar file
 !
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: npoin

 !Local variables
 TYPE (pair2_db), POINTER :: pair
 TYPE (surf2_db), POINTER :: surf
 INTEGER (kind=4) :: ipair,isurf,i,j
 LOGICAL :: area

 READ(51) nsurf,npair,oldis,ctime,wear       !control variables

 IF( wear )THEN
   ALLOCATE( wwear(npoin) )
   READ(51) (wwear(i),i=1,npoin)
 END IF

 !...  Read Pair Data Base

 CALL ini_cont2(headp,tailp)              !initializes pair list

 DO  ipair=1,npair
   CALL new_pair2 (pair)
   READ(51)  pair%pname , pair%master, pair%slave , &
             pair%imast , pair%islav , pair%indcon, pair%ncnod , &
             pair%mtsur , pair%slsur , pair%mtdof , pair%sldof , &
             pair%freq  , &
             pair%npenal, pair%tpenal, pair%static, pair%kinet , &
             pair%cutoff, pair%gapinc, pair%hcont , pair%econt , &
             pair%hv    , pair%rec   , pair%start , &
             pair%end   , &
             pair%prev  , pair%press , pair%cpress, pair%auto  , &
             pair%wrink
   ALLOCATE ( pair%issdb(nisdb,pair%ncnod) )
   ALLOCATE ( pair%rssdb(nrsdb,pair%ncnod) )
   READ(51) ((pair%issdb(i,j),i=1,nisdb),j=1,pair%ncnod)
   READ(51) ((pair%rssdb(i,j),i=1,nrsdb),j=1,pair%ncnod)
   IF(pair%press) THEN
     ALLOCATE( pair%presn(pair%ncnod) )
     READ(51) (pair%presn(j),j=1,pair%ncnod)
   END IF
   IF(pair%wrink) THEN
     ALLOCATE( pair%mingp(pair%ncnod) )
     READ(51) (pair%mingp(j),j=1,pair%ncnod)
   END IF
   CALL add_pair2 (pair, headp, tailp)
 END DO

 !...  Read Surface Database

 CALL ini_srf2 (shead,stail)              !initializes surface list

 DO isurf=1,nsurf
   CALL new_surf2 (surf)
   READ(51)  surf%sname, surf%cxc,    surf%bottom, surf%auto,   &
             surf%ncnod, surf%curved, surf%nsegm,  surf%iwrit,  &
             surf%imcod, surf%iscod,  surf%press,  surf%density, area
   IF( surf%iscod )THEN
     ALLOCATE ( surf%lcnod(surf%ncnod) )
     READ(51) (surf%lcnod(i),i=1,surf%ncnod)
   END IF
   IF(surf%nsegm > 0)THEN
     ALLOCATE ( surf%lcseg(2,surf%nsegm) )
     READ(51) ((surf%lcseg(i,j),i=1,2),j=1,surf%nsegm)
   END IF
   IF( surf%imcod )THEN
     ALLOCATE ( surf%nhseg(2,surf%nsegm) , surf%xc(2,surf%nsegm) )
     READ(51) ((surf%nhseg(i,j),i=1,2),j=1,surf%nsegm)
     READ(51) ((surf%xc(i,j),i=1,2),j=1,surf%nsegm)
   END IF
   IF(surf%bottom )THEN
     ALLOCATE ( surf%lcseb(2,surf%nsegm) , surf%nhseb(2,surf%nsegm) )
     READ(51) ((surf%lcseb(i,j),i=1,2),j=1,surf%nsegm)
     READ(51) ((surf%nhseb(i,j),i=1,2),j=1,surf%nsegm)
   END IF

   IF((surf%curved .OR. surf%auto) .AND. surf%iscod ) THEN
     ALLOCATE(surf%tn(2,surf%ncnod))
     READ(51) ((surf%tn(i,j),i=1,2),j=1,surf%ncnod)
   END IF
!   IF(surf%curved .OR. surf%auto ) THEN
!     ALLOCATE(surf%nr(npoin))
!     READ(51) (surf%nr(j),j=1,npoin)
!   END IF
   IF( surf%curved  .AND. surf%imcod ) THEN
     ALLOCATE(surf%cu(surf%nsegm))
     READ(51) (surf%cu(j),j=1,surf%nsegm)
   END IF
   IF( surf%auto  ) THEN
     ALLOCATE(surf%tc(2,surf%nsegm))
     READ(51) ((surf%tc(i,j),i=1,2),j=1,surf%nsegm)
   END IF
   IF( area )  THEN
     ALLOCATE(surf%area(surf%ncnod))
     READ(51)  (surf%area(i),i=1,surf%ncnod)
   END IF
   CALL add2_srf (surf, shead, stail)

 END DO

 RETURN

 END SUBROUTINE crest2
