SUBROUTINE rsdbe3 (pair, itask)
  ! dumping / restar of a pair
  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  CHARACTER (LEN=*), INTENT (IN) :: itask
  TYPE (db3pair_db) :: pair

  INTEGER(kind=4) :: nrtm, nmn, mst, ncrosi, ndbpoin
  INTEGER (kind=4) :: i,j
  !-----------------------------------------------

  IF (TRIM(itask) == 'DUMPIN') THEN
     WRITE (50,ERR=9999) pair%iparm,pair%rparm,pair%pname,pair%sname
     ndbpoin = pair%iparm (1)
     nrtm = pair%iparm (2)
     nmn = pair%iparm (3)
     mst = pair%iparm (4)
     ncrosi = pair%iparm (7)
     WRITE (50,ERR=9999) (pair%nstnew(i),i=1,nrtm)
     WRITE (50,ERR=9999) (pair%nstold(i),i=1,nrtm)
     WRITE (50,ERR=9999) ((pair%ldbs(i,j),i=1,2),j=1,nrtm)
     WRITE (50,ERR=9999) ((pair%irect(i,j),i=1,4),j=1,nrtm)
     WRITE (50,ERR=9999) ((pair%fdat(i,j),i=1,3),j=1,nrtm)
     WRITE (50,err=9999) ((pair%lnd(i,j),i=1,2),j=1,ndbpoin-1)
     WRITE (50,err=9999) (pair%daz(i),i=1,nrtm)
     WRITE (50,ERR=9999) ((pair%xl(i,j),i=1,2),j=1,ndbpoin)
     WRITE (50,ERR=9999) (pair%nseg(i),i=1,nmn+1)
     WRITE (50,ERR=9999) (pair%msr(i),i=1,nmn)
     WRITE (50,ERR=9999) (pair%lmsr(i),i=1,mst)
     WRITE (50,ERR=9999) (pair%nstini(i),i=1,MAX(1,ncrosi))
     WRITE (50,ERR=9999) ((pair%fdaini(i,j),i=1,2),j=1,MAX(1,ncrosi))

  ELSE IF (TRIM(itask) == 'RESTAR') THEN

     READ (51) pair%iparm,pair%rparm,pair%pname,pair%sname
     ndbpoin = pair%iparm (1)
     nrtm = pair%iparm (2)
     nmn = pair%iparm (3)
     mst = pair%iparm (4)
     ncrosi = pair%iparm (7)
     DO i = 1,5
       CALL alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi ,i)
     END DO
     READ (51) (pair%nstnew(i),i=1,nrtm)
     READ (51) (pair%nstold(i),i=1,nrtm)
     READ (51) ((pair%ldbs(i,j),i=1,2),j=1,nrtm)
     READ (51) ((pair%irect(i,j),i=1,4),j=1,nrtm)
     READ (51) ((pair%fdat(i,j),i=1,3),j=1,nrtm)
     READ (51) ((pair%lnd(i,j),i=1,2),j=1,ndbpoin-1)
     READ (51) (pair%daz(i),i=1,nrtm)
     READ (51) ((pair%xl(i,j),i=1,2),j=1,ndbpoin)
     READ (51) (pair%nseg(i),i=1,nmn+1)
     READ (51) (pair%msr(i),i=1,nmn)
     READ (51) (pair%lmsr(i),i=1,mst)
     READ (51) (pair%nstini(i),i=1,MAX(1,ncrosi))
     READ (51) ((pair%fdaini(i,j),i=1,2),j=1,MAX(1,ncrosi))

  END IF

RETURN
 9999 CALL runen2('')
END SUBROUTINE rsdbe3
