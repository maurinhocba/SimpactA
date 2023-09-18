 SUBROUTINE prdbln(flag)
 !.... write drawbead lines for visualization
 !.... open file for forces with surfaces

  USE ctrl_db, ONLY: ndime
  USE c_input, ONLY : openfi
  IMPLICIT NONE
  !   d u m m y  v a r i a b l e s
   INTEGER(kind=4) :: flag
  !   l o c a l   v a r i a b l e s

  INTEGER ::  i,j,nn,idbea,ns
  TYPE (db3pair_db), POINTER :: pair

   !.... Open file por post-process


  nn = 6        !position of IDBFOUT code
  numtdf = 0

  pair => headp

  DO i = 1, ndbea
     IF (pair%iparm(nn)==1) numtdf = numtdf + 1
     pair => pair%next
  END DO

  IF (numtdf > 0) THEN
     CALL openfi(42,flag=flag)
     IF (flag==0) THEN
       WRITE (42,ERR=9999) numtdf, ndime
       pair => headp
       DO i = 1, ndbea
         IF (pair%iparm(nn)==1)  WRITE (42,ERR=9999) pair%pname
         pair => pair%next
       END DO
     END IF

     !.... WRITE sheet for post-proccessing
     CALL openfi (43,flag=flag)
     IF( flag /= 1 )THEN
       WRITE (43,ERR=9999) ndbea
       pair => headp
       DO idbea=1,ndbea
         nn = pair%iparm(5)
         ns = pair%iparm(2)
         WRITE (43,ERR=9999) pair%pname,pair%sname,nn,ns,pair%iparm(8) 
         WRITE (43,ERR=9999) pair%irect(1:nn,1:ns)
         pair => pair%next
       END DO
     END IF

     !.... WRITE Drawb lines for post-proccessing
     nn = 2
     IF( flag /= 1)THEN
       CALL openfi(18,flag=1)        !open file
       pair => headp
       DO idbea=1,ndbea
         ns =  pair%iparm(1)-1
         WRITE(18,err=9999) ns ,pair%pname !NSegs Drawbead name
         WRITE(18,err=9999) nn             !Nodes by element and Number of pair
         WRITE(18,err=9999) (pair%lnd(1:nn,j),j=1,ns)   !Drawbeads conectivities
         pair => pair%next
       END DO
       CLOSE (18)                      !close file
     END IF

  END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE prdbln
