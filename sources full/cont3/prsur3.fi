 SUBROUTINE prsur3(flag)
 !.... write rigid surface data for visualization
 !.... open file for contact forces between surfaces
 IMPLICIT NONE

   INTEGER(kind=4) :: flag
   !local variables
   INTEGER (kind=4) :: i,j,isurf,nsegm
   TYPE (surf3_db), POINTER :: surf
   TYPE (pair3_db), POINTER :: pair
   LOGICAL :: openp

   !.... Open file por post-process

   CALL openfi(41,flag=flag)    !open file for contact forces
   IF( flag /= 1)THEN
     WRITE(41,ERR=9999) npair,3  !number of forces and vector size
     pair => headp               !point to first pair
     DO i=1,npair                !for each pair
       WRITE(41,ERR=9999) i,pair%pname    !contact pair order and name
       pair => pair%next         !point to next pair
     END DO
   END IF

   !.... WRITE contact surface for post-proccessing
   IF( flag /= 1)THEN
     CALL openfi(18,flag = 1)        !open file
     surf => shead                   !pointer to first surface
     DO isurf=1,nsurf                !for each surface
       nsegm  = surf%nsegm           !number of segments defining the surface
       IF (surf%iwrit .AND. nsegm > 0 ) THEN  !if visualization desired
         WRITE(18,ERR=9999) nsegm,surf%sname  ! No of segments
         WRITE(18,ERR=9999) 3                 ! No of nodes per segment
         WRITE(18,ERR=9999) surf%lcseg(1:3,1:nsegm) !connectivities
       END IF
       surf => surf%next             !pointer to next surface
     END DO
     CLOSE (18)                      !close file
   END IF

   !.... Update global flags BOTTOM and TOP
   openp = .FALSE.          !for surface pressure
   pair => headp             !point to first pair
   DO i=1,npair              !for each pair
     IF (pair%press .OR. pair%wrink) THEN    !if surface pressure or wrinkles
       IF (.NOT.openp) THEN
         CALL openfi(44,flag=flag)    !file for surface data
         openp = .TRUE.
       END IF
       IF( flag /= 1)THEN
         surf => shead                 !pointer to first surface
         isurf = 1
         DO
           IF (isurf == pair%islav) EXIT
           isurf = isurf+1
           surf => surf%next
         END DO

         WRITE(44,ERR=9999) surf%sname,pair%press,pair%wrink,pair%cpress !name, press and wrinkles
         WRITE(44,ERR=9999) surf%ncnod,surf%nsegm            !No of nodes & segments in the surface
         WRITE(44,ERR=9999) (surf%lcnod(j),j=1,surf%ncnod)   !internal numeration of nodes
         WRITE(44,ERR=9999) (surf%lcseg(:,j),j=1,surf%nsegm) !surface connectivities
       END IF
     END IF
     pair => pair%next       !point to next pair
   END DO
   IF (openp .AND. flag /= 1) WRITE(44,ERR=9999) 'NONAME                        ',.FALSE.,.FALSE.,.FALSE. !mark to end

   IF( wear )THEN    !if friction work
     CALL openfi(45,flag=flag)
     IF( flag /= 1)THEN
       surf => shead                 !pointer to first surface
       DO j=1,nsurf
         IF (surf%nsegm > 0) THEN
           WRITE(45,ERR=9999) surf%nsegm,surf%iscod            !No of segments in the surface
           WRITE(45,ERR=9999) (surf%lcseg(:,i),i=1,surf%nsegm) !surface connectivities
         END IF
         surf => surf%next
       END DO
       WRITE(45,ERR=9999) 0,.FALSE.  !mark to end, associated to nsegm
     END IF
   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE prsur3
