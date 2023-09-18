 SUBROUTINE trarot (ecdis,ircon)
!
!  this routines controls Rigid body translation and rotation
!  of elements sets as strategy
!
 USE param_db,ONLY: mnam,mich
 USE ctrl_db, ONLY: ndime, npoin,neulr,nstra,text
 USE meshmo_db, ONLY : nodset,numpo,strnam
 USE nsets_db
 USE c_input
 USE npo_db, ONLY : label,coora,coorc,euler
 USE gvar_db, ONLY : actchk

 IMPLICIT NONE
 ! dummy arguments
   INTEGER  (kind=4), INTENT(IN) ::     ecdis     ! control DOF for output
   INTEGER  (kind=4), INTENT(IN OUT) :: ircon     ! counter of restart files
 ! Local variables
 LOGICAL :: found
 INTEGER (kind=4) :: chnode,i,k,rnod
 REAL (kind=8) :: tvec(ndime),rvec(ndime),rang,rpoi(ndime),q(ndime),r(ndime), &
                  lb(ndime,ndime),eul(9)
 TYPE (nset), POINTER :: ns
 CHARACTER (len=mnam) :: sname
 !Functions
 CHARACTER(len=mich):: inttoch   !Puts an integer into a string

 INTERFACE
   INCLUDE 'elemnt.h'
   INCLUDE 'outdyn.h'
 END INTERFACE

 WRITE(*,"(/,5X,'***  BEGIN STAGE ',A,': ',A,/)") TRIM(inttoch(nstra,0)),text
 WRITE(55,"(/,5X,'***  BEGIN STAGE ',A,': ',A,/)",ERR=9999) TRIM(inttoch(nstra,0)),text
 CALL timing(11,1)              !mesh modifications

 WRITE (*, 1002        )        !send a message to screen
 WRITE (lures,1002,ERR=9999)     !write a message to output file
 1002 FORMAT (4x, '--- T R A N S L A T I N G   &   R O T A T I N G--',//)

 IF( .NOT.actchk )THEN   ! write initial data for post-processing
   CALL wrtpos(0)
 END IF


 DO                                     !loop over element_sets to modify
   CALL listen ('TRAROT')               !read a line
   IF(exists('ENDSTR'))EXIT             !exit loop
   !          read translation vector
   IF ( exists('TRAVEC',i)) THEN        !if key-word exists
     DO k=1,ndime
       tvec(k) = param(i)
       i = i+1
     END DO
   ELSE                                 !no translation vector
     tvec = 0d0
   END IF
   !          read rotation vector (axis)
   rpoi = 0d0
   IF( ndime == 3 )THEN                 !for 3-D problems
     IF ( exists('ROTVEC',i)) THEN      !if key word exists
       DO k=1,ndime
         rvec(k) = param(i)
         i = i+1
       END DO
       CALL vecuni(ndime,rvec,rang)     !unit vector and rotation angle
       IF( rang /= 0d0 ) THEN           !if rotation vector is not null
         IF( exists('ROTANG',i)) rang = param(i)  !read rotation angle
         IF( exists('ROTPOI',i))THEN    !read coordinates of fixed point
           DO k=1,ndime
             rpoi(k) = param(i)
             i = i+1
           END DO
           rnod = 0
         ELSE IF( exists('ROTNOD',i) )THEN  !or read fixed node label
           rnod = INT(param(i))
           rnod = chnode(rnod)
           rpoi = coora(:,rnod)
         ELSE
           CALL runend('TRAROT: a fixed point must be defined')
         END IF
       END IF
     ELSE
       rvec = 0d0
     END IF
   ELSE         !for 2-D problems (rvec = Z)
     rang = 0d0
     IF( exists('ROTANG',i)) rang = param(i)
     IF( exists('ROTPOI',i))THEN    !read coordinates of the fixed point
       DO k=1,ndime
         rpoi(k) = param(i)
         i = i+1
       END DO
       rnod = 0
     ELSE IF( exists('ROTNOD',i))THEN !or read the fixed node label
       rnod = INT(param(i))
       rnod = chnode(rnod)
       rpoi = coora(:,rnod)
     ELSE IF( rang /= 0 )THEN                         !check
       CALL runend('TRAROT: a fixed point must be defined')
     END IF

   END IF
   rang = rang*ATAN(1d0)/45d0              !angle from degrees to radians
   IF( ndime == 3 )rvec = rvec*rang        !rotation vector

   ! compute rotation matrix
   IF( ndime == 3 )THEN
     CALL expot8(rvec(1),lb(1,1))
   ELSE
     lb = RESHAPE((/ COS(rang), -SIN(rang), SIN(rang), COS(rang) /), &
                  (/ ndime,ndime /))
   END IF
                          !modify associated sets
   coorc = coora          !initializes
   DO
     CALL listen('TRAROT')              !read a line
     IF(exists('ENDSET'))EXIT           !exit
     IF ( exists('ELMSET',k)) THEN      !an element set
       sname = get_name(posin=k,stype='ESET')        !Element set name
       WRITE(lures,"(/'element set name to translate and rotate: ',A)",ERR=9999) TRIM(sname)
       ! node set definition based on an element set
       CALL elemnt ('NODSET', name=sname, flag2=found) !get nodes in ELM_SET
       IF (.NOT.found) CALL runend('TRAROT:ELEMENT SET NOT FOUND       ')

     ELSE IF (exists('SET   ',k)) THEN  !a node set
       sname = get_name(posin=k,stype='NSET')         !node set name
       WRITE(lures,"(/'nodes set name to translate and rotate: ',A)",ERR=9999) TRIM(sname)
       ! node set definition based on an node set
       CALL nsdb_search (sname, found, ns)          !search in list of sets
       IF (found) THEN                      !set exists, position is NS
         numpo = get_length (ns)            !number of nodes in the set
         ALLOCATE (nodset(numpo))           !get memory
         CALL ns_head (ns)                  !go to top of list
         DO i =1, numpo                     !for each node in the list
           nodset(i) = get_label(ns)        !get node label
           CALL ns_next (ns)                !go to next node
         END DO
       ELSE                                 !else SET name not found => error
         CALL runend('TRAROT: NODE  SET NOT FOUND       ')
       END IF
     ELSE
       CALL runend ('trarot: ELM_SET or SET name expected.')
     END IF

     IF( actchk ) CYCLE !do not modify coordinates
     ! modify coordinates
     DO i=1,numpo                        !for each node in the set
       k = chnode(nodset(i))                     !global node
       q = coora(:,k) - rpoi             !vector from fixed point
       r = MATMUL(lb,q)                  !rotated vector
       coorc(:,k) = rpoi + tvec + r      !new coordinates
     END DO

     IF( neulr > 0 )THEN                 !if Euler angles exists
       IF( ndime == 3 )  CALL expot8(rvec(1),lb(1,1))    !rotation matrix

       DO i=1,numpo                     !for each node in the set
         k = chnode(nodset(i))                  !global node
         IF(ndime == 3) THEN            !for 3-D problems
           eul = euler(:,k)
           CALL proma1(euler(1,k),lb(1,1),eul(1),3,3,3)
         ELSE                           !for 2-D problems
           euler(1,k) = euler(1,k) + rang
         END IF
       END DO
     END IF

     DEALLOCATE (nodset)         !release memory
   END DO
   coora = coorc                       !same actual and previous coordinates
 END DO
 IF( .NOT.actchk )THEN
   ! write data and results for post-processing
   CALL outdyn(.FALSE.,.TRUE., ecdis)
   CALL dumpin(ircon,.TRUE.,.FALSE.,.FALSE.)
 END IF
 CALL timing(11,2)              !mesh modifications

 WRITE(*,"(5X,'***  END OF STAGE ',A,': ',A,/)") TRIM(inttoch(nstra,0)),text
 WRITE(55,"(5X,'***  END OF STAGE ',A,': ',A,/)",ERR=9999) TRIM(inttoch(nstra,0)),text

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE trarot
