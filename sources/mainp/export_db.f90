 MODULE export_db

  ! module to export element set data

  USE param_db,ONLY: mnam
  IMPLICIT NONE
  SAVE

  INTEGER (kind=4) :: neset !Number of Export SETs

  TYPE export_set
    CHARACTER(len=mnam) :: name, &  !label of the element set to write
                           fname    !filename to write data
    LOGICAL :: ascii,            &  !
               cooro,            &  !
               matsc,            &  !
               bound,            &  !
               displ,            &  !
               inter                !
  END TYPE export_set

  TYPE (export_set), ALLOCATABLE :: exp_d(:)

 CONTAINS

 SUBROUTINE exp_data ( )
 !
 !   read data sets for future use
 !
 USE param_db,ONLY:  mnam
 USE c_input
 USE esets_db, ONLY : search_name

 IMPLICIT NONE
 CHARACTER(len=mnam) :: name     !label of the element set to write
 LOGICAL :: found         !flag
 INTEGER (kind=4) :: is,ln,stype

 CALL listen('EXPDAT')            !read a line
 IF( .NOT.exists('EXPORT') )THEN  !see if Export data is present
   backs = .TRUE.                 !back one line
   neset = 0                      !no sets to export
   IF(ALLOCATED(exp_d))DEALLOCATE(exp_d)  !release space
   RETURN                         !back to main stream
 END IF

 neset = INT(param(1))            !number of sets to export
 IF( neset == 0 )neset = 20       !no data, set to maximum
 IF(ALLOCATED(exp_d))DEALLOCATE(exp_d) !if previous data exists, nullify it
 ALLOCATE(exp_d(neset))           !get memory for data

 is = 0                           !initializes
 stype = 2                        !element set (for search_name)
 DO  ! loop for each Element Set to Export
   CALL listen('EXPORT')         !read a line
   IF( exists('ENDEXP') )EXIT    !exit loop

   name = get_name('ELMSET',found,'!Element set to Export:')     !set name
   ln = LEN_TRIM(name)                          !lenght of string

   CALL search_name (name, found, stype=stype)
   IF( .NOT.found ) THEN
     WRITE(lures,"('  element set to export not found ',a)",ERR=9999) name(1:ln)
     CALL runend('EXPORT: Elm Set not found ')
   END IF

   is = is + 1                   !increment counter
   exp_d(is)%name = name         !element set name
   CALL listen('EXPORT')         !read a line (data line)
   exp_d(is)%fname = get_name('FILE  ',found, '!Export File Name:')  !file name
   exp_d(is)%ascii = exists('TYPE  ') .AND. exists('ASCII')          !ASCII or binary
   exp_d(is)%cooro = .NOT.(exists('GEOMET').AND.exists('STAGE ')) ! Original or actual coordinates will be written
   exp_d(is)%matsc = exists('MATERI')    ! material/section properties will be written
   exp_d(is)%bound = exists('BOUNDA')    ! boundary conditions will be written
   exp_d(is)%displ = exists('DISPLA')    ! displacements will be written
   exp_d(is)%inter = exists('INTERN')    ! internal variables will be written
 END DO
 neset = is                              !keep number of export sets
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE exp_data

 SUBROUTINE dump_exp_data ( )
 ! dumps data in restar file
 IMPLICIT NONE
 INTEGER (kind=4) :: i

 WRITE(50,ERR=9999) neset
 IF( neset > 0 )THEN
   DO i=1,neset
     WRITE(50,ERR=9999) exp_d(i)%name,exp_d(i)%fname,exp_d(i)%ascii,exp_d(i)%cooro, &
                        exp_d(i)%matsc,exp_d(i)%bound,exp_d(i)%displ,exp_d(i)%inter
   END DO
 END IF
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump_exp_data

 SUBROUTINE rest_exp_data ( )
 ! reads data in restar file
 IMPLICIT NONE
 INTEGER (kind=4) :: i

 READ(51) neset
 IF( neset > 0 )THEN
   ALLOCATE(exp_d(neset))
   DO i=1,neset
     READ(51) exp_d(i)%name,exp_d(i)%fname,exp_d(i)%ascii,exp_d(i)%cooro, &
              exp_d(i)%matsc,exp_d(i)%bound,exp_d(i)%displ,exp_d(i)%inter
   END DO
 END IF

 RETURN
 END SUBROUTINE rest_exp_data

 SUBROUTINE export ( )
 !
 !   Export element data set for future use
 !
 USE param_db,ONLY: mnam
 USE c_input
 USE meshmo_db, ONLY : nodset,numpo
 USE ctrl_db, ONLY: ndime, neulr, ndofn
 USE esets_db, ONLY: nelms,rot_free
 USE npo_db, ONLY : label,coord,coora,euler,ifpre,iffix,eule0,naeul
 USE mat_dba

 IMPLICIT NONE
 REAL (kind=8), PARAMETER :: facto = 1.74532925199433d-02
 CHARACTER(len=mnam) :: name, &  !label of the element set to write
                        fname    !filename to write data
 LOGICAL :: found,ascii,cooro,matsc,bound,displ,inter,rest  !flags
 INTEGER (kind=4) :: is,i,j,k,l,n,chnode,ln,neu,nm,ns,bb(7),nd
 INTEGER (kind=4) :: auxi(100)
 LOGICAL, ALLOCATABLE :: secs(:),mats(:)
 REAL(kind=8) :: angle(3)
 TYPE(mater), POINTER :: mat
 TYPE(section), POINTER :: sect
 TYPE(curve), POINTER :: cur
 TYPE(postv), POINTER :: postp

 INTERFACE
   INCLUDE 'angeul.h'
   INCLUDE 'elemnt.h'
 END INTERFACE


 DO is=1,neset ! loop for each Element Set to Export

   name = exp_d(is)%name
   ln = LEN_TRIM(name)                          !lenght of string
   ! extract a set of nodes  set==>nodset(numpo)
   CALL elemnt ('NODSET', name=name, flag2=found)
   ! recover internal numeration
   DO i=1,numpo       !number of nodes in the set
     n = chnode(nodset(i)) !internal number
     nodset(i) = n         !keep internal number for easier use
   END DO

   fname = exp_d(is)%fname ! file name
   ascii = exp_d(is)%ascii ! ascii of binary
   cooro = exp_d(is)%cooro ! Original or actual coordinates will be written
   matsc = exp_d(is)%matsc ! material/section properties will be written
   bound = exp_d(is)%bound ! boundary conditions will be written
   displ = exp_d(is)%displ ! displacements will be written
   inter = exp_d(is)%inter ! internal variables will be written

   neu = neulr
   IF( neu == 9 )neu = 3
   ! Open file and write coordinates
   IF( ascii )THEN
     CALL openfi(10,fname)
     WRITE(10,"('GEOMETRY_DEFINITION',/,'GENERAL:  GSCALE= 1.00',/ &
              & 'SET ',a)",ERR=9999) name(1:ln)
     DO i=1,numpo               !for each node in the set
       n = nodset(i)            !internal number
       IF(neulr > 0)THEN          !If local system exists
         IF( naeul(i) )THEN
           IF(neulr == 1)THEN
             angle(1) = euler(1,n)/facto     !transform angle to degrees
           ELSE IF( neulr == 9 ) THEN
             CALL angeul(euler(1:9,n),angle) !regenerates Euler angles (in degrees)
           END IF
         ELSE
           angle(1:neu) = 0d0
         END IF
         IF(cooro)THEN  !original coordinates
           WRITE(10,"(i8,9e13.5)",ERR=9999) label(n),coord(:,n), &
                 eule0(1:neu,n)/facto,angle(1:neu)
         ELSE           !coordinates at the end of the stage
           WRITE(10,"(i8,9e13.5)",ERR=9999) label(n),coora(:,n), &
                 eule0(1:neu,n)/facto,angle(1:neu)
         END IF
       ELSE
         IF(cooro)THEN  !original coordinates
           WRITE(10,"(i10,3e15.6)",ERR=9999) label(n),coord(:,n)
         ELSE           !coordinates at the end of the stage
           WRITE(10,"(i10,3e15.6)",ERR=9999) label(n),coora(:,n)
         END IF
       END IF
     END DO
     WRITE(10,"('END_SET ',a,/,'END_GEOMETRY_DEFINITION')",ERR=9999) name(1:ln)
   ELSE  !binary output
     CALL openfi(46,fname)
     WRITE(46,ERR=9999) name
     WRITE(46,ERR=9999) cooro,matsc,bound,displ,inter
     WRITE(46,ERR=9999) numpo,ndime,neulr
     DO i=1,numpo       !number of nodes in the set
       n = nodset(i) !internal number
       IF(neulr > 0)THEN    !Euler angles
         IF( naeul(n) )THEN
           IF(neulr == 1)THEN
             angle(1) = euler(1,n)                  !angles in radians
           ELSE IF( neulr == 9 ) THEN
             CALL angeul(euler(1:9,n),angle,.TRUE.) !regenerates Euler angles (in radians)
           END IF
         ELSE
           angle(1:neu) = 0d0
         END IF
         IF(cooro)THEN  !original coordinates
           WRITE(46,ERR=9999) label(n),coord(:,n),eule0(1:neu,n),angle(1:neu)
         ELSE           !coordinates at the end of the stage
           WRITE(46,ERR=9999) label(n),coora(:,n),eule0(1:neu,n),angle(1:neu)
         END IF
         IF( displ)  WRITE(46,ERR=9999) coora(:,n)-coord(:,n)  !Write displacements
       ELSE
         IF(cooro)THEN  !original coordinates
           WRITE(46,ERR=9999) label(n),coord(:,n)
         ELSE           !coordinates at the end of the stage
           WRITE(46,ERR=9999) label(n),coora(:,n)
         END IF
         IF( displ)  WRITE(46,ERR=9999) coora(:,n)-coord(:,n)  !Write displacements
       END IF
     END DO
   END IF

   ! materials and sections
   IF( matsc ) THEN
     ! first search in the element data set the sections used => AUXI
     auxi = -1   !initializes to none
     CALL elemnt ('SECDAT', name=name, ivect=auxi, flag2=found )  !search
     ALLOCATE( secs(nusect), mats(numats) )          !auxiliar arrays
     secs = .FALSE.        !initializes
     mats = .FALSE.        !initializes
     i = 1
     DO                        !for each value in AUXI
       k = auxi(i)             !section label
       IF( k < 0 )EXIT         !if not a valid section label (exit)
       secs(k) = .TRUE.
       sect => psecs(k)%p      !point to section
       k = sect%mtbas%matno         !associated material label
       mat => mhead                !point to first material
       l = 1                       !initializes
       DO
         IF(mat%matno == k )THEN     !if material label found
           mats(l) = .TRUE.          !set to write
           EXIT                      !exit loop
         END IF
         l = l+1                     !increase counter
         mat => mat%next             !point to next material
       END DO
       EXIT                        !exit loop
       i = i + 1
     END DO

     !  export materials and sections to disk
     IF( ascii )THEN
       WRITE(10,"('MATERIAL_DEFINITION',/)",ERR=9999)
       mat => mhead                    !point to first material
       DO i=1,numats                   !for each possible material
         IF( .NOT.mats(i)) CYCLE       !material not used in the set
         ! ************
         ! this is an ackward task due to the large number of possibilities
         ! ************
         mat => mat%next               !point to next material
       END DO
       WRITE(10,"('END_MATERIAL_DEFINITION',/)",ERR=9999)

       WRITE(10,"('SECTION_DEFINITION',/)",ERR=9999)
       sect => shead                   !point to first section
       DO i=1,nusect                   !for each possible section
         IF( .NOT.secs(i)) CYCLE       !section not used in the set
         ! ************
         ! this is an easier task (possibilities are more restricted)
         ! ************
         sect => sect%next
       END DO
       WRITE(10,"('END_SECTION_DEFINITION',/)",ERR=9999)

     ELSE   !for binary file just download the materials and sections
            ! as it is done in DUMP_MATE
       nm = 0                     !count the materials for reference
       DO i=1,numats
         IF( mats(i) )nm = nm+1
       END DO
       ns = 0                     !count the sections for reference
       DO i=1,nusect
         IF( secs(i) )ns = ns+1
       END DO
       WRITE(46,ERR=9999) nm,ns   !number of materials and sections

       mat => mhead     !point to first material
       DO i=1,numats    !for each material
         IF( .NOT.mats(i)) CYCLE  !not used in the set
         WRITE (46,ERR=9999) mat%matno, mat%mtype, mat%matdef !write all the control variables
         WRITE (46,ERR=9999) mat%prope, mat%propp, mat%props  !write all the properties
         ! write associated curves if exists
         IF( mat%matdef(12) > 0 )THEN
           cur => mat%chead        !point to first curve
           DO j=1,mat%matdef(12)   !for each curve
             WRITE(46,ERR=9999) cur%np       !number of points in the curve
             WRITE(46,ERR=9999) cur%val      !curve values
             cur => cur%next       !point to next curve
           END DO
         END IF
         mat => mat%next
       END DO

       sect => shead    !point to first section
       DO i=1,nusect   !for each section
         IF( .NOT.secs(i)) CYCLE  !not used in the set
         WRITE (46,ERR=9999) sect%secno, sect%secty, sect%mabas, sect%secdef !write all the control variables
         WRITE (46,ERR=9999) sect%iprop, sect%rprop !write all the properties
         IF( sect%secdef(4) > 0)THEN
           postp => sect%postp
           DO j=1,sect%secdef(4)
             WRITE (46,ERR=9999) postp%type,postp%dim,postp%name
             postp => postp%next
           END DO
         END IF
         IF( sect%secdef(5) >= 3)THEN  !if the section may have FLC
           IF( sect%iprop(3) > 0)THEN  !If the FLC label is greater than 0
             CALL wrflc(sect%iprop(3),46)  !write FLC curve on file 46
           END IF
         END IF
         sect => sect%next
       END DO
     END IF
   END IF

   ! element connectivities and other Element set data
   IF( ascii ) THEN  !  istop = 0  => print connectivities only
     WRITE(10,"('SET_DEFINITION')",ERR=9999)
     CALL elemnt ('EXPORT', name=name, flag1=ascii, flag2=found, istop=0)
     WRITE(10,"('  END_ELEMENT_DEFINITON',/,'END_SET_DEFINITON')",ERR=9999)
   ELSE              ! binary files (all elemental information)
     WRITE(46,ERR=9999) numpo     !number of nodes in the set
     WRITE(46,ERR=9999) (label(nodset(i)),i=1,numpo)
     k = 0               !do NOT include internal variables
     IF( inter ) k = 1   !include internal variables
     CALL elemnt ('EXPORT', name=name, flag1=ascii, flag2=found, istop=k)
   END IF

   ! boundary conditions
   !  The following approach is based on IFPRE and IFFIX only
   IF( bound )THEN
     nd = ndofn
     IF( rot_free )nd = nd+1
     IF( ascii )THEN
       WRITE(10,"('KINEMATIC_CONDITIONS',/,'BOUNDARY_CONDITIONS')",ERR=9999)
     ELSE
       ns = 0          !initializes counter of restrained nodes
       CALL openfi(47) !scratch binary file
     END IF
     !                loop over the nodes in the set
     DO i=1,numpo     !for each node
       n = nodset(i)  !internal node
       rest =  .FALSE.             !initializes
       bb = 0                      !initializes
       DO j=1,ndofn                          !for each dof
         IF( ifpre(j,n) <= 0 )THEN  !nn should be used instead
           bb(j) = 1      !set to constrained
           rest = .TRUE.
         END IF
       END DO
       IF( rot_free )THEN
         bb(nd) = iffix(n)         !
         IF( iffix(n) > 0 )rest = .TRUE.
       END IF
       IF( rest )THEN                          !if constrained
         IF( ascii )THEN                  !for ASCII file
           WRITE(10,"(i10,4x,7i2)",ERR=9999) label(n),bb(1:nd)
         ELSE                             !for binary file print to auxiliar file
           ns = ns+1                      !increase counter
           WRITE(47,ERR=9999) label(n),bb(1:nd)
         END IF
       END IF
     END DO
     IF( ascii )THEN
       WRITE(10,"('END_BOUNDARY_CONDITIONS',/,'END_KINEMATIC_CONDITIONS')",ERR=9999)
     ELSE                       !for binary file
       REWIND(47)               !go to first
       WRITE(46,ERR=9999) numpo,ndime,neulr     !control variables
       WRITE(46,ERR=9999) (label(nodset(i)),i=1,numpo)
       WRITE(46,ERR=9999) ns,nd
       DO i=1,ns
         READ (47) n,bb(1:nd)
         WRITE(46,ERR=9999) n,bb(1:nd)
       END DO
       CLOSE(47)
     END IF
   END IF

   ! initial conditions
   IF(( displ .OR. inter ) .AND. ascii)THEN
     WRITE(10,"('INITIAL_CONDITIONS')",ERR=9999)
     ! initial displacements
     IF( displ )THEN
       WRITE(10,"('DISPLACEMENTS')",ERR=9999)
       DO i=1,numpo       !number of nodes in the set
         n = nodset(i) !internal number
         WRITE(10,"(i10,3e15.6)",ERR=9999) label(n),coora(:,n)-coord(:,n)
       END DO
       WRITE(10,"('END_DISPLACEMENTS')",ERR=9999)
     END IF

     ! internal variables
     IF( inter )THEN
       WRITE(10,"('GAUSS_INITIALIZATION ',a)",ERR=9999) name(1:ln)
       !          istop = 1  => print internal variables
       CALL elemnt ('EXPORT', name=name, flag1=ascii, flag2=found, istop=1)
       WRITE(10,"('END_GAUSS_INITIALIZATION')",ERR=9999)
     END IF
     WRITE(10,"('END_INITIAL_CONDITIONS')",ERR=9999)
   END IF
   IF( ascii )THEN
     CLOSE(10)
   ELSE
     CLOSE(46)
   END IF
 END DO
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE export
 END MODULE export_db
