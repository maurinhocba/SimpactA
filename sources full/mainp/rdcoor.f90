 SUBROUTINE rdcoor(actio)
 !-----------------------------------------------------------------------
 !     READ nodal coordinates & nodal systems
 !-----------------------------------------------------------------------
 USE param_db,ONLY: mnam,max_npoin_2d,max_npoin_3d
 USE ctrl_db,ONLY:  ndime, ndofn, nrotd, neulr, npoin, npoio, echo_chnode, nload, therm, ndoft
 USE outp_db, ONLY : updlon_outp, iwrit
 USE c_input
 USE ndinf_db
 USE nsets_db
 USE esets_db, ONLY : add_name
 USE npo_db !, ONLY: label,oldlb,euler,coord,coorc,coora,velnp,tempe,dtemp,resid
 USE cms_db, ONLY : nconm,updlon_cm
 USE damp_db, ONLY : updlon_damp
 USE loa_db, ONLY : updlon_lo
 USE kinc_db, ONLY : naris,nvelr
 USE sms_db, ONLY : selective_mass_scaling, updlon_sms
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*),INTENT(IN OUT):: actio !NEW, NSTRA0

   !--- Local variables
   CHARACTER(len=mnam) :: nsname
   INTEGER (kind=4):: k, n, nn, nrads, neu
   REAL (kind=8) :: fi, raddi, rm(9), theta, x(6), scale
   REAL (kind=8), PARAMETER :: facto = 1.7453292519943295769236907684886d-02   !pi/180

   TYPE (ndata) :: nd          ! nodal data
   TYPE (list_ndata) :: ndi    ! nodal information, temporary database
   TYPE (nset), POINTER :: ns  ! set of nodes

   INTERFACE
     INCLUDE 'angeul.h'
     INCLUDE 'delnod.h'
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
   END INTERFACE


   neu = 2*ndime - 3
   ! read initial control card  and print header
   CALL listen('RDCOOR')
   IF (.NOT.exists('GEOMET')) THEN
     IF (npoin == 0) THEN
       CALL runend('RDCOOR: GEOMETRY_DEFINITION EXPECTED')
     ELSE                          ! old problem and nothing to read
       backs = .TRUE.              !one line back
       npoio = npoin               ! keep number of points in previous mesh
       DEALLOCATE (oldlb)          ! release memory
       ALLOCATE (oldlb(npoio))     ! get memory
       DO n=1,npoin
         oldlb(n) = label(n)       ! store old labels for future reference
       END DO
       RETURN                      ! exit routine
     END IF
   END IF

   ! If this data set exists for actio == 'NSTRA0' => actio = 'NSTRA2'

   CALL ln_ini (ndi)  ! initialize ndi dbase

   IF (TRIM(actio) /= 'NEW') THEN

     ! Regenerates Euler angles if Nodal systems exists
     IF( neulr > 0 )THEN
       IF(ndime == 2)THEN
         euler = euler/facto        !transform angles to degrees
       ELSE IF( ndime == 3 ) THEN
         DO n=1,npoin                 !for each existing point
           rm = euler(1:9,n)          !Local nodal system
           CALL angeul(rm,euler(1:3,n),.TRUE.) !regenerates Euler angles (in degrees)
         END DO
       END IF
     END IF

     DO n = 1, npoin                !for each existing point
       nd%label = label(n)                  !store nodal label
       nd%coord(1:ndime) = coord(1:ndime,n) !store original coordinates
       nd%coorc(1:ndime) = coorc(1:ndime,n) !store present coordinates
       nd%coora(1:ndime) = coora(1:ndime,n) !store previous coordinates
       IF (neulr > 0) THEN
         nd%euler(1:neu,1) = eule0(1:neu,n)  !initial Euler ang.
         nd%euler(1:neu,2) = euler(1:neu,n)  !previous Euler ang.
       END IF
       nd%veloc(1:ndofn) = velnp(1:ndofn,n) !store nodal velocities
       IF (therm) nd%temp(1,1:ndoft) = tempe(:,n)  !store nodal temperature
       IF (therm) nd%temp(2,1:ndoft) = dtemp(:,n)  !store nodal old temperature
       CALL add_nd_at_end (ndi, nd)         !add to list (labels are SORTED)
     END DO

     npoio = npoin                          !keep previous number of nodes
     DEALLOCATE (oldlb)                     !forget
     ALLOCATE (oldlb(npoio))                !get memory
     DO k=1,npoio
       oldlb(k) = label(k)                  !keep previous label
     END DO
     DEALLOCATE ( coora, coorc, coord, euler, velnp, label) !release mem
     IF( ASSOCIATED (eule0 ))DEALLOCATE ( eule0) !release mem

     CALL delnod(iwrit, ndi)                !delete nodes if required

     npoin = get_npoin (ndi)                !present number of nodes
     WRITE (lures,"(/,'Number of nodes kept from the previous strategy: ',i6/ &
                  &  70('-'),/)",ERR=9999) npoin

   END IF

   CALL listen('RDCOOR')           !read a card
   IF (exists('GENERA')) THEN      !control card expected before node information
     WRITE (lures,"(/,80('-'),//,'  G E O M E T R Y   DEFINITION'/)",ERR=9999)
     scale=getrea('GSCALE',1d0,' SCALE Factor for Coordinates .....')
     nrads=getint('NRADS ',  0,' Coordinates Sys 0:CART 1:CYL 2:SPH')
   ELSE                            !previous card is compulsory, else STOP
     CALL runend('RDCOOR: GENERAL CARD EXPECTED      ')
   END IF

   IF(iwrit == 1) WRITE(lures,"(//'   ECHO of the Coordinate Data ',/ &
    & 5x,' Node      X',9X,'Y',9X,'Z',7X,'Alfa      Beta      Gamma')",ERR=9999)

   !                        initialize values for a new problem
   IF (TRIM(actio) == 'NEW')  CALL nsdb_ini      !initializes temporary D_B

   n = get_npoin (ndi)             !actual number of points in D_B

   DO                        !  GLOBAL LOOP
     CALL listen('RDCOOR')             !read a card
     IF(exists('ENDGEO')) EXIT         !check if end control card read
     IF(exists('SET   ',k))THEN          !check if set control card read
       nsname = get_name(posin=k,stype='NSET')        !set name
       CALL ns_ini (ns)                  !initializes list of nodes (set)
       CALL put_nsname (ns, nsname)      !assign name set
       CALL add_name (nsname,1)          !add to names list
       WRITE(lures,"('=============  nodes in set ',A,' ============' )",ERR=9999) &
             TRIM(nsname)
       DO
         CALL listen('RDCOOR')             !read a card
         IF( exists('ENDSET') )THEN        !end of set read
           CALL add_ns_at_end (ns)           !add set of nodes to the list of sets
           EXIT                              !cycle GLOBAL LOOP
         END IF
         n = n + 1                           !increase number of nodes and test MAX
         IF( (n > max_npoin_3d .AND. ndime == 3 ) .OR. (npoin > max_npoin_2d .AND. ndime == 2 )) &
         CALL runend('RDCOOR: Max. number of nodes EXCEDED !!!!!!')
         nn = INT(param(1))                  !node label
         IF( nn <= 0 )CALL runend('NODE label must be positive integers')
         IF( iwrit == 1) THEN                !If echo
           SELECT CASE (neulr)               !select how and what to echo
           CASE (0)
             WRITE(lures,"(i10,3f10.5)",ERR=9999) nn,param(2:ndime+1)
           CASE (1)
             WRITE(lures,"(i10,2f10.5,10x,f10.5)",ERR=9999) nn,param(2:4)
           CASE (9)
             WRITE(lures,"(i10,6f10.5)",ERR=9999) nn,param(2:7)
           END SELECT
         END IF

         CALL vecasi(6,param(2),x)           !for convenience only
         IF(neulr == 1) THEN                 !EULER angles in a 2-D problem
           rm(1) = x(3)*facto                !angle is the third parameter
         ELSE IF(neulr == 9) THEN            !EULER angles in a 3-D problem
           rm = 0d0                          !initializes Rotation matrix
           rm(1:3) = x(4:6)*facto            !first-second-third Euler angles
         END IF
         IF(nrads == 1) THEN                 !cylindrical coordinates (2-D or 3-D)
         ! ***    change cylindrical coordinates to cartesian
           raddi = x(1)                      !radius
           theta = x(2)*facto                !angle in radians
           x(1) = raddi*SIN(theta)
           x(2) = raddi*COS(theta)
         ELSE IF (nrads == 2) THEN           !Spherical coordinates (3-D problems)
         ! ***    change spherical coordinates to cartesian
           raddi = x(1)                      !radius
           theta = facto*x(2)                !first angle in radians
           fi    = facto*x(3)                !second angle in radians
           x(1) = raddi*SIN(fi)*SIN(theta)
           x(2) = raddi*SIN(fi)*COS(theta)
           x(3) = raddi*COS(fi)
         END IF

         x(1:ndime) = x(1:ndime)*scale       !scale coordinates

         nd%label = nn                       !assign read data to the list
         nd%coord(1:ndime) = x(1:ndime)      !orig. coord
         nd%coorc(1:ndime) = x(1:ndime)      !actual coord
         nd%coora(1:ndime) = x(1:ndime)      !previous coord
         IF (neulr > 0) nd%euler(1:neu,1) = rm(1:neu)  !euler angles
         IF (neulr > 0) nd%euler(1:neu,2) = rm(1:neu)  !euler angles
         nd%veloc(1:ndofn) = 0d0             !initializes velocity
         IF (therm) nd%temp = 0d0            !initializes temperature

         CALL l_insert (ndi, nd)             !add node to list of nodes
         CALL add_at_end (ns, nn)            !add node to SET
       END DO  !loop to read nodes in a set
     ELSE                            !error in data input
       CALL runend('RDCOOR: SET card for nodes expected')
     END IF

   END DO  !loop to read sets of nodes

   npoin = n                           !new number of nodes
              ! get memory for NPO_DB
   IF( (npoin <= max_npoin_3d .AND. ndime == 3 ) .OR. (npoin <= max_npoin_2d .AND. ndime == 2 )) THEN
     ALLOCATE (coord(ndime,npoin), coora(ndime,npoin), coorc(ndime,npoin),  &
               label(npoin), velnp(ndofn,npoin))
     IF( ndofn == 8 )THEN
       ALLOCATE( psia(2,npoin),psic(2,npoin))
       psia = 0d0
       psic = 0d0
     ELSE
       NULLIFY( psia ,psic )
     END IF
   ELSE
     STOP ' number of nodes exceeded for Academic version'
   END IF

   IF (therm) THEN ! allocation of nodal vectors for thermal analysis
     ALLOCATE ( tempe(ndoft,npoin), dtemp(ndoft,npoin) )
     tempe = 0d0
     dtemp = 0d0
   END IF

   IF(neulr > 0)THEN                   !IF nodal systems
     ALLOCATE(eule0(neu,npoin), euler(neulr,npoin))      !get memory
     euler = 0d0                       !initializes
   ELSE
     ALLOCATE( euler (1,1) )           !to avoid a null pointer
   END IF

   CALL l_head (ndi)                   !go to top of the list of nodes
   DO n=1,npoin                        !for each node in the list
     nd = get_ndata (ndi)              !get nodal data and assign
     label(n) = nd%label
     coord(1:ndime,n) = nd%coord(1:ndime)
     coorc(1:ndime,n) = nd%coorc(1:ndime)
     coora(1:ndime,n) = nd%coora(1:ndime)
     velnp(1:ndofn,n) = nd%veloc(1:ndofn)
     IF (neulr > 0) THEN
       eule0(1:neu,n) = nd%euler(1:neu,1)
       euler(1:neu,n) = nd%euler(1:neu,2)
     END IF
     IF (therm) tempe(:,n) = nd%temp(1,1:ndoft)
     IF (therm) dtemp(:,n) = nd%temp(2,1:ndoft)
     CALL l_next (ndi,k)                 !next node in the list
   END DO

   CALL dalloc_list_ndata (ndi)  ! deallocate all the temporary database

   IF (TRIM(actio) == 'NEW') THEN   !for a NEW problem
     npoio = npoin               !old values are present values
     ALLOCATE (oldlb(npoio))
     DO n=1,npoin
       oldlb(n) = label(n)
     END DO
   ELSE                          !for subsequent strategies
     actio = 'NSTRA2'            !modify action for rest of the TASKS
     ! updating local node numbers in existing problem
     echo_chnode = .FALSE.
     CALL elemnt ('UPDLON')               ! elements/segments/nodes
     CALL updlon_kc (naris,nvelr,ndime)   ! KINEMATIC conditions
     IF(nconm)CALL updlon_cm (nrotd)      ! CONCENTRATED MASSES
     IF(selective_mass_scaling)CALL updlon_sms ( ) ! Selective Mass Scaling
     CALL updlon_outp (oldlb)             ! output nodes
     CALL updlon_damp ()                  ! DAMPING nodes
     CALL updlon_lo (nload,oldlb)         ! LOADS
     CALL contac ('UPDLON',0)             ! CONTACT nodes
     echo_chnode = .TRUE.
     DEALLOCATE (oldlb)
     ALLOCATE (oldlb(npoin))
     DO n=1,npoin
       oldlb(n) = label(n)
     END DO
   END IF

 RETURN
  9999 CALL runen2('')
 END SUBROUTINE rdcoor
