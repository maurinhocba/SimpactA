 SUBROUTINE imcoor (nimpf,lab0,impda,els_name,actio)
 !-----------------------------------------------------------------------
 !     READ nodal coordinates & nodal systems from binary files
 !-----------------------------------------------------------------------
 USE param_db,ONLY: mnam,max_npoin_2d,max_npoin_3d
 USE ctrl_db,ONLY: ndime, ndofn, neulr, npoin, npoio, echo_chnode, nload, therm, ndoft
 USE outp_db, ONLY : updlon_outp, iwrit
 USE c_input
 USE ndinf_db
 USE nsets_db
 USE esets_db, ONLY : add_name
 USE npo_db
 USE cms_db, ONLY : nconm,updlon_cm
 USE damp_db, ONLY : updlon_damp
 USE loa_db, ONLY : updlon_lo
 USE kinc_db, ONLY : naris,nvelr
 IMPLICIT NONE
   !--- Dummy arguments
   INTEGER (kind=4), INTENT(IN) :: nimpf,lab0(nimpf)
   LOGICAL, INTENT(IN) :: impda(9,nimpf)
   CHARACTER(len=mnam), INTENT(IN) :: els_name(nimpf)
   CHARACTER(len=*),INTENT(IN OUT):: actio


   !--- Local variables
   CHARACTER(len=mnam) :: nsname
   INTEGER (kind=4):: i,n,nn,idime,ieulr,np,nv,nf,is,ir,neu
   REAL (kind=8) :: rm(9), x(9),d(3)
   LOGICAL :: cooro,displ,overw,found,chlab
   REAL (kind=8), PARAMETER :: facto = 1.7453292519943295769236907684886d-02   !pi/180

   TYPE (ndata) :: nd          ! nodal data
   TYPE (list_ndata) :: ndi    ! nodal information, temporary database
   TYPE (nset), POINTER :: ns  ! set of nodes

   INTERFACE
     INCLUDE 'angeul.h'
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
   END INTERFACE

   neu = 2*ndime - 3
   chlab = .FALSE.
   CALL ln_ini (ndi)  ! initialize ndi dbase

   IF (TRIM(actio) /= 'NEW') THEN
     ! Regenerates Euler angles if Nodal systems exists
     IF( neulr > 0 )THEN
       IF( ndime == 3 ) THEN
         DO n=1,npoin                 !for each existing point
           rm = euler(1:9,n)          !Local nodal system
           CALL angeul(rm,euler(1:3,n),.TRUE.) !regenerates Euler angles (in Rads)
           euler(4:6,n) = 0d0         !no values
         END DO
       END IF
     END IF

     DO n = 1, npoin                !for each existing point
       nd%label = label(n)                  !store nodal label
       nd%coord(1:ndime) = coord(1:ndime,n) !store original coordinates
       nd%coorc(1:ndime) = coorc(1:ndime,n) !store present coordinates
       nd%coora(1:ndime) = coora(1:ndime,n) !store previous coordinates
       IF (neulr > 0) THEN
         nd%euler(1:neu,1) = eule0(1:neu,n)  !initial Euler ang. (Rads)
         nd%euler(1:neu,2) = euler(1:neu,n)  !previous Euler ang.  (Rads)
       END IF
       nd%veloc(1:ndofn) = velnp(1:ndofn,n) !store nodal velocities
       IF (therm) nd%temp(1,1:ndoft) = tempe(:,n)    !store nodal temperature
       IF (therm) nd%temp(2,1:ndoft) = dtemp(:,n)    !store nodal incremental temperature
       CALL add_nd_at_end (ndi, nd)         !add to list (labels are SORTED)
     END DO

     npoio = npoin                          !keep previous number of nodes
     DEALLOCATE (oldlb)                     !forget
     ALLOCATE (oldlb(npoio))                !get memory
     oldlb = label                          !keep previous label
     DEALLOCATE ( coora, coorc, coord, euler, velnp, label) !release mem
     IF( ASSOCIATED (eule0 ))DEALLOCATE ( eule0) !release mem

     npoin = get_npoin (ndi)                !present number of nodes

     WRITE (lures,"(/,'Number of nodes kept from the previous strategy: ',i6/ &
                  &  70('-'),/)",ERR=9999) npoin
   END IF
   n = npoin                              !actual number of points in D_B

   IF(iwrit == 1) WRITE(lures,"(//'   ECHO of the Coordinate Data ',/ &
    & 5x,' Node      X',9X,'Y',9X,'Z',7X,'Alfa      Beta      Gamma')",ERR=9999)

   DO is=1,nimpf            !  GLOBAL LOOP over each imported element set
     nf = 70+is
     READ(nf)np,idime,ieulr
     IF( idime /= ndime )THEN
       WRITE(lures,"(' Ndime of present problem is:',i5)",ERR=9999) ndime
       WRITE(lures,"(' While imported set ',a,' has ndime:',i5)",ERR=9999) nsname,idime
       CALL runend('IMPORT: NDIME must be the same')
     END IF
     IF( ieulr > neulr )THEN  !think about it
       neulr = 1
     END IF
     nv = ndime
     IF( ieulr > 0 ) nv = nv + 2*neu
     nsname = els_name(is)  !element set name
     cooro = impda(4,is)    !original or actual
     displ = impda(7,is)    !displacement included
     overw = impda(9,is)    !overwrite previous coordinates

     IF( overw )THEN        !if nodes exist and must be overwrited
       DO i=1,np                          !loop over all the nodes
         READ(nf)nn,x(1:nv)               !node label and coordinates
         IF( displ ) READ(nf) d(1:ndime)  !displacements
         !IF(iwrit==1)WRITE(lures,"(i10,3f10.5)",ERR=9999) nn,x(1:nv)
         IF(ieulr == 1) THEN                 !EULER angles in a 2-D problem
           rm(1:2) = x(3:4)                  !angle(rad) is the third parameter
         ELSE IF(ieulr == 9) THEN            !EULER angles in a 3-D problem
           rm(1:6) = x(4:9)                  !Euler angles (rads)
         END IF
         IF( neulr > ieulr ) rm = 0d0        !
         CALL l_search (ndi,nn,found)        ! search node in data-base
         IF(.NOT.found )THEN  !troubles
           WRITE(*,9000)TRIM(nsname),nn
           WRITE(lures,9000)TRIM(nsname),nn
           WRITE(55,9000)TRIM(nsname),nn
         END IF
         IF( cooro )THEN
           ndi%posic%data%coord(1:ndime) = x(1:ndime)      !orig. coord
           IF( displ )THEN
              ndi%posic%data%coora(1:ndime) = x(1:ndime)+d(1:ndime)      !actual coord
           ELSE
              ndi%posic%data%coora(1:ndime) = x(1:ndime)      !actual coord
           END IF
         ELSE
           ndi%posic%data%coora(1:ndime) = x(1:ndime)      !stage coord
           IF( displ )THEN
              ndi%posic%data%coord(1:ndime) = x(1:ndime)-d(1:ndime)      !actual coord
           ELSE
              ndi%posic%data%coord(1:ndime) = x(1:ndime)   !stage coord
           END IF
         END IF
         ndi%posic%data%coorc(1:ndime) = x(1:ndime)      !actual coord
         ndi%posic%data%veloc(1:ndofn) = 0d0             !initializes velocities

         IF (neulr > 0) THEN
           ndi%posic%data%euler(1:neu,1) = rm(1:neu)        !initial angles (rad)
           ndi%posic%data%euler(1:neu,2) = rm(neu+1:2*neu)  !present angles (rad)
         END IF
       END DO

     ELSE  !if new nodes are included

       chlab = .TRUE.
       IF( impda(2,is) )THEN   !change node labels
         IF( impda(3,is) )THEN !SEQUEntial
           ir = 2
         ELSE
           ir = 3
         END IF
       ELSE                      !keep node labels
         ir = 1
       END IF
       CALL ns_ini (ns)                  !initializes list of nodes (set)
       CALL put_nsname (ns, nsname)      !assign name set
       CALL add_name (nsname,1)          !add to names list
       WRITE(lures,"('=============  nodes in set ',A,' ============' )",ERR=9999) &
             TRIM(nsname)
       i = n + np                          !increase number of nodes and test MAX
       !IF( i > maxve) CALL runend('IMCOOR: Max. number of nodes EXCEDED !!!!!!')
       DO i=1,np                          !loop over all the nodes
         n = n+1
         READ(nf)nn,x(1:nv)               !node label and coordinates
         IF( displ ) READ(nf) d(1:ndime)  !displacements
         IF(iwrit==1)WRITE(lures,"(i10,9f10.5)",ERR=9999) nn,x(1:nv)
         IF(ieulr == 1) THEN                 !EULER angles in a 2-D problem
           rm(1:2) = x(3:4)                  !angle (rad) is the third parameter
         ELSE IF(ieulr == 9) THEN            !EULER angles (rad) in a 3-D problem
           rm(1:6) = x(4:9)                  !Euler angles
         END IF
         IF( neulr > ieulr ) rm = 0d0        !

         SELECT CASE (ir)            !select a label
         CASE (1)
           nd%label = nn             !assign read data
         CASE (2)
           nd%label = i  + lab0(is)  !consecutive values
         CASE (3)
           nd%label = nn + lab0(is)  !add a constant value
         END SELECT
         nn = nd%label

         IF( cooro )THEN
           nd%coord(1:ndime) = x(1:ndime)      !orig. coord
           IF( displ )THEN
              nd%coora(1:ndime) = x(1:ndime)+d(1:ndime)      !actual coord
           ELSE
              nd%coora(1:ndime) = x(1:ndime)      !actual coord
           END IF
         ELSE
           nd%coora(1:ndime) = x(1:ndime)      !stage coord
           IF( displ )THEN
              nd%coord(1:ndime) = x(1:ndime)-d(1:ndime)      !actual coord
           ELSE
              nd%coord(1:ndime) = x(1:ndime)   !stage coord
           END IF
         END IF
         nd%coorc(1:ndime) = x(1:ndime)      !actual coord

         IF (neulr > 0) THEN
           nd%euler(1:neu,1) = rm(1:neu)        !initial angles(rad)
           nd%euler(1:neu,2) = rm(neu+1:2*neu)  !present angles(rad)
         END IF
         nd%veloc(1:ndofn) = 0d0             !initializes velocity
         IF (therm) nd%temp = 0d0            !initializes temperature

         CALL l_insert (ndi, nd)             !add node to list of nodes
         CALL add_at_end (ns, nn)  !add node to SET
       END DO

       CALL add_ns_at_end (ns)         !add set of nodes to the list of sets
     END IF
   END DO

   npoin = n                           !new number of nodes
              ! get memory for NPO_DB
   IF( (npoin <= max_npoin_3d .AND. ndime == 3 ) .OR. (npoin <= max_npoin_2d .AND. ndime == 2 )) THEN
     ALLOCATE (coord(ndime,npoin), coora(ndime,npoin), coorc(ndime,npoin),  &
               label(npoin), velnp(ndofn,npoin))
   ELSE
      STOP ' number of nodes exceeded for Academic version'
   END IF
   IF (therm) THEN ! allocation of nodal vectors for thermal analysis
     ALLOCATE ( tempe(1:ndoft,npoin), dtemp(1:ndoft,npoin) )
     tempe = 0d0 ; dtemp = 0d0
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
       !CALL inrotm(nd%euler(1:neu,2),euler(1:9,n))
     END IF
     IF (therm) tempe(:,n) = nd%temp(1,1:ndoft)
     IF (therm) dtemp(:,n) = nd%temp(2,1:ndoft)
     CALL l_next (ndi)                 !next node in the list
   END DO
   CALL dalloc_list_ndata (ndi)  ! deallocate all the temporary database

   IF (TRIM(actio) == 'NEW') THEN   !for a NEW problem
     npoio = npoin               !old values are present values
     ALLOCATE (oldlb(npoio))
     DO i=1,npoio
       oldlb(i) = label(i)
     END DO
     actio = 'NSTRA2'
   ELSE IF( chlab )THEN  !for subsequent strategies
     actio = 'NSTRA2'
     ! updating local node numbers in existing problem
     echo_chnode = .FALSE.
     CALL elemnt ('UPDLON')               ! elements/segments/nodes
     CALL updlon_kc (naris,nvelr,ndime)   ! KINEMATIC conditions
     IF(nconm)CALL updlon_cm (ndofn)      ! CONCENTRATED MASSES
     CALL updlon_outp (oldlb)             ! output nodes
     CALL updlon_damp ()                  ! DAMPING nodes
     CALL updlon_lo (nload,oldlb)         ! LOADS
     CALL contac ('UPDLON',0)             ! CONTACT nodes (CONT 2 & 3 & 4 only)
     echo_chnode = .TRUE.
     DEALLOCATE (oldlb)
     ALLOCATE (oldlb(npoin))
     DO i=1,npoio
       oldlb(i) = label(i)
     END DO
   END IF

 RETURN
  9999 CALL runen2('')
  9000 FORMAT(' imported node as OVERWRITE does not exists, ELM_SET ',a,' NODE ',i7)
 END SUBROUTINE imcoor
