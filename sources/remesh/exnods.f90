 SUBROUTINE exnods(refset,flag)
 !-----------------------------------------------------------------------
 !     Modifies and Restore nodal coordinates
 !-----------------------------------------------------------------------
 USE param_db,ONLY: mnam,max_npoin_2d,max_npoin_3d
 USE ctrl_db,ONLY: ndime, ndofn, neulr, nload, npoin, npoio, neq
 USE esets_db,ONLY: esets, eset, add_name
 USE lispa0,ONLY:  lures
 USE ndinf_db
 USE npo_db
 USE nsets_db
 USE outp_db
 USE meshmo_db
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*),INTENT(IN):: refset   !label of the element set to modify
   LOGICAL,INTENT(IN) :: flag

   !--- Local variables
   !     dbase used as auxiliar array to store coordinates
   !INTEGER (kind=4):: i, n, iset, nn, k, chnode, npkept, nold
   INTEGER(kind=4):: i, j, n, nn, lc(1), mdeln, neu
   INTEGER(kind=4),ALLOCATABLE :: lbind(:),lbarr(:)
   !REAL(kind=8)::  x(6)
   REAL(kind=8), POINTER :: cordi(:,:)
   LOGICAL  :: found,oldns
   TYPE(ndata):: nd
   TYPE(list_ndata):: ndi    ! nodal information, temporary database
   TYPE(nset),POINTER:: ns   ! pointer to set of nodes

   INTERFACE
     INCLUDE 'deloln.h'
     INCLUDE 'verifd.h'
     INCLUDE 'nodfnd.h'
   END INTERFACE


   IF (flag) THEN
     cordi => coora ! updated lagrangian
   ELSE
     cordi => coord ! total lagrangian ??
   END IF
   ! check points for requested output
   IF( nreqd > 0 )CALL verifd(nreqd,nprqd,numpo,ndime,numpn,cordi,coorn, &
                              nodset,label)
   IF( nreqa > 0 )CALL verifd(nreqa,nprqa,numpo,ndime,numpn,cordi,coorn, &
                              nodset,label)
   IF( nreqv > 0 )CALL verifd(nreqv,nprqv,numpo,ndime,numpn,cordi,coorn, &
                              nodset,label)
   IF( nreql > 0 )CALL verifd(nreql,nprql,numpo,ndime,numpn,cordi,coorn, &
                              nodset,label)
   IF( nreqc > 0 )CALL verifd(nreqc,nprqc,numpo,ndime,numpn,cordi,coorn, &
                              nodset,label)
   NULLIFY(cordi)
   !IF( ncdis > 10)THEN
   !  nn = 1
   !  list(1) = ncdis/10
   !  n = MOD(ncdis,10)
   !  CALL verifd(nn,list,numpo,ndime,numpn,coord,coorn, &
   !                      nodset,label)
   !  IF( nn < 0 )ncdis = 10*list(1)+n
   !END IF

   ! store previous nodal information into dbase NDI

   CALL ln_ini(ndi)  ! initialize ndi dbase

   neu = 2*ndime - 3
   DO n = 1, npoin           !for each previous point
     nd%label = label(n)     !node label
     IF (flag) THEN   !solid
       nd%coord(1:ndime) = coora(1:ndime,n)  ! coora is used as coord after rmesh
     ELSE
       nd%coord(1:ndime) = coord(1:ndime,n)  ! coord is used as coord after rmesh
     END IF
     nd%coorc(1:ndime) = coora(1:ndime,n)  ! coora is used as last converged
     nd%coora(1:ndime) = coora(1:ndime,n)  ! coora is coora of course
     IF (neulr > 0) THEN
      IF( neulr > 0 )THEN
        nd%euler(1:neu,1) = eule0(1:neu,n)  !initial euler angles
        IF(ndime == 2)THEN
          nd%euler(1,2) = euler(1,n)        !local systems
        ELSE IF( ndime == 3 ) THEN
          CALL angeul(euler(1:3,n),nd%euler(1:neu,2),.TRUE.) !regenerates Euler angles (in degrees)
        END IF
      END IF
     END IF
     nd%veloc(1:ndofn) = velnp(1:ndofn,n)  ! velocities
     CALL add_nd_at_end(ndi, nd)           ! add to list
   END DO

   npoio = npoin                           ! keep previous number of nodes
   DEALLOCATE(oldlb)                       ! release memory of OLD LABELS
   ALLOCATE(oldlb(npoio))                  !
   oldlb = label                           ! copy LABELS into OLD LABELS
   ! release memory of previous nodal arrays
   DEALLOCATE( coora, coorc, coord, coors, resid, euler, velnp,     &
               label )
   !    delete old nodes in remeshed part from list NDI
   CALL deloln(iwrit, numpo, nodset, maxnn, ndi)
   maxnn = MAX(maxnn,n_first)       !numfst = set label of the 1st node
   mdeln = MAXVAL(nodset(:))        !maximum label in deleted list
   IF(r_elm_zone) maxnn = MAX(maxnn,mdeln)
   maxnn = ((maxnn/1000)+1)*1000    !I prefer it this way (FF)
   n = get_npoin (ndi)    !number of remaining nodes
   WRITE (lures,"(/,'Number of nodes kept from the previous strategy: ',i6/ &
                &   70('-'),/)",ERR=9999) n

   ! add the new set to the list of sets
   CALL nsdb_search (refset, oldns, ns)  !search in the list of sets
   IF( .NOT.oldns )THEN
     CALL ns_ini(ns)                  !initializes list of nodes (set)
     CALL put_nsname(ns, refset)      !assign name set
     CALL add_name (refset,1)         !add to names list
   END IF

   IF (iwrit == 1) WRITE(lures,"(//'   ECHO of the Coordinate Data ',/ &
      & 5x,' Node      X',14X,'Y',14X,'Z',14X,'Alfa      Beta      Gamma')",ERR=9999)

   IF(r_elm_zone)THEN
     ALLOCATE(nodlb(2,nelbnd)) ! nodes labels in zone remeshing
     nodlb = 0 !initialize
   END IF
   j = 0
   DO i = 1, numpn     !for each new node
     n = n + 1         !increase number of total nodes (remaining + new)
     ! to allow us a control for a limited version
     IF( (n >= max_npoin_3d .AND. ndime == 3 ) .OR. (n >= max_npoin_2d .AND. ndime == 2 )) &
       STOP ' number of nodes exceeded for Academic version'
     nn = maxnn + i    !label for new nodes
     IF(r_elm_zone)THEN  !restore old label in boundary of remeshed zone
       CALL nodfnd(lnbnd,nlbnd,cdbnd,oldlb,nn,coorn(:,i),found)
       IF(found)THEN ! if node pertain to zone boundary
         j = j + 1
         nodlb(:,j) = (/(maxnn+i),nn/) ! save this (new and old label)
       END IF
     END IF
     IF (iwrit == 1) WRITE(lures,"(i10,3e15.5)",ERR=9999) nn, coorn(1:ndime,i)
     CALL add_at_end(ns, nn)  !add node to SET

     nd%label = nn                           !node label
     nd%coord(1:ndime) = coorn(1:ndime,i)    !original coordinates
     IF (flag) THEN
       nd%coora(1:ndime) = coorn(1:ndime,i)  !all the same
     ELSE
       nd%coora(1:ndime) = coran(1:ndime,i)  !different
     END IF
     nd%coorc(1:ndime) = nd%coora(1:ndime)   !coorc = coora as usual
     IF (neulr > 0)  nd%euler = 0d0
     nd%veloc(1:ndofn) = v_new(1:ndofn,i)    !interpolated velocities
     CALL l_insert (ndi, nd)                 !add node to list of nodes (FF)
     !CALL add_nd_at_end(ndi, nd)            !add node to list
   END DO
   IF(.NOT.oldns)CALL add_ns_at_end(ns)         !add set of nodes to the list of sets
   npoin = n                                 !update number of nodes

   !      get memory for new nodal arrays
   ALLOCATE(coord(ndime,npoin), label(npoin), coors(ndime,npoin),          &
            coora(ndime,npoin), coorc(ndime,npoin), resid(ndofn,npoin),    &
            velnp(ndofn,npoin))
   velnp = 0d0                 !initializes velocities
   resid = 0d0                 !initializes internal forces

   IF (neulr > 0) THEN           !for nodal systems
     DEALLOCATE(eule0)
     ALLOCATE(eule0(neu,npoin),euler(neulr,npoin))
   ELSE
     ALLOCATE( euler (1,1) )
   END IF

   IF(r_elm_zone)THEN
     ALLOCATE(lbarr(npoin),lbind(npoin)) ! reserve memory
     CALL l_head(ndi)            !point to first point in the list
     DO n=1,npoin                !for each point
       nd = get_ndata(ndi)       !associate data to auxiliar ND
       lbarr(n) = nd%label       !node label
       CALL l_next(ndi)          !point to next point in NDI list
     END DO
     CALL sortlb(npoin,lbarr,lbind) ! node label sorting and indexing
   END IF

   CALL l_head(ndi)            !point to first point in the list
   DO n=1,npoin                !for each point
     nd = get_ndata(ndi)       !associate data to auxiliar ND
     j = n                     !original index
     IF(r_elm_zone)THEN
       lc = MAXLOC(lbind,MASK= lbind == n)  !sorted index for
       j  = lc(1)                           !zone remeshing
     END IF
     label(j) = nd%label                        !node label
     coord(1:ndime,j) = nd%coord(1:ndime)       !coordinates
     coorc(1:ndime,j) = nd%coorc(1:ndime)
     coors(1:ndime,j) = nd%coorc(1:ndime)       ! discuss it!!!!!
     coora(1:ndime,j) = nd%coora(1:ndime)
     velnp(1:ndofn,j) = nd%veloc(1:ndofn)       !velocities
     IF (neulr > 0) THEN
       eule0(1:neu,n) = nd%euler(1:neu,1)    !initial angles
       euler(1:neulr,n) = nd%euler(1:neu,2)  !present angles
     END IF
     CALL l_next(ndi)                      !point to next point in NDI list
   END DO

   ! check points for requested output
   IF (nreqd < 0 .OR. r_elm_zone) CALL verife(nreqd,nprqd,maxnn,r_elm_zone,nodlb,oldlb)
   IF (nreqa < 0 .OR. r_elm_zone) CALL verife(nreqa,nprqa,maxnn,r_elm_zone,nodlb,oldlb)
   IF (nreqv < 0 .OR. r_elm_zone) CALL verife(nreqv,nprqv,maxnn,r_elm_zone,nodlb,oldlb)
   IF (nreql < 0 .OR. r_elm_zone) CALL verife(nreql,nprql,maxnn,r_elm_zone,nodlb,oldlb)
   IF (nreqc < 0 .OR. r_elm_zone) CALL verife(nreqc,nprqc,maxnn,r_elm_zone,nodlb,oldlb)
   !IF (ncdis < 0) THEN
   !  ncdis = -ncdis
   !  nn = ncdis/10
   !  n = MOD(ncdis,10)
   !  ncdis = 10*(nn+maxnn)+n
   !END IF

   IF(r_elm_zone) DEALLOCATE(lbarr,lbind) ! reserve memory
   CALL dalloc_list_ndata(ndi)  ! deallocate all the temporary database

 RETURN
  9999 CALL runen2('')
 END SUBROUTINE exnods
