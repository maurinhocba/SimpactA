 SUBROUTINE cscf3a(lcseg,nhseg,ncnod,lcnod,xc,            &
                   issdb,rssdb,coors,coorm,disma,ffdis,   &
                   fcont,surtf,emass,indco,cprop,         &
                   dtime,nsegm,freq,                      &
                   cursl,curms,tn,cu,press,presn,         &
                   wear,wwear,wrink,mingp,cmpfr,          &
                   area,mtdof,sldof,td)

 ! perform search for contact interface
 ! and compute contact forces for a given master/slave pair
 !
 !  PARALLEL VERSION
 !  commented with different number of !
 !  !! are standard parallelized loop
 !  !!! for parallelized sections
 !  !!!! not checked for large numbers of nodes

 ! !$ USE omp_lib
 USE gvar_db, ONLY : updiv
 IMPLICIT NONE
 !   Dummy arguments
 LOGICAL, INTENT(IN) :: cursl,curms,press,wear,wrink,cmpfr
 INTEGER (kind=4), INTENT(IN) :: ncnod,indco,nsegm,freq,   &
                   lcseg(:,:),nhseg(:,:),lcnod(:),mtdof,sldof
 INTEGER (kind=4), INTENT(IN OUT) :: issdb(:,:)
 REAL (kind=8), INTENT(IN) :: coors(:,:),coorm(:,:),emass(:,:),disma,  &
                          ffdis,dtime,cprop(:),xc(:,:),tn(:,:),cu(:,:),td
 REAL (kind=8), POINTER :: area(:)
 REAL (kind=8), INTENT(IN OUT) :: rssdb(:,:),fcont(:,:),surtf(:),      &
                                  presn(:),wwear(:),mingp(:)

 !  local variables
 INTEGER (kind=4) :: icnod,np,ia,nearn,isn,nnp,npc,order(ncnod),n,node,first,every
 REAL (kind=8) :: prdat(6),eresf(3,0:3),lwear(0:3)

 !Auxiliar variables for parallel code
 !!!! !$ INTEGER(kind=4), PARAMETER :: maxno=10240,minno= 8192    !2 threads
 !!!! !$ INTEGER(kind=4):: nloop,iloop,irest
 ! !$ INTEGER(kind=4),PARAMETER:: MAXLOCK=2048               !number of keys
 ! !$ INTEGER(kind=4),SAVE :: lock1(0:MAXLOCK-1),lock2 !keys
 ! !$ LOGICAL,SAVE         :: inilock=.FALSE.           !initialization flag

 !Initialize exclusion variables for parallel code (only first time)
 ! !$ IF (.NOT.inilock) THEN                   !first time only
 ! !$   DO icnod=0,MAXLOCK-1                       !for each key
 ! !$     CALL omp_init_lock(lock1(icnod))          !unlock
 ! !$   END DO
 ! !$   CALL omp_init_lock(lock2)              !unlock
 ! !$   inilock = .TRUE.                         !change flag
 ! !$ END IF

 !first loop to identify nodes to be checked
 npc = 0               !initializes nodes to be checked
 order = 0

 first = 1
 every = 1

 !!! !$OMP PARALLEL                                                &
 !!! !$OMP SHARED(npc,order,ffdis,ncnod,lcnod,issdb,freq)          &
 !!! !$OMP PRIVATE(ia,icnod,isn,nearn,nnp,np,first,every)
 !!! !$ first = 1+omp_get_thread_num()
 !!! !$ every = omp_get_num_threads()

 ! !$OMP DO  &                                                                              !!
 ! !$OMP PRIVATE(icnod,isn,nearn,nnp,np,ia)                                                 !!

 DO icnod=first,ncnod,every

   !... step 1 ... look for the nearest segment

   isn   = lcnod(icnod)   ! global number of the current slave node
   nearn = issdb(1,icnod) ! previous nearest segment
   nnp   =-issdb(3,icnod) ! previous time of effective penetration

   IF( nearn < 0 )THEN           !If a search is advisable
     nearn = -nearn
     IF(nearn > 1000000 )THEN        !proj. seems not necessary yet
       np = nearn / 1000000          !remaining steps until proj.
       ia = ffdis*np                   !updates remaining steps
       nearn = MOD(nearn, 1000000) + ia * 1000000  !composed value
       IF(nearn > 1000000)THEN
         issdb(1,icnod) = -nearn + 1000000      !update nnear
         issdb(3,icnod) = issdb(3,icnod) - 1    !increase counter
         CYCLE !check another node
       ELSE
         issdb(1,icnod) = -nearn                !update nnear
       END IF
     END IF
   END IF
 !   !$ CALL omp_set_lock(lock2)    !Start blocking access in shared variables
   npc = npc + 1
   order(npc) = icnod
 !   !$ CALL omp_unset_lock(lock2)  !End blocking access in shared variables

 END DO
 ! !$OMP END DO                                                                             !!
 !!! !$OMP END PARALLEL

!.... loop over all slave nodes to be checked
 first = 1
 every = 1

 ! !$OMP PARALLEL                                                                   &
 ! !$OMP SHARED(npc,order,lcnod,issdb,freq,xc,nsegm,coors,coorm,lcseg,rssdb,nhseg,  &
 ! !$OMP        disma,cprop,curms,cu,wrink,mingp,cmpfr,indco,dtime,emass,cursl,tn,  &
 ! !$OMP        press,wear,mtdof,sldof,area,td,presn,surtf,fcont,updiv,wwear)       &
 ! !$OMP PRIVATE(ia,first,every,icnod,isn,nearn,nnp,prdat,eresf,lwear,n,node)
 ! !$ first = 1+omp_get_thread_num()
 ! !$ every = omp_get_num_threads()
 !!!! this (next) version does not work
 !!!! !$ ia      = (npc+omp_get_num_threads()-1)/omp_get_num_threads()
 !!!! !$ first   = MIN(npc,1+omp_get_thread_num()*ia)     !Find the first element in the list
 !!!! !$ every   = MIN(npc,ia*(1+omp_get_thread_num()))   !Find the last element in the list

 !! !$OMP DO  &  !the other version is faster
 !! !$OMP PRIVATE(ia,icnod,isn,nearn,nnp,prdat,eresf,lwear,n,node)

 DO ia=first,npc,every
 !!!! DO ia=first,every !desagregates by zones

   icnod = order(ia)
   !... step 1 ... look for the nearest segment

   isn   = lcnod(icnod)   ! global number of the current slave node
   nearn = issdb(1,icnod) ! previous nearest segment
   nnp   =-issdb(3,icnod) ! previous time of effective penetration

   IF(nearn == 0 .OR. nnp > freq ) THEN   !if no previous or too many steps
     !global search in ALL master segment the nearest segment to ISN
     CALL nears3(xc,nsegm,nearn,coors(:,isn))
     issdb(3,icnod) = 0   ! re-initializes counter

   ELSE IF( nearn < 0)THEN           !If a search is advisable
     nearn = -nearn
   END IF

   !.... step 2 .... project slave node onto the master surface

   !....   three-dimensional 3-node master segments
   CALL projt3(coors(:,isn),nearn,coorm,lcseg,rssdb(:,icnod),   &
        nhseg,issdb(:,icnod),prdat,disma,cprop(4),cprop(5),curms,cu)
   IF( wrink )THEN
     IF( rssdb(3,icnod) < mingp(icnod) )THEN
       IF( rssdb(3,icnod) > 0d0 .OR. prdat(1) < 0d0) &
         mingp(icnod) = rssdb(3,icnod)
     END IF
   END IF

   IF( prdat(1) == 0d0 )CYCLE   ! no penetration
   IF( cmpfr )THEN

     !.... step 3 .... Compute contact force

     CALL force3(isn,issdb(:,icnod),rssdb(:,icnod),lcseg,indco,cprop, &
                 coorm,prdat,dtime,emass,eresf,cursl,tn,              &
                 icnod,press,wear,lwear,mtdof,sldof,area,td)

     IF( press ) presn(icnod) = presn(icnod) + prdat(1)

 !     !$ CALL omp_set_lock(lock2)    !Start blocking access in shared variables
     surtf = surtf - eresf(1:3,0)*dtime
 !     !$ CALL omp_unset_lock(lock2)  !End blocking access in shared variables

     ! Add forces into global contact forces
     IF (indco == 0 .OR. indco == 1) &             !forces on SLAVE
       fcont(1:3,isn) = fcont(1:3,isn) - eresf(1:3,0)  !
     IF (indco == 0 .OR. indco == 2) THEN          !forces on MASTER
       DO n=1,3                 !for each node
         node = lcseg(n,nearn)  !global number
 !         !$ CALL omp_set_lock(lock1(MOD(node,MAXLOCK)))    !Start blocking access in shared variables
         fcont(1:3,node) = fcont(1:3,node) - eresf(1:3,n)
 !         !$ CALL omp_unset_lock(lock1(MOD(node,MAXLOCK)))  !End blocking access in shared variables
       END DO
     END IF

     ! add to friction work (wearing measure)
     IF (wear .AND. issdb(3,icnod) == 2 .AND. updiv ) THEN
       wwear(isn) = wwear(isn) + lwear(0)
       DO n=1,3                 !for each node
         node = lcseg(n,nearn)  !global number
 !         !$ CALL omp_set_lock(lock1(MOD(node,MAXLOCK)))    !Start blocking access in shared variables
         wwear(node) = wwear(node) + lwear(n)
 !         !$ CALL omp_unset_lock(lock1(MOD(node,MAXLOCK)))  !End blocking access in shared variables
       END DO
     END IF

   ELSE
     rssdb(3,icnod) = prdat(1)
   END IF

  END DO
 !! !$OMP END DO
 ! !$OMP END PARALLEL                                                                         !!!
 RETURN
 END SUBROUTINE cscf3a
