 SUBROUTINE cscf2a(lcseg,nhseg,ncnod,lcnod,xc,            &
                   issdb,rssdb,coors,coorm,disma,ffdis,   &
                   fcont,surtf,emass,indco,cprop,         &
                   dtime,nsegm,freq,                      &
                   cursl,curms,tn,cu,press,presn,         &
                   wear,wwear,wrink,mingp,cmpfr,          &
                   area,mtdof,sldof,td)

 ! perform search for contact interface
 ! and compute contact forces for a given master/slave pair

 IMPLICIT NONE
 !   Dummy arguments
 LOGICAL, INTENT(IN) :: cursl,curms,press,wear,wrink,cmpfr
 INTEGER (kind=4), INTENT(IN) :: ncnod,indco,nsegm,freq,   &
                   lcseg(:,:),nhseg(:,:),lcnod(:),mtdof,sldof
 INTEGER (kind=4), INTENT(IN OUT) :: issdb(:,:)
 REAL (kind=8), INTENT(IN) :: coors(:,:),coorm(:,:),emass(:,:),disma,  &
                          ffdis,dtime,cprop(:),xc(:,:),tn(:,:),cu(:),td
 REAL (kind=8), POINTER :: area(:)
 REAL (kind=8), INTENT(IN OUT) :: rssdb(:,:),fcont(:,:),surtf(:),      &
                                  presn(:),wwear(:),mingp(:)

 !  local variables
 INTEGER (kind=4) :: icnod,np,ia,nearn,isn,nnp
 REAL (kind=8) :: prdat(6)

 !.... loop over all slave nodes
 DO icnod = 1, ncnod

   !... step 1 ... look for the nearest segment

   isn = lcnod(icnod) ! global number of the current slave node
   nearn = issdb(1,icnod) !previous nearest segment
   nnp   =-issdb(3,icnod) !previous time of effective penetration

   IF(nearn == 0 .OR. nnp > freq ) THEN   !if no previous or too many steps
     !global search in ALL master segment the nearest segment to ISN
     CALL nearst(xc,nsegm,nearn,coors(:,isn))
     issdb(3,icnod) = 0   ! re-initializes counter

   ELSE IF( nearn < 0)THEN           !If a search is advisable
     nearn = -nearn
     IF(nearn > 1000000 )THEN        !proj. seems not necessary yet
       np = nearn / 1000000          !remaining steps until proj.
       ia = ffdis*np                   !updates remaining steps
       nearn = MOD(nearn, 1000000) + ia * 1000000  !composed value
       IF(nearn > 1000000)THEN
         issdb(1,icnod) = -nearn + 1000000      !update nnear
         issdb(3,icnod) = issdb(3,icnod) - 1    !increase counter
         CYCLE !check another node
       END IF
     END IF
   END IF

   !.... step 2 .... project slave node onto the master surface

   !....   two-dimensional 2-node master segments
   CALL projt2(coors(:,isn),nearn,coorm,lcseg,rssdb(:,icnod), &
        nhseg,issdb(:,icnod),prdat,disma,cprop(4),cprop(5),curms,cu)
   IF( wrink )THEN
     IF( rssdb(2,icnod) < mingp(icnod) )THEN
       IF( rssdb(2,icnod) > 0d0 .OR. prdat(1) < 0d0) &
         mingp(icnod) = rssdb(2,icnod)
     END IF
   END IF

   IF( prdat(1) == 0d0 )CYCLE   ! no penetration
   IF( cmpfr )THEN

     !.... step 3 .... Compute contact force

     CALL force2(isn,issdb(:,icnod),rssdb(:,icnod),lcseg,indco,cprop, &
                 coorm,prdat,dtime,emass,fcont,surtf,cursl,tn,icnod,  &
                 press,wear,wwear,mtdof,sldof,area,td)

     IF( press ) presn(icnod) = presn(icnod) + prdat(1)
   ELSE
     rssdb(2,icnod) = prdat(1)
   END IF

 END DO
 RETURN
 END SUBROUTINE cscf2a
