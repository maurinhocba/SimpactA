 SUBROUTINE cscf1a(lcseg,ncnod,lcnod,issdb,rssdb,coord,          &
                   fcont,surtf,emass,indco,cprop,dtime,nsegm,    &  !cursl,curms,tn,cu,
                   press,presn,velnp,cmpfr)

 ! perform search for contact interface
 ! and compute contact forces for a given master/slave pair

 IMPLICIT NONE
 !   Dummy arguments
 LOGICAL, INTENT(IN) :: press,cmpfr  !cursl,curms,
 INTEGER (kind=4), INTENT(IN) :: ncnod,indco,nsegm,   &
                                 lcseg(:,:),lcnod(:)
 INTEGER (kind=4), INTENT(IN OUT) :: issdb(:,:)
 REAL (kind=8), INTENT(IN) :: coord(:,:),emass(:,:),dtime,cprop(:),velnp(:,:)  !,xc(:,:),tn(:,:),cu(:)
 REAL (kind=8), INTENT(IN OUT) :: rssdb(:,:),fcont(:,:),surtf(:),      &
                                  presn(:)

 !  local variables
 INTEGER (kind=4) :: icnod,nearn,justn,isn
 REAL (kind=8) :: prdat(9),dd

 justn = 1
 !.... loop over all slave nodes
 dd = (cprop(5) - cprop(6))/2d0   !half the diameter difference
 DO icnod = 1, ncnod

   !... step 1 ... look for the nearest segment

   isn = lcnod(icnod) ! global number of the current slave node
   nearn = issdb(1,icnod) !previous nearest segment
   !nnp   =-issdb(3,icnod) !previous time of effective penetration

   IF( nearn == 0 ) THEN   !if no previous
     !global search in ALL master segment the nearest segment to ISN
     nearn = justn
     CALL nearst(lcseg,nsegm,nearn,coord(:,isn),coord)
     justn = nearn

   END IF

   !.... step 2 .... project slave node onto the master surface

   !....   three-dimensional 2-node master segments
   CALL projt1(coord(:,isn),nearn,coord,lcseg,rssdb(:,icnod), &
               issdb(:,icnod),prdat,dd,nsegm)     !,curms,cu)

   IF( prdat(1) >= 0d0 )CYCLE   ! no penetration
   IF( cmpfr )THEN

     !.... step 3 .... Compute contact force

     CALL force1(isn,issdb(:,icnod),rssdb(:,icnod),lcseg,indco,cprop, &
                 coord,prdat,dtime,emass,fcont,velnp,surtf,press)  !cursl,tn

     IF( press ) presn(icnod) = presn(icnod) + prdat(1)
   ELSE
     rssdb(2,icnod) = prdat(1)
   END IF

 END DO
 RETURN
 END SUBROUTINE cscf1a
