 SUBROUTINE csrfc1(ipair,coord,pair,master,slave, &
                   dtime,emass,fcont,velnp,cmpfr)

 !.... searching procedure for a given master-slave skin-tube 3-d  pair

 IMPLICIT NONE
 !     Dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ipair  !pair position
 REAL (kind=8), INTENT(IN) :: dtime,      & !time increment
                              emass(:,:), & ! nodal masses
                              coord(:,:), & ! coordinates
                              velnp(:,:)    ! velocities
 REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)  !contact forces
 TYPE (pair1_db), POINTER :: pair             !contact pair
 TYPE (surf1_db), POINTER :: master,slave     !contact surfaces
 LOGICAL, INTENT(IN) :: cmpfr
 !     local variables
 INTEGER (kind=4) ncnod,nsegm,indco
 REAL (kind=8) cprop(6)
 INTEGER (kind=4), POINTER :: lcseg(:,:)
 !REAL (kind=8), POINTER :: tn(:,:),cu(:)
 !LOGICAL :: cursl, curms

 !.... set parameters
 nsegm = master%nsegm     !No of nodes per segment in master surface
 ncnod = pair%ncnod       !No of nodes in slave surface
 !.... recall pair contac parameters
 indco = pair%indcon       !type of contact
 cprop(1) = pair%npenal    !normal penalty parameter
 cprop(2) = pair%tpenal    !tangential penalty parameter
 cprop(3) = pair%static    !Static friction coefficient
 cprop(4) = pair%kinet     !Kinetic friction coefficient
 cprop(5) = master%diam    !skin diameter
 cprop(6) = slave%diam     !tube diameter

 lcseg => master%lcseg

 !cursl = slave%curved
 !IF( cursl )tn => slave%tn
 !curms = master%curved
 !IF( curms )cu => master%cu

 !.... perform contact searching and compute contact forces
 !... . for a given master-slave pair

 CALL cscf1a(lcseg,ncnod,slave%lcnod,pair%issdb,pair%rssdb,coord,   &
             fcont,surtf(:,ipair),emass,indco,cprop,dtime,nsegm,    &  !cursl,curms,tn,cu,
             pair%press,pair%presn,velnp,cmpfr)

 RETURN
 END SUBROUTINE csrfc1
