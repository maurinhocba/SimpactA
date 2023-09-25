 SUBROUTINE csrfc3(ipair,coors,coorm,disma,ffdis,pair,master,slave, &
                   dtime,emass,fcont,cmpfr,td)

!.... searching procedure for a given master-slave 3-d surfaces pair

 IMPLICIT NONE
!     Dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ipair  !pair position
 REAL (kind=8), INTENT(IN) :: disma,   &!10 times the max. displ. expected
                              dtime,   &!time increment
                              ffdis,   &!displacement increment factor
                              td,      &!thermal dtime
                              emass(:,:), &! nodal masses
                              coors(:,:), &! coordinates of slave surface
                              coorm(:,:)   ! coordinates of master surface
 REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)  !contact forces
 TYPE (pair3_db), POINTER :: pair             !contact pair
 TYPE (surf3_db), POINTER :: master,slave     !contact surfaces
 LOGICAL, INTENT(IN) :: cmpfr
!     local variables
 INTEGER (kind=4) ncnod,nsegm,indco,i
 REAL (kind=8) cprop(10)
 INTEGER (kind=4), POINTER :: lcseg(:,:),nhseg(:,:)
 REAL (kind=8), POINTER :: tn(:,:),cu(:,:)
 LOGICAL :: cursl, curms

!.... set parameters
 nsegm = master%nsegm     !No of nodes per segment in master surface
 ncnod = pair%ncnod       !No of nodes in slave surface
!.... recall pair contac parameters
 indco = pair%indcon       !type of contact
 cprop(1) = pair%npenal    !normal penalty parameter
 cprop(2) = pair%tpenal    !tangential penalty parameter
 cprop(3) = pair%static    !friction coefficient
 cprop(4) = pair%cutoff    !maximum penetration expected
 cprop(5) = pair%gapinc    !maximum incremental penetration expected
 cprop(6) = pair%kinet     !kinetic friction coefficient
 cprop(7) = pair%hcont     ! h factor for heat conduction
 cprop(8) = pair%econt     ! e exponent for heat conduction
 cprop(9) = pair%hv        ! Vicker's hardness
 cprop(10)= pair%rec       ! Relative efusivity (slave)
 !IF( .NOT. cmpfr ) cprop(4) = 3d0*cprop(4)   !maximum incremental penetration expected

 IF( pair%mtsur < 0)THEN ! if master surface uses the bottom surface
   lcseg => master%lcseb
   nhseg => master%nhseb
 ELSE                      ! else use top surface
   lcseg => master%lcseg
   nhseg => master%nhseg
 END IF

!.... Recompute coordinates of center of triangles if necessary (MASTER)
 IF( master%cxc )THEN
   DO i=1,nsegm
     master%xc(:,i) = (coorm(:,lcseg(1,i)) + coorm(:,lcseg(2,i)) + &
                       coorm(:,lcseg(3,i)) )/3d0
   END DO
   master%cxc = .FALSE.     !set to recomputed
 END IF

 cursl = slave%curved
 IF( cursl )tn => slave%tn
 curms = master%curved
 IF( curms )cu => master%cu

!.... perform contact searching and compute contact forces
!... . for a given master-slave pair

 CALL cscf3a(lcseg,nhseg,ncnod,slave%lcnod,master%xc,               &
             pair%issdb,pair%rssdb,coors,coorm,disma,ffdis,         &
             fcont,surtf(:,ipair),emass,indco,cprop,                &
             dtime,nsegm,pair%freq,                                 &
             cursl,curms,tn,cu,pair%press,pair%presn,               &
             wear,wwear,pair%wrink,pair%mingp,cmpfr,                &
             slave%area,pair%mtdof,pair%sldof,td)

 RETURN
 END SUBROUTINE csrfc3
