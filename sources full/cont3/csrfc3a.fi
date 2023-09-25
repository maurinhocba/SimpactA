 SUBROUTINE csrfc3a(ipair,coora,disma,ffdis,pair,surf, &
                    dtime,emass,fcont,cmpfr,td)

 !.... searching procedure for a given self contac 2-d surface

 IMPLICIT NONE
 !     Dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ipair  !pair position
 REAL (kind=8), INTENT(IN) :: disma,   &!10 times the max. displ. expected
                              dtime,   &!time increment
                              ffdis,   &!displacement increment factor
                              td,      &!thermal dtime
                              emass(:,:), &! nodal masses
                              coora(:,:)   ! coordinates of surface
 REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)  !contact forces
 TYPE (pair3_db), POINTER :: pair             !contact pair
 TYPE (surf3_db), POINTER :: surf             !contact surface
 LOGICAL, INTENT(IN) :: cmpfr
 !     local variables
 INTEGER (kind=4) ncnod,nsegm
 REAL (kind=8) cprop(10),maxdi
 INTEGER (kind=4), POINTER :: lcseg(:,:),nhseg(:,:)

 !.... set parameters
 nsegm = surf%nsegm       !No of nodes per segment in surface
 ncnod = pair%ncnod       !No of nodes in surface
 !.... recall pair contac parameters
 cprop(1) = pair%npenal    !normal penalty parameter
 cprop(2) = pair%tpenal    !tangential penalty parameter
 cprop(3) = pair%static    !Static friction coefficient
 cprop(4) = pair%cutoff    !maximum penetration expected
 cprop(5) = pair%gapinc    !maximum incremental penetration expected
 cprop(6) = pair%kinet     !Kinetic friction coefficient
 cprop(7) = pair%hcont     ! h factor for heat conduction
 cprop(8) = pair%econt     ! e exponent for heat conduction
 cprop(9) = pair%hv        ! Vicker's hardness
 cprop(10)= pair%rec       ! Relative efusivity (slave)
 !IF( .NOT. cmpfr ) cprop(4) = 3d0*cprop(4)   !maximum incremental penetration expected

 IF( pair%mtsur < 0)THEN ! if surface uses the bottom surface
   lcseg => surf%lcseb
   nhseg => surf%nhseb
 ELSE                    ! else uses top surface
   lcseg => surf%lcseg
   nhseg => surf%nhseg
 END IF

 !.... Recompute coordinates of center of triangles if necessary
 IF( surf%cxc ) CALL cxctc3(coora,surf)

 maxdi = (REAL(pair%freq)/5d0*disma)

 !.... perform contact searching and compute contact forces

 CALL cscf3aa(lcseg,nhseg,ncnod,surf%lcnod,surf%xc,surf%tc,    &
              pair%issdb,pair%rssdb,coora,disma,ffdis,maxdi,   &
              fcont,surtf(:,ipair), emass,cprop,               &
              dtime,nsegm,pair%freq,surf%cxc,                  &
              surf%tn,pair%press,pair%presn,                   &
              pair%wrink,pair%mingp,cmpfr,                     &
              surf%area,pair%mtdof,pair%sldof,td)
 surf%cxc = .FALSE.     !set to recomputed

 RETURN
 END SUBROUTINE csrfc3a
