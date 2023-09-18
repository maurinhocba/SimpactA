 SUBROUTINE cscf3aa(lcseg,nhseg,ncnod,lcnod,xc,tc,        &
                   issdb,rssdb,coora,disma,ffdis,maxdi,   &
                   fcont,surtf,emass,cprop,               &
                   dtime,nsegm,freq,cxc,            &
                   tn,press,presn,                        &
                   wrink,mingp,cmpfr,                     &
                   area,mtdof,sldof,td)

 ! perform search for contact interface
 ! and compute contact forces for a given self-contact surface
 ! NOT PARALLELIZED, only adapted to use new version of FORCE3 routine

 ! !$ USE omp_lib
 USE gvar_db, ONLY : updiv
 IMPLICIT NONE
 !   Dummy arguments
 LOGICAL, INTENT(IN) :: press,   & !map of nodal pressure required
                        wrink,   & !minimum gap information required
                        cxc,     & !surface geometry recomputed
                        cmpfr      !compute contact forces
 INTEGER (kind=4), INTENT(IN) :: ncnod,     & !number of nodes in the surface
                                 nsegm,     & !number of segments in the surface
                                 freq,      & !frequency of geometry update
                                 lcseg(:,:),& !segment connectivities
                                 nhseg(:,:),& !neighbour segments
                                 lcnod(:),  & !list of nodes
                                 mtdof,sldof
 INTEGER (kind=4), INTENT(IN OUT) :: issdb(:,:)   !integer data-base
 REAL (kind=8), INTENT(IN) :: coora(:,:),   & !present coordinates
                              emass(:,:),   & !nodal mass
                              disma,        & !maximum displacement in 10 steps
                              tc(:,:),      & !outward normals at segment centers
                              ffdis,        & !factor to update distance
                              td,           & !thermal dtime
                              dtime,        & !time increment for penalty (forces)
                              cprop(:),     & !properties of contact
                              xc(:,:),      & !segment centers
                              tn(:,:),      & !outward normals at nodes
                              maxdi           !maximum displacement between updates
 REAL (kind=8), POINTER :: area(:)
 REAL (kind=8), INTENT(IN OUT) :: rssdb(:,:), & !real data base
                                  fcont(:,:), & !contact forces
                                  surtf(:),   & !contact forces
                                  presn(:),   & !nodal pressure
                                  mingp(:)      !minimum gaps

 !  local variables
 INTEGER (kind=4) :: icnod,np,ia,nearn,isn,nnp,n,node
 REAL (kind=8) :: prdat(6),           & !contact data at node
                  eresf(3,0:3),       & !patch contact forces
                  lwear(0:3),         & !auxiliar
                  cutoff                !maximum penetration for search
 REAL(kind=8), POINTER :: cu(:,:)


 INTERFACE
   INCLUDE 'nearsa.h'
 END INTERFACE


  cutoff = 3d0*cprop(4)  !maximum penetration to discard nodes in search

  DO icnod=1,ncnod

   !... step 1 ... look for the nearest segment

   isn = lcnod(icnod)         ! global number of the current slave node
   nearn = issdb(1,icnod)     !previous nearest segment
   nnp   =-issdb(3,icnod)     !previous time of effective penetration

   IF(( nearn <= 0 .OR. nnp > freq) .AND. cxc ) THEN   !if no penetration and surface recomputed
     !global search in ALL master segment the nearest segment to ISN
     CALL nearsa(xc,nsegm,nearn,coora(:,isn),tc,tn(:,icnod),maxdi,cutoff)
     IF( nearn == 0 ) THEN                      !no candidate found
       issdb(1,icnod) = 0 !nearn                !update nnear
       issdb(3,icnod) = issdb(3,icnod) - 1      !increase counter
       CYCLE !check another node
     END IF
   ELSE IF( nearn == 0 )THEN  !no candidate cycle
     CYCLE
   ELSE IF( nearn < 0)THEN    !a previous candidate
     nearn = -nearn
     IF(nearn >= 1000000 )THEN        !proj. seems not necessary yet
       np = nearn / 1000000           !remaining steps until proj.
       ia = ffdis*np                  !updates remaining steps
       nearn = MOD(nearn, 1000000) + ia * 1000000  !composed value
       IF(nearn >= 1000000)THEN
         issdb(1,icnod) = -nearn + 1000000      !update nnear
         issdb(3,icnod) = issdb(3,icnod) - 1    !increase counter
         CYCLE !check another node
       END IF
     END IF
   END IF

   !.... step 2 .... project slave node onto the master surface

   !....   two-dimensional 2-node master segments
   CALL projt3(coora(:,isn),nearn,coora,lcseg,rssdb(:,icnod),   &
        nhseg,issdb(:,icnod),prdat,disma,cprop(4),cprop(5),.FALSE.,cu)
   IF( wrink )THEN   !IF wrinkles control
     IF( rssdb(3,icnod) < mingp(icnod) )THEN
       IF( rssdb(3,icnod) > 0d0 .OR. prdat(1) < 0d0) &
         mingp(icnod) = rssdb(3,icnod)
     END IF
   END IF

   IF( prdat(1) == 0d0 )CYCLE   ! no penetration
   IF( cmpfr )THEN

     !.... step 3 .... Compute contact force

     CALL force3(isn,issdb(:,icnod),rssdb(:,icnod),lcseg,   0 ,cprop, &
                 coora,prdat,dtime,emass,eresf,.FALSE.,tn,              &
                 icnod,press,.FALSE.,lwear,mtdof,sldof,area,td)

     !IF nodal pressure map
     IF( press ) presn(icnod) = presn(icnod) + prdat(1)

     surtf = surtf + ABS(eresf(1:3,0))*dtime

     ! Add forces into global contact forces
     fcont(1:3,isn) = fcont(1:3,isn) - eresf(1:3,0)  !
     DO n=1,3                 !for each node
       node = lcseg(n,nearn)  !global number
       fcont(1:3,node) = fcont(1:3,node) - eresf(1:3,n)
     END DO

   ELSE
     rssdb(3,icnod) = prdat(1)
   END IF

 END DO
 RETURN
 END SUBROUTINE cscf3aa
