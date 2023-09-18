 MODULE static_db
   IMPLICIT NONE
   LOGICAL :: sstop,          & !stop if no convergence
              spback=.FALSE., & !elastic springback strategy
              autod             !use automatic damping
   INTEGER (kind=4) :: sstep, & !number of steps
                       sitep, & !actual step
                       nides, & !number of steps to introduce prescribed displacement
                       miter, & !maximum number of iterations
                       nperi, & !minimum number of iteration for each step
                       ipred, & !prediction order 1-3
                       ittol, & !tolerance type
                       pmin, pmax, & !minimum & maximun values for damping period
                       lovec, & !load-velocity control set (+load -velocity)
                       str_type, & !strategy type
                       jstep, & !control step
                       nss,   & ! number of iteration to average
                       nds,   & ! number of steps to perform in velocities
                       pcont, & !control point
                       nldof    !loaded DOF with larger displacement/load ratio
   REAL (kind=8) :: tol,      & !tolerance
                    tolm,     & !maximum tolerance
                    rfnor,    & !residual force norm for comparison
                    fdamp,    & !maximum factor to change damping period
                    ltime,    & !time when standard strategy failed
                    olamb,lambd,llamb,   & !factor for load or prescribed velocities
                    kem,      & !kinetic energy for comparison
                    dn0         !displacement increment in NLDOF
   REAL (kind=8), POINTER :: disax(:,:), & !incremental displacements and load factors
                             locsy(:,:), & !local systems at the beginning of each step
                             smass(:,:), & !(ndofn,npoin) modified nodal mass
                             resid_o(:,:), & !(ndime,npoin) equivalent nodal forces at previous step
                             eqres(:)      !residual forces at each active DOF

 CONTAINS
   SUBROUTINE dump_static (neq,ndime,npoin)
   IMPLICIT NONE
   INTEGER(kind=4), INTENT(IN) :: neq,ndime,npoin

   INTEGER(kind=4) :: i,j

   WRITE(50,ERR=9) sstop,autod,sstep,sitep,nides,miter,nperi,ipred,ittol,pmin,pmax,lovec, &
                   str_type,jstep,nss,nds,pcont,nldof
   WRITE(50,ERR=9) tol,tolm,rfnor,fdamp,ltime,olamb,lambd,llamb,ltime,kem,dn0
   WRITE(50,ERR=9) ((disax(i,j),i=1,neq+1),j=1,3)
   WRITE(50,ERR=9) ((resid_o(i,j),i=1,ndime),j=1,npoin)

   RETURN
   9 CALL runen2('')

   END SUBROUTINE dump_static

   SUBROUTINE rest_static (neq,ndime,npoin)
   IMPLICIT NONE
   INTEGER(kind=4), INTENT(IN) :: neq,ndime,npoin

   INTEGER(kind=4) :: i,j

   READ(51)  sstop,autod,sstep,sitep,nides,miter,nperi,ipred,ittol,pmin,pmax,lovec, &
             str_type,jstep,nss,nds,pcont,nldof
   READ(51)  tol,tolm,rfnor,fdamp,ltime,olamb,lambd,llamb,ltime,kem,dn0
   ALLOCATE (disax(neq+1,3))
   READ(51) ((disax(i,j),i=1,neq+1),j=1,3)
   ALLOCATE (resid_o(ndime,npoin))
   READ(51) ((resid_o(i,j),i=1,ndime),j=1,npoin)
   RETURN
   END SUBROUTINE rest_static

 END MODULE static_db
