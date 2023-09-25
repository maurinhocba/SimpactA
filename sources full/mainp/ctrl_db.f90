 MODULE ctrl_db
 ! global control paramters
 USE param_db, ONLY: midn, mttl, mnam
 IMPLICIT NONE

   CHARACTER (len=midn) :: ptype='EXPLCT' ! problem type
   CHARACTER (len=mttl) :: text           ! problem title

   ! size (dimension) parameters
   INTEGER (kind=4) ::  &
             ndime,        & !problem dimension
             ndofn,        & !number of DOFs per node
             nrotd,        & !number of DOFs per node, exclusing extra DOFs
             neulr,        & !if local systems used
             ntype,        & !problem type for 2-D 1:plane stress 2:plane strain 3:axilsymmetric
             nload=0,      & !number of load sets
             nheat,        & !number of heat sets
             npoin=0,      & !number of points in the mesh
             npoio=0         !number of points in previous mesh

   INTEGER (kind=4) :: neq,    & !number of active equations
                       neqold    !number of equations at previous step
   INTEGER (kind=4) :: neqt,   & !number of active temperature equations
                       neqto     !number of temperature equations at prev step
   LOGICAL          :: lumped = .TRUE.  !use lumped mass matrix

   ! stepping control parameters
   INTEGER (kind=4)  :: istep=0,   & ! step number
                        nstra=1,   & ! counter of strategies read from data file
                                     ! it can be istra /= nstra due to remeshing
                        ncdlt        ! Frequency for DT computation
   ! time integration parameters
   REAL (kind=8) :: &
                    dtuser,    & ! Time Increment = 0.-> Auto calc
                    dtscal,    & ! SCAle fact. for DT in auto mode
                    begtm,     & ! Time of Begining of the new strategy
                    endtm,     & ! Final Time
                    dtime,     & ! time step
                    dtold=0d0, & ! "old" time step (in the previous step)
                    ttime=0d0, & ! accumulated time of the process
                    vefac=0d0, & ! velocity smoothing factor for viscous behavior
                    mscal        ! Mass SCALing factor (affecting dtime)

   ! contact/drawbead codes
   INTEGER (kind=4), PARAMETER :: mconts = 4 ! number of impl. contact algors.
                                             ! 2-3: contacts, 1: drawbeads
   INTEGER (kind=4) :: nconp(mconts) = 0,  & ! numbers of contact pairs
                       numct = 0,          & ! number of contact algors. used
                       ndimc = 0,          & ! number of contact DOFs
                       ctype(mconts) = 0     ! contact types (algors.) used
   LOGICAL :: cactive(mconts) = .FALSE.,   & ! active contacs
              bottom = .FALSE.,            & ! bottom surface is used
              top    = .FALSE.               ! top surface is used

   ! ground accelerogram parameters (disabled?)
   INTEGER (kind=4) :: &
                    ifixd,    &  !direction of ground acceler.
                    ifunc,    &  !flag of ground acceleration (0-active)
                    nacce        !number of points in the accelerogram
   REAL (kind=8) :: &
                    dtrec,    &  ! Time increment of Accelerogram
                    tbega,    &  ! Time increment of Accelerogram
                    a0=1d0,   &  ! acceleration factor
                    xbase(3)=0d0 ! base acceleration, velocity & displacement

 !   thermal analysis parameters
   LOGICAL       :: &
                    itemp= .FALSE.,  &  ! compute deformations due to temperature mod
                    therm= .FALSE.,  &  ! THERMal analysis code code
                    ascal= .FALSE.      ! automatic scaling factor for Mass Matrix

   INTEGER (kind=4) :: &
                    ndoft      !number of temperature DOFs per node

   REAL (kind=8) :: &
                    tdtuse,  &  ! Time increment for heat transfer analysis
                    tdtime,  &  ! Time increment for heat transfer analysis
                    tdtsca,  &  ! Time scaling factor for thermal analysis
                    tmscal,  &  ! scaling factor for capacity matrix
                    tref,    &  ! Reference temperature
                    tscal       ! Time scaling factor for thermal analysis
 
 !   aerodynamical analysis
   LOGICAL :: aero

   ! other
   LOGICAL :: echo_chnode = .TRUE.
   LOGICAL :: initial_displacements = .FALSE.
   INTEGER (kind=4)   ::  memo        !size of auxiliar array for input

 CONTAINS

    SUBROUTINE dump_ctrl
    !   dumps database of control parameters
    IMPLICIT NONE

      WRITE (50,ERR=9999) ndime, ndofn, nrotd, neulr, ntype, nload, npoin, npoio
      WRITE (50,ERR=9999) neq, memo, lumped
      WRITE (50,ERR=9999) ptype, text
      WRITE (50,ERR=9999) istep, nstra, ncdlt
      WRITE (50,ERR=9999) dtuser, dtscal, begtm, endtm, dtime, dtold, ttime, vefac, mscal
      WRITE (50,ERR=9999) nconp, numct, ndimc, ctype, cactive, bottom, top
      WRITE (50,ERR=9999) ifixd, ifunc, nacce, dtrec, tbega, a0, xbase
      WRITE (50,ERR=9999) neqt, neqto, itemp, therm,  ndoft, &
                          tdtuse, tdtsca, tmscal, tref, tscal
      WRITE (50,ERR=9999) aero

    RETURN
    9999 CALL runen2('')
    END SUBROUTINE dump_ctrl

    SUBROUTINE rest_ctrl

      !   restores database of control parameters

      READ (51) ndime, ndofn, nrotd, neulr, ntype, nload, npoin, npoio
      READ (51) neq, memo, lumped
      READ (51) ptype, text
      READ (51) istep, nstra, ncdlt
      READ (51) dtuser, dtscal, begtm, endtm, dtime, dtold, ttime, vefac, mscal
      READ (51) nconp, numct, ndimc, ctype, cactive, bottom, top
      READ (51) ifixd, ifunc, nacce, dtrec, tbega, a0, xbase
      READ (51) neqt, neqto, itemp,  therm,  ndoft,  &
                tdtuse, tdtsca, tmscal,  tref, tscal
      READ (51) aero

    END SUBROUTINE rest_ctrl

 END MODULE ctrl_db
