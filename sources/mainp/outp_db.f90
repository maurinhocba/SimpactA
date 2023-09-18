 MODULE outp_db
   !store data to control global results output
   !read at CONTOL and INTIME
   USE param_db,ONLY: mnam
   IMPLICIT NONE
   !*** slist and slave_list created to compute total internal forces
   !*** over master nodes for output with NREQL
   TYPE slist
     INTEGER (kind=4) :: node,dof,pos
     TYPE (slist), POINTER :: next
   END TYPE slist

   TYPE slave_list
     INTEGER (kind=4) :: nvalues
     INTEGER (kind=4), ALLOCATABLE :: deps(:,:) !(3,nvalues) dof,node,pos
     TYPE (slave_list), POINTER :: next
   END TYPE slave_list

   TYPE (slave_list), POINTER :: sl_head,sl_tail
   !***

   CHARACTER (len=mnam) :: cname='' ! curve name (label)

   ! output parameters
   INTEGER (kind=4) :: iwrit ! Print-out code (flag)     -     0:NO 1:Yes
   INTEGER (kind=4) :: iener, & !to print energy values, old, not in use
                       nreqa=0, & !number of nodes with accelerations
                       nreqc=0, & !number of nodes with contact forces
                       nreqd=0, & !number of nodes with displacements
                       nreql=0, & !number of nodes with internal forces
                       nreqv=0, & !number of nodes with velocities
                       nreqp=0, & !number of nodes with pressures
                       nreqt=0, & !number of nodes with temperatures
                       ncdis=0  !point to control output
   INTEGER (kind=4),POINTER :: nprqa(:), & !nodes to print accelerations
                               nprqc(:), & !nodes to print contact forces
                               nprqd(:), & !nodes to print displacements
                               nprql(:), & !nodes to print internal forces
                               nprqv(:), & !nodes to print velocities
                               nprqp(:), & !nodes to print pressures
                               nprqt(:)    !nodes to print temperatures
   REAL (kind=8) :: toutd   !frequency to print nodal values
   REAL (kind=8),POINTER :: energ(:),  & !variables associated to energy
                            toutp(:),  & !parameters to control output
                            res(:,:)     !internal forces at selected points
   CHARACTER (len=1) :: postype = 'T'    !Postprocess by time, displacement or curve value
   LOGICAL :: lastst = .FALSE.
   LOGICAL :: thickc = .FALSE.
   REAL (kind=8) :: thick1=1d0, thick2=1d0
   REAL     (kind=8) :: sumat  ! total mass

   ! timing information
   REAL (kind=8) :: time(40)=0d0,  & ! elapsed time for different tasks
                    timed(20)=0d0, & ! elapsed time for different tasks (development)
                    cpui             ! initial system time (for comparison)

 CONTAINS
   SUBROUTINE dump_outp
   IMPLICIT NONE
   INTEGER(Kind=4) :: i,j

   j = INT(toutp(1)) + 1
   WRITE(50,ERR=9999) iwrit,iener, nreqa, nreqc, nreqd, nreql, nreqv, nreqp, nreqt, ncdis, j, cname
   WRITE(50,ERR=9999) nprqa, nprqc,  nprqd, nprql, nprqv, nprqp, nprqt
   WRITE(50,ERR=9999) energ, sumat, toutd, (toutp(i),i=1,j)
   WRITE (50,ERR=9999) time, timed, cpui
   RETURN
   9999 CALL runen2('')
   END SUBROUTINE dump_outp

   SUBROUTINE rest_outp
   IMPLICIT NONE
   INTEGER(Kind=4) :: i,j

   READ(51) iwrit,iener, nreqa, nreqc, nreqd, nreql, nreqv, nreqp, nreqt, ncdis, j, cname

   ALLOCATE( nprqa(MAX(1,nreqa)), nprqc(MAX(1,nreqc)),           &
             nprqd(MAX(1,nreqd)), nprql(MAX(1,nreql)),           &
             nprqv(MAX(1,nreqv)), nprqp(MAX(1,nreqp)),           &
             nprqt(MAX(1,nreqt)),energ(6*iener), toutp(j) )

   READ(51) nprqa, nprqc, nprqd, nprql, nprqv, nprqp , nprqt
   READ(51) energ, sumat, toutd, (toutp(i),i=1,j)
   IF ( nreql > 0 ) CALL cmp_slist ( )
   READ (51) time, timed, cpui

   END SUBROUTINE rest_outp

   SUBROUTINE updlon_outp(oldlb)

   IMPLICIT NONE

   INTEGER(Kind=4), POINTER :: oldlb(:)

   INTEGER(Kind=4) :: i,j,lab,chnode
   LOGICAL :: flag


   ! accelerations
   j = 0
   flag = nreqa > 0
   DO i=1,nreqa
     lab = oldlb(nprqa(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j = j+1
       nprqa(j) = lab
     END IF
   END DO
   nreqa = j
   IF( flag .AND. nreqa == 0 )DEALLOCATE(nprqa)

   ! contact forces
   j = 0
   flag = nreqc > 0
   DO i=1,nreqc
     lab = oldlb(nprqc(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j = j+1
       nprqc(j) = lab
     END IF
   END DO
   nreqc = j
   IF( flag .AND. nreqc == 0 )DEALLOCATE(nprqc)

   ! displacements
   j = 0
   flag = nreqd > 0
   DO i=1,nreqd
     lab = oldlb(nprqd(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j = j+1
       nprqd(j) = lab
     END IF
   END DO
   nreqd = j
   IF( flag .AND. nreqd == 0 )DEALLOCATE(nprqd)

   ! internal forces
   j = 0
   flag = nreql > 0
   DO i=1,nreql
     lab = oldlb(nprql(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j = j+1
       nprql(j) = lab
     END IF
   END DO
   nreql = j
   IF( flag .AND. nreql == 0 )DEALLOCATE(nprql)

   ! velocities
   j = 0
   flag = nreqv > 0
   DO i=1,nreqv
     lab = oldlb(nprqv(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j = j+1
       nprqv(j) = lab
     END IF
   END DO
   nreqv = j
   IF( flag .AND. nreqv == 0 )DEALLOCATE(nprqv)

   ! temperatures
   j = 0
   flag = nreqt > 0
   DO i=1,nreqt
     lab = oldlb(nprqt(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j=j+1
       nprqt(j) = lab
     END IF
   END DO
   nreqt = j
   IF( flag .AND. nreqt == 0 )DEALLOCATE(nprqt)

   ! pressure
   j = 0
   flag = nreqp > 0
   DO i=1,nreqp
     lab = oldlb(nprqp(i))
     lab = chnode(lab)
     IF( lab > 0 )THEN
       j=j+1
       nprqp(j) = lab
     END IF
   END DO
   nreqp = j
   IF( flag .AND. nreqp == 0 )DEALLOCATE(nprqp)

   END SUBROUTINE updlon_outp


 END MODULE outp_db
