 SUBROUTINE rdsldn (iwrit,actio)

 !This routine reads Sliding nodes

 USE param_db,ONLY: mnam
 USE c_input
 USE nsld_db
 IMPLICIT NONE
 CHARACTER(len=*),INTENT(INOUT):: actio
 INTEGER :: iwrit

 ! Local
 TYPE (nsld_nod), POINTER :: nsld,head,tail
 !TYPE (nset), POINTER :: ns
 INTEGER :: slave,master,i,j,n
 LOGICAL :: dummy,autom,nosld,norot

 IF (iwrit == 1) WRITE(lures,"(//,3x,'Slave/Master sliding Data',/)",ERR=9)

 !Initialize empty list
 CALL ini_sld(head,tail)
 IF (TRIM(actio) == 'NEW' .OR. .NOT.exists('ADD   ')) THEN
   nnsld = 0                       !initializes list
 ELSE IF (exists('ADD   ')) THEN
   IF (iwrit == 1) WRITE(lures,                                 &
      "(//,3x,'Kept Slave/Master sliding nodes Data from the previous strategy',/)",ERR=9)
   DO n=1,nnsld                       !for previous data
     ALLOCATE (nsld)                 !get memory
     nsld%slave = nsldat(1,n)        !slave node
     nsld%master = nsldat(2,n)       !master node
     nsld%autom  = nsldfg(1,n)       !AUTOMatic flag
     nsld%nosld  = nsldfg(2,n)       !NO SLiDing flag
     nsld%norot  = nsldfg(3,n)       !NO ROTation flag
     CALL add_sld (nsld, head, tail)  !add to list
   END DO
 END IF

 IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'
 IF (ASSOCIATED (nsldat) ) DEALLOCATE (nsldat, nsldfg) !release memory

 IF (iwrit == 1) WRITE(lures,"(/'   Slave Node     Master Node  AUTOM  NO-SLD  NO-ROT')",ERR=9)

 !Loop to read dat and add them to the list
 DO
   CALL listen('READP0')             !read a line
   IF (exists('ENDSLI')) EXIT        !end of data found => exit
   dummy = exists('SLAVE ',i) .AND. exists('MASTER',j)
   IF( dummy ) THEN
     slave = INT(param(i))
     master= INT(param(j))
     autom = exists('AUTOM')
     nosld = exists('NOSLD')
     norot = exists('NOROT')
   ELSE
     slave = INT(param(1))
     master= INT(param(2))
     autom = INT(param(3)) == 1
     nosld = INT(param(4)) == 1
     norot = INT(param(5)) == 1
     !master = getint('MASTER',0,'!-Number of a Master node .........') !get master node
     !slave  = getint('SLAVE ',0,'!Number of a Slave node ...........')  !get slave node
   END IF
   IF( nosld .AND. norot ) CALL runend('Both NO-SLD & NO-ROT are not possible, use Dependant nodes instead')
   IF( slave /= 0 .AND. master /= 0 )THEN
     ALLOCATE (nsld)                    !get memory for a new constraint
     nsld%slave  = slave                !assign
     nsld%master = master               !assign
     nsld%autom  = autom                !assign
     nsld%nosld  = nosld                !assign
     nsld%norot  = norot                !assign
     IF (iwrit == 1) WRITE(lures, '(2i10,3(L7,1x))',ERR=9) slave,master, autom, nosld, norot
     nnsld=nnsld+1                      !increase number of dependencies
     CALL add_sld( nsld, head, tail )   !add to list
   ELSE          !set not found => ERROR
     WRITE (lures, "(' Invalid data for SLIDING node/s')",ERR=9)
     CALL runend('Kinematic constraints: invalid data ')
   END IF

 END DO

 !Store data in an array
 ALLOCATE (nsldat(2,nnsld),nsldfg(3,nnsld))       !get memory for array
 CALL store_sld(head)                !transfer from list to array
                                     !and release memory
 CALL ini_sld(head,tail)             !nullifies pointer

 RETURN
 9 CALL runen2(' error while writing data to the disk')
 END SUBROUTINE rdsldn
