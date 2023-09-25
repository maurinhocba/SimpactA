 SUBROUTINE getnod(ncnod,ncnxx,nsegm,nnseg,lcnod,lcseg)

 !.... set list of nodes from connectivities for a surface
 !....     add to existing list
 !
 !.... input
 !....   ncnxx = Maximum number of nodes admited in the master surface considered
 !....   nsegm = number of segments of the master surface considered
 !....   nnseg = number of nodes per segment of the surface considered
 !....   lcseg(inseg,icseg) = global node number for the local element node
 !....                        [inseg] of the segment [icseg]
 !.... output
 !....   ncnod = number of nodes of the surface considered
 !....   lcnod(ncnod) = global node number for the local node [icnod] of the
 !....                  surface considered

 IMPLICIT NONE
 !     arguments
 INTEGER(kind=4),INTENT(IN) :: nsegm,nnseg,lcseg(:,:),ncnxx
 INTEGER(kind=4),INTENT(IN OUT) :: ncnod,lcnod(:)
 !     local variables
 INTEGER (kind=4)  jm,icseg,inseg

 DO icseg = 1,nsegm
   DO inseg = 1,nnseg
     jm = lcseg(inseg,icseg)
     IF(ALL(lcnod(1:ncnod) /= jm))THEN
       ncnod = ncnod+1
       IF( ncnod > ncnxx )THEN
         WRITE(55,"('GETNOD: Size:',i6)")ncnxx
         WRITE(55,"('Procesed:',i6,' of ',i6)")icseg,nsegm
         WRITE(55,"(10i8)")(lcnod(jm),jm=1,ncnxx)
         CALL RUNEND('GETNOD: Increase Auxiliar Memory !')
       END IF
       lcnod(ncnod) = jm
     END IF
   END DO
 END DO
 CALL sortnd (lcnod(1), ncnod)
 RETURN
 END SUBROUTINE getnod
