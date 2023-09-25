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
END SUBROUTINE getnod
