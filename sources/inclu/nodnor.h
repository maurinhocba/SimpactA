 SUBROUTINE nodnor( nelem,nnode,lnods,coord,eule0,euler)
 !******************************************************************************
 !
 !****this routine compute the position angles of the nodal coordinate systems
 !
 !     input:   nelem:         number of elements
 !              nnode:         number of nodes per element
 !              coord:         nodal coordinates
 !              lnods:         element conectivities
 !              eule0:  /= 0   angles of the nodes. not modified
 !                      == 0   nodes where to compute eule0 angles
 !
 !     output:  eule0:  angles of the nodes.
 !
 !******************************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER :: nd = 3
 INTEGER (kind=4), INTENT(IN) :: nelem,nnode, lnods(:,:)
 REAL (kind=8), INTENT(IN) :: coord(:,:)
 REAL (kind=8), INTENT(IN OUT) :: eule0(:,:),euler(:,:)
 END SUBROUTINE nodnor
