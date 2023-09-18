 SUBROUTINE csdata(ncnod,nsegm,ncnxx,nsexx,nnseg,lcnod,lcseg,iwrit, &
                   coord,sname)

 !.... READ node numbers and segment connectivities

 USE surf_db
 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: nnseg,iwrit,ncnxx,nsexx
 INTEGER (kind=4), INTENT(IN OUT) :: ncnod,nsegm
 INTEGER (kind=4), INTENT(OUT) :: lcnod(:),lcseg(:,:)
! REAL(kind=8), INTENT(IN) :: coord(:,:)
 REAL(kind=8), POINTER :: coord(:,:)
 CHARACTER (len=*), INTENT(IN) :: sname
 END SUBROUTINE csdata
