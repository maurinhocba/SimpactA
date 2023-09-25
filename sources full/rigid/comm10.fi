 SUBROUTINE comm10(indic,nelem,ntype,nnode,nmast,elsnam,elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nelem,ntype,nnode,nmast
 TYPE (ele10_set) :: elset

 IF(indic == 0) THEN
    elset%sname = elsnam
    elset%nelem = nelem
    elset%nnode = nnode
    elset%nmast = nmast
    elset%ntype = ntype
  ELSE
    elsnam = elset%sname
    nelem = elset%nelem
    nnode = elset%nnode
    nmast = elset%nmast
    ntype = elset%ntype
  END IF

 END SUBROUTINE comm10
