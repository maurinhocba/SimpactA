 SUBROUTINE commv1(indic,nnode,nelem,nreqs,narch,elsnam,elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nnode,nelem,nreqs,narch
 TYPE (ele01_set) :: elset

 IF(indic == 0) THEN
    elset%nnode = nnode
    elset%sname = elsnam
    elset%nelem = nelem
    elset%nreqs = nreqs
    elset%narch = narch
  ELSE
    nnode = elset%nnode
    nelem = elset%nelem
    nreqs = elset%nreqs
    narch = elset%narch
    elsnam = elset%sname
  END IF

 END SUBROUTINE commv1
