 SUBROUTINE commv7 (indic, nelem, nreqs, narch, stype, ngaus, elsnam, nnode, nstre, elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nelem, nreqs, narch, ngaus, stype, nnode, nstre
 TYPE (ele07_set) :: elset

 IF(indic == 0) THEN
    elset%sname  = elsnam
    elset%nelem  = nelem
    elset%nreqs  = nreqs
    elset%narch  = narch
    elset%stype  = stype
    elset%ngaus  = ngaus
    elset%nnode  = nnode
    elset%nstre  = nstre
 ELSE
    nelem  = elset%nelem
    elsnam = elset%sname
    nreqs   = elset%nreqs
    narch   = elset%narch
    stype   = elset%stype
    ngaus   = elset%ngaus
    nnode   = elset%nnode
    nstre   = elset%nstre
 END IF

 END SUBROUTINE commv7
