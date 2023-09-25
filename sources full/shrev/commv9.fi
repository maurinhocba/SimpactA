 SUBROUTINE commv9 (indic, nelem, nnode, ngaus, nstre, axesc, &
                    nreqs, narch, elsnam, elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic, nelem, nnode, ngaus, nstre, axesc, &
&                    nreqs, narch
 TYPE (ele09_set) :: elset

 IF(indic == 0) THEN
    elset%sname  = elsnam
    elset%nelem  = nelem
    elset%nnode  = nnode
    elset%ngaus  = ngaus
    elset%nstre  = nstre
    elset%axesc  = axesc
    elset%nreqs  = nreqs
    elset%narch  = narch
 ELSE
    elsnam = elset%sname
    nelem = elset%nelem
    nnode = elset%nnode
    ngaus = elset%ngaus
    nstre = elset%nstre
    axesc = elset%axesc
    nreqs = elset%nreqs
    narch = elset%narch
 END IF

 END SUBROUTINE commv9
