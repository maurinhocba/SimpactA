 SUBROUTINE comm05(indic,nnode,nelem,nreqs,narch,elsnam,elset,ngaus,quad)

 ! parameters of present set

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nnode,nelem,nreqs,narch,ngaus
 LOGICAL :: quad
 TYPE (ele05_set) :: elset

 IF(indic == 0) THEN
    elset%sname = elsnam         !set name
    elset%nnode = nnode          !number of nodes per element
    elset%nelem = nelem          !number of elements
    elset%nreqs = nreqs          !number of required elements for output
    elset%narch = narch          !file for output
    elset%ngaus = ngaus          !number of integration points
    elset%quad  = quad           !quadratic approach
  ELSE
    elsnam = elset%sname
    nnode  = elset%nnode
    nelem  = elset%nelem
    nreqs  = elset%nreqs
    narch  = elset%narch
    ngaus  = elset%ngaus
    quad   = elset%quad
  END IF

 RETURN
 END SUBROUTINE comm05
