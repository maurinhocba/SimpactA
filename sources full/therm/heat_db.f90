 MODULE heat_db
   USE param_db,ONLY: midn, mnam
   IMPLICIT NONE

   ! Derived type for a list containing nodal sources
   TYPE heat_nod
     INTEGER (kind=4) :: node            !node label
     REAL (kind=8) :: source(3)          !point source
     TYPE (heat_nod), POINTER :: next    !pointer to next node
   END TYPE heat_nod

   ! Derived type for a list containing edge heat (2D solid problem types only)
   TYPE edg_nod
     LOGICAL :: pos            ! .FALSE. = Bottom  .TRUE. = Top (surface) - not relevant in solids
     INTEGER (kind=4) :: nn2            ! 2 or 3, No of nodes defining the edge
     INTEGER (kind=4) :: lnod2(3)       ! node labels
     REAL (kind=8) :: flow2(3)          ! nodal press components
     TYPE (edg_nod), POINTER :: next    ! pointer to next segment
   END TYPE edg_nod

   ! Derived type for a list containing surface heat
   TYPE srf_nod
     LOGICAL :: pos            ! .FALSE. = Bottom  .TRUE. = Top (surface)
     INTEGER (kind=4) :: nn3            ! No of nodes defining the surface segm.
     INTEGER (kind=4) :: lnod3(6)       ! node labels
     REAL (kind=8) :: flow(6)           ! nodal flow
     TYPE (srf_nod), POINTER :: next    ! pointer to next segment
   END TYPE srf_nod

   ! Derived type for storing heat data (nodal,edges,surfaces)
   TYPE heat_set
     CHARACTER (len=mnam) :: lbl ! set label
     REAL (kind=8):: factor  ! scale factor
     INTEGER (kind=4) :: ntype !solid 2D type problem (segments/points)
     LOGICAL :: solid          !.TRUE. for solid problems (.FALSE. is shell)
     ! point heat
     INTEGER (kind=4) :: ipsour
     TYPE (heat_nod), POINTER :: headn, tailn
     ! data for edge heat
     INTEGER (kind=4) :: nedge
     TYPE (edg_nod), POINTER :: heade, taile
     ! data for surface heat
     INTEGER (kind=4) :: nsurf
     TYPE (srf_nod), POINTER :: heads, tails

     TYPE (heat_set), POINTER :: next
   END TYPE heat_set

   ! Derived type for a list containing convection and/or radiation surface or segments
   TYPE srf_con
     CHARACTER (len=mnam) :: sname ! set name
     LOGICAL :: conve, &       ! .TRUE. convection on this surface
                radia, &       ! .TRUE. radiation on this surface
                solid, &       ! .TRUE. for solid 2D and 3D problem
                carea          ! .TRUE. compute element area
     INTEGER (kind=4) :: nn, & ! No of nodes defining the surface segm.
                         nelem ! Number of surfaces defining the surface
     INTEGER (kind=4), POINTER :: lnods(:,:)  ! connectivities
     REAL (kind=8 ):: hf(6)     ! convection/radiation factors
                      !(1) hconc convection factor central surface (hconv for solids)
                      !(2) hcont convection factor top surface
                      !(3) hconb convection factor bottom surface
                      !(4) hradc radiation factor central surface  (hrads for solids)
                      !(5) hradt radiation factor top surface
                      !(6) hradb radiation factor bottom surface
     REAL (kind=8) :: tf(6)     ! convection/radiation temperatures
                      !(1) tambc ambient temperature central surface   (tambn for solids)
                      !(2) tambt ambient temperature top surface
                      !(3) tambb ambient temperature bottom surface
                      !(4) tradc radiation temperature central surface (trads for solids)
                      !(5) tradt radiation temperature top surface
                      !(6) tradb radiation temperature bottom surface
     INTEGER(kind=4) :: tcur(6) ! curve position for temperatures (pos 1 & 4 used in solids)
     CHARACTER(len=mnam) :: lbl(6) ! set label                    (pos 1 & 4 used in solids)
     REAL(kind=8), POINTER :: area(:)
     INTEGER(kind=4) :: afreq, &  ! steps to update area segments
                        mdof      ! number of DOFs per node
     TYPE (srf_con), POINTER :: next    ! pointer to next segment
   END TYPE srf_con

   ! variables for surfaces/segments/nodes with prescribed heat flow
   TYPE (heat_set), POINTER :: headhs, tailhs !first and last pointers

   ! variables for convection/radiation surfaces/segments
   INTEGER (kind=4) :: ncsur = 0 !number of convection/radiation surfaces
   TYPE (srf_con), POINTER :: headsc, tailsc !first and last pointers for surf or segm

 CONTAINS

   SUBROUTINE alloc_heats (heats)
     !allocates and initializes a heat set
     TYPE (heat_set), POINTER :: heats

     ALLOCATE (heats)
     CALL ini_heat (heats%headn, heats%tailn)
     CALL ini_edg (heats%heade, heats%taile)
     heats%nedge = 0
     CALL ini_srf (heats%heads, heats%tails)

   END SUBROUTINE alloc_heats

   SUBROUTINE ini_heat (head, tail)
     !initialize a list

     !Dummy arguments
     TYPE (heat_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_heat

   SUBROUTINE add_heat (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (heat_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a segment to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_heat

   ! routines for edge heats (2D solid problem types only)

   SUBROUTINE ini_edg (head, tail)
     !initialize a list

     !Dummy arguments
     TYPE (edg_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_edg

   SUBROUTINE add_edg (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (edg_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a segment to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_edg

   ! rotines for surface heats
   SUBROUTINE ini_srf (head, tail)
     !initialize a list

     !Dummy arguments
     TYPE (srf_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_srf

   SUBROUTINE add_srf (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (srf_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a segment to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_srf

   SUBROUTINE ins_srf (new, top)
     !insert a surface to the list
     TYPE (srf_nod), POINTER :: new, top

     IF (ASSOCIATED (top%next)) THEN
       new%next = top%next
       top%next => new
     ELSE
         NULLIFY(new%next)
         top%next => new
       ENDIF
   END SUBROUTINE ins_srf


   ! functions for heat sets

   SUBROUTINE ini_heats (head, tail)
     !initialize a list

     !Dummy arguments
     TYPE (heat_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_heats

   SUBROUTINE add_heats (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (heat_set), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a segment to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_heats

   SUBROUTINE srch_heats (head, anter, posic, lbl, found)
     !This subroutine searches for a heat set identified with a heat curve
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: lbl ! Load set label
     TYPE (heat_set), POINTER :: head, anter, posic

     found = .FALSE.
     NULLIFY (posic,anter)
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN
       posic => head
       DO
         IF(posic%lbl == lbl) THEN
           found = .TRUE.
           EXIT
         END IF
         IF (ASSOCIATED(posic%next) ) THEN
           anter => posic
           posic => posic%next
         ELSE
           EXIT
         END IF
       END DO
     ENDIF
     IF (.NOT.found) NULLIFY (posic,anter)
   END SUBROUTINE srch_heats

   SUBROUTINE del_heats (head, tail, anter, posic)
     !deletes a heat set pointed with posic
     TYPE (heat_set), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     ! if posic == tail
     IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
     CALL dalloc_heats (posic)
     NULLIFY (anter)
   END SUBROUTINE del_heats

   SUBROUTINE dalloc_heats (heats)
     ! deallocates a heat set
     TYPE (heat_set), POINTER :: heats
     TYPE (heat_nod), POINTER :: heat, heatux
     TYPE (edg_nod), POINTER :: edg, edgaux
     TYPE (srf_nod), POINTER :: srf, srfaux

     heat => heats%headn
     DO
       IF (.NOT.ASSOCIATED (heat) ) EXIT
       heatux => heat%next
       DEALLOCATE (heat)
       heat => heatux
     END DO

     edg => heats%heade
     DO
       IF (.NOT.ASSOCIATED (edg) ) EXIT
       edgaux => edg%next
       DEALLOCATE (edg)
       edg => edgaux
     END DO

     srf => heats%heads
     DO
       IF (.NOT.ASSOCIATED (srf) ) EXIT
       srfaux => srf%next
       DEALLOCATE (srf)
       srf => srfaux
     END DO

     DEALLOCATE (heats)
   END SUBROUTINE dalloc_heats

   !------------------------------------------------------------
   ! routines for convection/radiation surfaces

   SUBROUTINE ini_csr (head,tail)
     !initialize a surface with convection/radiation

     !Dummy arguments
     TYPE (srf_con), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_csr

   SUBROUTINE add_csr (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (srf_con), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a segment to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_csr

   SUBROUTINE ins_csr (new, top)
     !insert a surface to the list
     TYPE (srf_con), POINTER :: new, top

     IF (ASSOCIATED (top%next)) THEN
       new%next = top%next
       top%next => new
     ELSE
         NULLIFY(new%next)
         top%next => new
       ENDIF
   END SUBROUTINE ins_csr

  SUBROUTINE srch_csr(head, anter, posic, name, found)
  !This subroutine searches for a surface named "name"
  IMPLICIT NONE
    !--- Dummy arguments
    LOGICAL:: found
    CHARACTER(len=*):: name ! set name
    TYPE(srf_con),POINTER:: head, anter, posic

    found = .FALSE.
    NULLIFY(posic,anter)
    !Check if a list is empty
    IF (ASSOCIATED(head)) THEN
      posic => head
      DO
        IF (TRIM(posic%sname) == TRIM(name)) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next)) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    ENDIF
    IF (.NOT.found) NULLIFY(posic,anter)

  RETURN
  END SUBROUTINE srch_csr

  SUBROUTINE del_csr(head, tail, anter, posic)
  !deletes a surface pointed with posic
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(srf_con),POINTER:: head, tail, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    ! if posic == tail
    IF (.NOT.ASSOCIATED(posic%next)) tail=>anter
    CALL dalloc_csr(posic)
    NULLIFY(anter)

  RETURN
  END SUBROUTINE del_csr

  SUBROUTINE dalloc_csr(surface)
  ! deallocates a surface
  IMPLICIT NONE
    !--- Dummy variables
    TYPE(srf_con),POINTER:: surface

    DEALLOCATE(surface%lnods,surface%area)
    DEALLOCATE(surface)

  RETURN
  END SUBROUTINE dalloc_csr


 END MODULE heat_db
