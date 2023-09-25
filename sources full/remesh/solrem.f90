SUBROUTINE solrem(refset)
!-----------------------------------------------------------------------
! remeshing routine for solid elements
!-----------------------------------------------------------------------
USE damp_db, ONLY:  updlon_damp
USE c_input, ONLY: openfi, ludat
USE npo_db,ONLY: coora, velnp, oldlb
USE meshmo_db,ONLY: itransf, numpn, coorn, v_new, strnd, p_new,   &
                    sg_old, sg_new, l_new, ne_new, elemp, nodset, headp,  &
                    tailp, delete_pline, &
                    nelbnd, cdbnd, lnbnd, nlbnd, fside, nodlb   !remeshing database
USE ctrl_db, ONLY: ntype,ndime,ndofn,nload,echo_chnode
USE kinc_db, ONLY : nvelr,naris
USE surf_db,ONLY: surfa, new_srf, dalloc_srf
USE param_db,ONLY: mlen
USE name_db,ONLY: output
USE SPR_Constants, ONLY: nvarg
USE cms_db, ONLY : nconm,updlon_cm
!USE outp_db, ONLY : updlon_outp
USE loa_db, ONLY : updlon_lo
!#if UNIX
!     USE ifport
!#else
   USE dflib
   USE dfport
!#endif
IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=*),INTENT(IN):: refset   !label of the element set to modify
  !--- Local variables
  LOGICAL          :: found, axisym
  CHARACTER (len=mlen) :: command
  CHARACTER(len=mlen), PARAMETER :: dirgid ='C:\Utils\GiD\Compack1.0.5\'
  CHARACTER(len=mlen):: giddir
  CHARACTER(len=mlen):: gidv = 'SIMPACT'   !len=32
  INTEGER (kind=4) :: lengn = 27
  INTEGER (kind=4) :: istop=0, ng, nelm, nnod, ludat_keep
  INTEGER (kind=4), POINTER :: ifx_b(:,:)   !  boundary fixities
  INTEGER(kind=4),POINTER:: l_old(:,:)      ! (nnode,nelem) old elements
  REAL    (kind=8) :: fact

  REAL(kind=8):: scalef

  INTERFACE
    INCLUDE 'contac.h'
    INCLUDE 'elemnt.h'
    INCLUDE 'get_ifxb.h'
    INCLUDE 'plines.h'
    INCLUDE 'pregid.h'
    INCLUDE 'qloop.h'
    INCLUDE 'rmctrln.h'
    INCLUDE 'rnmesh.h'
    INCLUDE 'transf.h'
    INCLUDE 'transf1.h'
    INCLUDE 'transf_g.h'
  END INTERFACE

  found = .FALSE.
  axisym = ntype == 3

  !Extract a set of nodes used in the definition of refset
  CALL elemnt('NODSET',name=refset,istop=istop, flag2=found)  !==>nodset(numpo)
  !  Get conectivities of the set
  CALL elemnt('SLNODS',name=refset,lnod=l_old,flag2=found) ! ==>  l_old(nnod,nelm)
  nelm = SIZE(l_old,2)   !Get number of elements in the set
  nnod = SIZE(l_old,1)   !Get number of nodes by element
  ! get external surface/line definition
  CALL new_srf(surfa)         !Get memory for surface
  CALL elemnt('BOUNDA', name=refset, flag2=found)  !==>surfa
  !--- Generation of a list of nodes in the boundary
  CALL rmctrln(surfa%nelem,surfa%head,nelbnd,cdbnd,lnbnd,nlbnd,fside)
                                    !==> lnbnd(nelbnd)       Conectivities of segments
                                    !==> nlbnd(nline)        Number of lines in boundary
                                    !==> cdbnd(ndime,nelbnd) Coordinates of nodes
  CALL dalloc_srf(surfa)   !Release 'surfa' from memory
  ! Eliminates loops
  CALL qloop(nelbnd,cdbnd,nlbnd)
  !get boundary conditions on the boundary line ==> ifx_d
  CALL get_ifxb(nelbnd,lnbnd,ifx_b)
  ! analyze contour fixity codes and vertices to find polylines
  CALL plines(cdbnd,nlbnd,fside,ifx_b)
  ! prepare GiD batch file
  fact = scalef     !keep scale factor
  scalef = 1d4      !modify factor previous to writing background mesh
  CALL pregid(scalef,nnod,nelbnd,cdbnd,nlbnd,ifx_b)  !write background mesh

  CALL delete_pline(headp,tailp)  !Release memory of the boundary arrays
  ! generation of a new mesh
  lengn = GETENVQQ(gidv, giddir)     !get enviromental variable ==> giddir
  IF( lengn <= 0 )THEN
    giddir = dirgid
    lengn = LEN_TRIM (giddir)
  END IF
!#if UNIX
!  command='".'//giddir(1:lengn)//'gid.exe" -n -b '//TRIM(output)//'.bgi'
!  istop = SYSTEM(TRIM(command))
!#else
  command='"'//giddir(1:lengn)//'gid.exe" -n2 -b '//TRIM(output)//'.bgi' ! -n... not work in GiD 9 and later for Windows
  istop = SYSTEM(TRIM(command))
!#endif

!#if UNIX
!  command = 'rm -rf .//'//TRIM(output)//'.gid'
!  istop = SYSTEM(TRIM(command))   ! delete GiD files (& directory)
!#else
  command = 'rmdir /S /Q .\'//TRIM(output)//'.gid'
  istop = SYSTEM(TRIM(command))   ! delete GiD files (& directory)
!#endif

  CALL openfi(54)         !open mesh file to read
  ludat_keep = ludat      !remember standard value of LUDAT
  ludat = 54              !changed temporarily to allow the use of listen with this unit

  !read new coordinates and connectivities
  CALL rnmesh(scalef,f_axsym=axisym)
  scalef = fact !restore value once coordinates have been read

  ! transfer all the nodal variables from the old mesh to the new one
  IF (nnod ==3)ng=1          !number of Gauss points (1 for triangles)
  IF (nnod ==4)ng=4          !number of Gauss points (4 for quads)
  IF (itransf == 2) THEN      ! transfer from nodal smoothed vars
    CALL transf(numpn, coorn, coora, l_old, velnp, v_new ) ! transfer nodal vars
    !Smoothing of Gaussian variables in the remeshed set of elements
    CALL elemnt('SMOOTH', name=refset, istop=istop, flag2=found) !==>strnd(nstre,npoin)
    nvarg = SIZE(strnd,1)     !number of Gauss variables to transfer
    ALLOCATE(sg_new(nvarg,ng*ne_new))  !get memory for new Gauss points
    CALL transf1(ne_new,nnod,ng,coorn,coora,l_old,l_new,sg_new,strnd) ! transfer gauss vars
  ELSE IF (itransf == 1) THEN  !transfer all GP variables using
                               !superconvergent patch recovery (SCPR) algorithm
    CALL transf(numpn, coorn, coora, l_old, velnp, v_new)     ! transfer nodal vars
    CALL elemnt('INTVAR', name=refset, istop=istop, flag2=found) !==>sg_old(nvarg,nelem)
    nvarg = SIZE(sg_old,1)     !number of Gauss variables to transfer
    ALLOCATE(sg_new(nvarg,ng*ne_new))  !get memory for new Gauss points
    CALL transf_g(numpn, coorn,  coora, l_old, sg_old, l_new, sg_new)  ! transfer GP vars
  END IF
  ! exchange nodes - delete nodes of the old mesh and add new nodes ==> modify nodal database
  CALL exnods(refset,.TRUE.)   !Updated lagrangian
  ! reallocates the element set after remeshing, compute Gaussian variables
  ! exchange mesh - deallocate old mesh and put new mesh,
  CALL elemnt('REALLO', name=refset, flag2=found)
  ! update boundary conditions, delete those associated with old mesh & add new
  CALL newfix( )
  CLOSE(54,STATUS='DELETE')      !close and delete new mesh data file
  ludat = ludat_keep    !restore read file to original

  !--- Release memory

  DEALLOCATE(coorn, l_old, ifx_b)
  DEALLOCATE(lnbnd,cdbnd,nlbnd,fside)
  IF (ASSOCIATED(nodlb))  DEALLOCATE(nodlb)
  IF (ASSOCIATED(l_new))  DEALLOCATE(l_new)
  IF (ASSOCIATED(nodset)) DEALLOCATE(nodset)
  IF (ASSOCIATED(elemp))  DEALLOCATE(elemp)
  IF (ASSOCIATED(strnd))  DEALLOCATE(strnd)
  IF (ASSOCIATED(sg_new)) DEALLOCATE(sg_new)
  IF (ASSOCIATED(sg_old)) DEALLOCATE(sg_old)
  IF (ASSOCIATED(v_new))  DEALLOCATE(v_new)
  IF (ASSOCIATED(p_new))  DEALLOCATE(p_new)

  ! updates data bases
  echo_chnode = .FALSE.
  CALL elemnt ('UPDLON')   ! elements/segments connectivities
  CALL updlon_kc (naris,nvelr,ndime)   ! KINEMATIC conditions
  CALL updlon_cm (ndofn)               ! CONCENTRATED MASSES
  !CALL updlon_outp (oldlb)            ! output nodes (unnecesary)
  CALL updlon_damp ()                  ! DAMPING nodes and sets
  CALL updlon_lo (nload,oldlb)         ! LOADS
  CALL contac (task='UPDLON',iwrit=0)             ! CONTACT nodes (CONT 2 & 3 & 4 only)
  echo_chnode = .TRUE.

  RETURN
END SUBROUTINE solrem
