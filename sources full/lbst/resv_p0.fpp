 LOGICAL :: large,   &! TRUE if log strains must be used
            ttti      ! TRUE if Trougth The Thickness Integration is necessary
 LOGICAL, ALLOCATABLE :: pflag(:) ! TRUE if plastic flow in the step

 LOGICAL :: newmt,   &! TRUE if material constant computation necessary
            natst,   &! TRUE for large strain analysis
            elast,   &! TRUE if material is strictly elastic
            shell,   &! TRUE if bending included
            plast,   &! TRUE if plasticity is to be considered
            shear,   &! TRUE if transverse shear stresses are considered in the constitutive relation
            visco,   &! TRUE if visco-elastic behavior
            coupled   ! TRUE if membrane-bending coupling for laminates

 INTEGER (kind=4) ielem, & !element number
                  i,j,k,n,l, &!different indexes
                  ilv,jlv,   &!indexes for layer varin
                  error    !error flag

 INTEGER (kind=4) isec,  & !associated material
               nl,nlayr, & !number of layers
                  mtype, & !associated material type
                  secty, & !section type
                  oldm,  & !old material label
                  mate,  & !material label
                  numpt, & !number of points in curve
                  osec,  & !associated material of previous element
                  nvar     !number of internal variables per layer

 REAL (kind=8) stres(3),  & !stresses (different measures)
               sigma(3),  & !stresses (different measures)
               stran(3),  & !C=U^2  also Log strains
               stral(3),  & !local strains
               strel(3),  & !local streses
               stine(10), & !t.t.t integrated stresses (forces and moments)
               r1,r2,     & !eigevenctor components in local system
               lb(3),     & !eigenvalues
               lc(3),     & !eigenvalues for rubbers
               thnew,     & !present thickness
               zk,        & !distance to mid surface
               aux,area1, & !auxiliar value
               !! t0,t1,j0,  & !thermical dilatation coeff
               s1,s2,efpst,strpl(3),dstpl(4),eql,stren
 REAL (kind=8) thick,     & !thickness (original)
               poiss,db,  & !poisson ratio, bending stiffness, elastic & plastic curvature
               c(4),gh,   & !Elastic constitutive matrix for plane stress
               dm(21),    & !Elastic integrated constitutive matrix
               alpha,     & !thermical dilatation coeff
               propi(13), & !Plastic Material properties
               chi(12),   & !Hill 48 coefficients
               chib(4),chid(6),  & !yield and potential coefficients
               deatht,    & !end time for plasticity
               minstr,    & !minimum strain to integrate trougth the thickness
               min_tr,    & !Minimum Allowable Thickness ratio
               max_tr       !Maximum Allowable Thickness ratio
 ! Gauss points throught the thickness
 REAL (kind=8), ALLOCATABLE :: thf(:),wei(:),shf(:)

 REAL (kind=8), ALLOCATABLE ::  thickl(:), zposl(:), rr(:,:)
 INTEGER (kind=4), ALLOCATABLE ::  lvari(:)
 REAL (kind=8), POINTER :: val(:,:)
 REAL (kind=8), ALLOCATABLE :: varin(:,:)                      !internal variables

 CHARACTER(len=1 ) letra  !for console input

 TYPE (section), POINTER :: sec         !pointer to a section data
 TYPE (mater), POINTER :: mat           !pointer to a material data

 TYPE (pmater), POINTER :: matsl(:)    !array of pointers to materials data

 INTERFACE
   INCLUDE 'gentbs.h'
   INCLUDE 'rubberps.h'
   INCLUDE 'core06ps.h'
   INCLUDE 'stre14.h'
   INCLUDE 'corr05.h'
   INCLUDE 'ud_shell_2.h'
 END INTERFACE

