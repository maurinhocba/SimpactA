 MODULE param_db
 ! Module where global PARAMETERS of Program are defined
 IMPLICIT NONE

   !--- Constants for character lenght
   INTEGER(kind=4),PARAMETER:: midn = 6 ,&  !Maximum length for identifier
                               milb = 6 ,&  !Maximum length for internal labels
                               mchr = 6 ,&  !Maximum length for ???
                               miln = 132,&  !Maximum length for lines of data input
                               mlin = 131,&  !Maximum length for general lines
                               mttl = 99,&  !Maximum length for problem title
                               mset = 15,&  !Maximum length for sets comments
                               mvar = 15,&  !Maximum length for variable name and components
                               mich = 11    !Maximum length of digits for integers

   INTEGER(kind=4),PARAMETER:: mlen = 128,& !Maximum length for input parameters and file names
                               mprg = 8,  & !Maximum length for program label
                               mnam = 30    !Maximum length for data label names

   INTEGER(kind=4),PARAMETER:: mstl = 128, & !Maximum length for long character variables
                               msts = 30      !Maximum length for medium character variables

!--- Constants to limit usage
!MS$if acad > 0
   !--- for academic version
   INTEGER(kind=4),PARAMETER:: max_npoin_2d =     400, & !Maximum number of nodes por 2D problems
                               max_npoin_3d =    1000, & !Maximum number of nodes for 3D problems
                               max_memo     =   524000   !Maximum value for memo in beginp
!MS$else
   !--- full version
   INTEGER(kind=4),PARAMETER:: max_npoin_2d =  100000, & !Maximum number of nodes por 2D problems
                               max_npoin_3d = 1000000, & !Maximum number of nodes for 3D problems
                               max_memo     = 13000000   !Maximum value for memo in beginp
!MS$endif
 END MODULE param_db
