 SUBROUTINE rdheat (iwrit,ndoft,nheat,actio)

 !This routine reads heat sets (nodal heats, edge and surface heats)

 USE ctrl_db, ONLY: ntype,ndime 
 USE c_input
 USE heat_db
 USE esets_db, ONLY : delete_name, add_name
 IMPLICIT NONE
 CHARACTER(len=6 ) :: actio
 INTEGER (kind=4) :: iwrit,ndoft,nheat

 END SUBROUTINE rdheat
