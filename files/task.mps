* SCIP STATISTICS
*   Problem name     : /mnt/d/Projects/Cutting-Plane-Generators/files/task.lp
*   Variables        : 2 (0 binary, 2 integer, 0 implicit integer, 0 continuous)
*   Constraints      : 2
NAME          /mnt/d/Projects/Cutting-Plane-Generators/files/task.lp
OBJSENSE
  MIN
ROWS
 N  Obj 
 L  a1 
 L  a2 
COLUMNS
    INTSTART  'MARKER'                            'INTORG'                           
    x         Obj                             -1  a1                               4 
    x         a2                               3 
    y         Obj                           -1.1  a1                               3 
    y         a2                               4 
    INTEND    'MARKER'                            'INTEND'                           
RHS
    RHS       a1                              10  a2                              10 
BOUNDS
 PL Bound     x                                  
 PL Bound     y                                  
ENDATA