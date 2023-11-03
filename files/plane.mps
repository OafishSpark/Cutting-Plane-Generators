* SCIP STATISTICS
*   Problem name     : /mnt/d/Projects/Cutting-Plane-Generators/files/plane.lp
*   Variables        : 2 (0 binary, 2 integer, 0 implicit integer, 0 continuous)
*   Constraints      : 2
NAME          /mnt/d/Projects/Cutting-Plane-Generators/files/plane.lp
OBJSENSE
  MIN
ROWS
 N  Obj 
 L  a1 
 L  a2 
COLUMNS
    INTSTART  'MARKER'                            'INTORG'                           
    x         Obj                             11  a1                             -11 
    x         a2                               1 
    y         Obj                          -13.1  a1                              13 
    y         a2                              10 
    INTEND    'MARKER'                            'INTEND'                           
RHS
    RHS       a1                             -10  a2                              90 
BOUNDS
 UP Bound     x                                9 
 PL Bound     y                                  
ENDATA