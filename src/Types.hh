#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;



// Enumeration of boundary conditions
typedef enum { NOSLIP = 0, SLIP = 1, OUTFLOW = 2, INFLOW = 3, PERIODIC = 4 } BcType;

// Enumeration of directions
typedef enum { NORTH, SOUTH, EAST, WEST } Direction;

// Enumeration of obstacle types
typedef enum { NO_OBSTACLE = 0, RECTANGLE = 1, CIRCLE = 2 } ObstacleType;


#endif //TYPES_HH
