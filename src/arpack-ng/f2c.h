/* f2c.h  --  Standard Fortran to C header file */
#pragma once

typedef struct { float r, i; } complex;
typedef struct { double r, i; } doublecomplex;
typedef int logical;

#define VOID void
#define TRUE_ (1)
#define FALSE_ (0)

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)
