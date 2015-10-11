// -----------------------------------------------------------------
// Includes and definitions
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <qmp.h>
#include <qcdlib.h>

#if QDP_PrecisionInt == 2 || QDP_Precision == 'D'
#define QOP_PrecisionInt 2
#else
#define QOP_PrecisionInt 1
#endif

#define PI 3.1415926535897932384626433832795
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Functions
void getPoint(int j, int ndim, int *pnt, int argc, char *argv[]);

void mult_gamma_lr(QDP_DiracPropagator *dest,
                   QDP_DiracPropagator *src, int gl, int gr);

double print_mesons(QDP_DiracPropagator *prop, int t0, int dir, char *pat);
double print_corrs(QDP_DiracPropagator **prop, QDP_ColorMatrix **u,
                   int t0, int dir, int ls, char *pat);
// -----------------------------------------------------------------
