// -----------------------------------------------------------------
// Includes and definitions
#include "utils.h"

#define realimag(i) ((i==0)||(i==7)||(i==11)||(i==13)||(i==14)||(i==15))
#define negate(i) (!((i==3)||(i==5)||(i==6)||(i==7)||(i==8)||(i==15)))
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function timeslice_func, with definition of QDP_Subset slices
static QDP_Subset *slices[4] = {NULL, NULL, NULL, NULL};
static int timeslice_func(int *coord, void *dum) {
  return coord[*((int *)dum)];
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up source point
void getPoint(int j, int ndim, int *pnt, int argc, char *argv[]) {
  int i;

  if (j == 0) {                 // No point given in input, use default
    for (i = 0; i < ndim; ++i)
      pnt[i] = 0;
  }
  else {
    if (isdigit(argv[j][0]))
      pnt[0] = atoi(argv[j]);
    else
      pnt[0] = atoi(argv[j] + 1);
    for (i = 1; i < ndim; ++i) {
      if ((++j < argc) && ((argv[j][0] == '-') || isdigit(argv[j][0])))
        pnt[i] = atoi(argv[j]);
      else
        pnt[i] = 0;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Multiply the input propagator by gamma matrices to the left and right
void mult_gamma_lr(QDP_DiracPropagator *dest,
                   QDP_DiracPropagator *src, int gl, int gr) {

  if ((gl == 0) && (gr == 0))
    QDP_P_eq_P(dest, src, QDP_all);

  else {
    QLA_DiracPropagator *s;
    QLA_DiracPropagator *d;
    s = QDP_expose_P(src);
    d = QDP_expose_P(dest);
    for (int i = 0; i < QDP_sites_on_node; i++) {
      QLA_DiracPropagator td;
      QLA_P_eq_gamma_times_P(&td, s + i, gl);
      QLA_P_eq_P_times_gamma(d + i, &td, gr);
    }
    QDP_reset_P(src);
    QDP_reset_P(dest);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out all 256 zero-momentum-projected meson correlators
double print_mesons(QDP_DiracPropagator *prop, int t0, int dir, char *pat) {
  int i, j, k, nt, t;
  double dtime;   // Time to print meson spectrum
  QDP_DiracPropagator *tp;
  char outfn[80];
  FILE *outfile = NULL;

  // Open file
  if (QDP_this_node == 0) {
    sprintf(outfn, "mesons%s", pat);
    outfile = fopen(outfn, "w");
  }
  dtime = -QDP_time();

  nt = QDP_coord_size(dir);
  if (!slices[dir])
    slices[dir] = QDP_create_subset(timeslice_func, &dir, sizeof(dir), nt);

  tp = QDP_create_P();

  QLA_Complex *z[16][16];   // All 256 possible meson correlators
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      z[i][j] = malloc(nt * sizeof(QLA_Complex));
      // gamma_5-hermiticity introduces factors of gamma_5
      mult_gamma_lr(tp, prop, 15 - i, 15 - j);
      QDP_c_eq_P_dot_P_multi(z[i][j], prop, tp, slices[dir], nt);
    }
  }

  for (k = 0; k < nt; k++) {
    t = (t0 + k) % nt;
    for (i = 0; i < 16; i++) {
      for (j = 0; j < 16; j++) {
        // Print only nonzero parts
        if (realimag(i) == realimag(j)) {
          fprintf0(outfile, "%i\t%i\t%i\t%i\t%.8g\n", t0, t, i, j,
                   QLA_real(z[i][j][t]));
        }
        else {
          fprintf0(outfile, "%i\t%i\t%i\t%i\t%.8g\n", t0, t, i, j,
                   QLA_imag(z[i][j][t]));
        }
      }
    }
  }

  // Clean up
  if (QDP_this_node == 0)
    fclose(outfile);

  QDP_destroy_P(tp);
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++)
      free(z[i][j]);
  }

  dtime += QDP_time();
  return dtime;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print zero-momentum-projected conserved--local temporal correlators
// to determine ZA (and optionally check ZV)
double print_corrs(QDP_DiracPropagator **prop, QDP_ColorMatrix **gauge,
                 int t0, int dir, int ls, char *pat) {

  int i, j, nt, t;
  double dtime;                 // Time to calculate and print correlators
  QDP_DiracPropagator *tpone, *tptwo, *tp;
  QLA_Complex *tc;
  QLA_Real *za, *zv;                // Zero-momentum-projected
  QDP_Shift shift;
  char projName[80];
  FILE *projFile = NULL;

  dtime = -QDP_time();
  nt = QDP_coord_size(dir);
  if (!slices[dir])
    slices[dir] = QDP_create_subset(timeslice_func, &dir, sizeof(dir), nt);

  tpone = QDP_create_P();
  tptwo = QDP_create_P();
  tp = QDP_create_P();

  tc = (QLA_Complex *) malloc(nt * sizeof(QLA_Complex));
  zv = (QLA_Real *) malloc(nt * sizeof(QLA_Real));
  za = (QLA_Real *) malloc(nt * sizeof(QLA_Real));
  for (i = 0; i < nt; i++) {
    za[i] = 0.0;
    zv[i] = 0.0;
  }

  // Calculate correlator of temporal component of conserved
  // axial (vector) current with local pseudoscalar (scalar)
  // Since realimag(7) = realimag(15) but realimag(8) != realimag(0),
  // imaginary part of za and real part of zv vanish (checked)
  shift = QDP_neighbor[dir];
  for (i = 0; i < ls; i++) {
    if (i == ls / 2) {    // Negate left side of axial correlator
      for (j = 0; j < nt; j++)
        za[j] *= -1;
    }

    // P_-4 term
    QDP_P_eq_sP(tpone, prop[i], shift, QDP_forward, QDP_all);
    QDP_P_eq_M_times_P(tptwo, gauge[3], tpone, QDP_all);
    QDP_P_eq_gamma_times_P(tpone, tptwo, 8, QDP_all);
    QDP_P_meq_P(tptwo, tpone, QDP_all);               // P_-4
    QDP_P_eq_gamma_times_P(tpone, tptwo, 15, QDP_all);
    QDP_c_eq_P_dot_P_multi(tc, prop[ls - 1 - i], tpone, slices[dir], nt);
    for (j = 0; j < nt; j++)
      QLA_r_peq_Re_c(za[j], tc[j]);
    QDP_P_eq_P_times_gamma(tp, tpone, 15, QDP_all);
    QDP_c_eq_P_dot_P_multi(tc, prop[ls - 1 - i], tp, slices[dir], nt);
    // Extra negative sign from negate(8) != negate(0)
    for (j = 0; j < nt; j++)
      QLA_r_meq_Im_c(zv[j], tc[j]);

    // P_4 term
    QDP_P_eq_Ma_times_P(tpone, gauge[3], prop[i], QDP_all);
    QDP_P_eq_gamma_times_P(tptwo, tpone, 8, QDP_all);
    QDP_P_peq_P(tpone, tptwo, QDP_all);               // P_4
    QDP_P_eq_gamma_times_P(tptwo, tpone, 15, QDP_all);
    QDP_P_eq_sP(tpone, prop[ls - 1 - i], shift, QDP_forward, QDP_all);
    QDP_c_eq_P_dot_P_multi(tc, tpone, tptwo, slices[dir], nt);
    for (j = 0; j < nt; j++)
      QLA_r_meq_Re_c(za[j], tc[j]);
    QDP_P_eq_P_times_gamma(tp, tptwo, 15, QDP_all);
    QDP_c_eq_P_dot_P_multi(tc, tpone, tp, slices[dir], nt);
    for (j = 0; j < nt; j++)
      QLA_r_peq_Im_c(zv[j], tc[j]);
  }

  // Open file for printing
  if (QDP_this_node == 0) {
    sprintf(projName, "mesonsc%s", pat);
    projFile = fopen(projName, "w");
  }
  // Normalize spin projections and print
  for (i = 0; i < nt; i++) {
    t = (t0 + i) % nt;
    fprintf0(projFile, "%i\t%i\t7\t15\t%.8g\n", t0, t, za[t] / 2);
    fprintf0(projFile, "%i\t%i\t8\t0\t%.8g\n", t0, t, zv[t] / 2);
  }

  // Clean up and close down
  if (QDP_this_node == 0)
    fclose(projFile);
  QDP_destroy_P(tpone);
  QDP_destroy_P(tptwo);
  QDP_destroy_P(tp);
  free(zv);
  free(za);

  dtime += QDP_time();
  return dtime;
}
// -----------------------------------------------------------------
