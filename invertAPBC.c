// -----------------------------------------------------------------
// Includes and definitions
#include "utils.h"

#define PRERES 5e-6   // Float inner residual
#define MAX_IT 10000
#define RESTART 500
#define MAX_RESTARTS 20
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function print_usage just prints lengthy usage info if necessary
void print_usage(char *argv0) {
  char *fmt = " %-9s %s\n";

  printf0("\n%s [options]\n", argv0);
  printf0("options:\n");
  printf0(fmt, "l file", "lattice file (required)");
  printf0(fmt, "L #", "Ls");
  printf0(fmt, "m #", "mass");
  printf0(fmt, "M #", "M5, in the range [0,2]");
  printf0(fmt, "n #", "number of smearings");
  printf0(fmt, "o pat", "output pattern (required)");
  printf0(fmt, "p #", "profile if nonzero");
  printf0(fmt, "r #", "residual");
  printf0(fmt, "s #", "smearing factor");
  printf0(fmt, "x #", "source point (four components)");
  printf0("\n");
  QDP_abort(1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Main function
int main(int argc, char *argv[]) {
  int i, j, mu, ndim, ls = 16, nsmear = 0, ic, is, prof = 0;
  int *pnt, *latsize;
  double res = 1e-8, invtime = 0.0, dtime = 0.0;
  double BCs[4] = {-1.0, -1.0, -1.0, -1.0};
  char *latfile = NULL, *outpat = NULL, outpat_mid[80], outprop[80];
  char latprec;
  QLA_Real m = 0.0, M5 = 1.8, rsmear = 0.0;
  QLA_Real plaq, splaq, tplaq, *plaqs;
  QDP_ColorMatrix **gauge;
  QDP_DiracFermion **out, **in, *td;
  QDP_DiracPropagator *wall, *mid, **prop;
  QCD_DW *dw;
  QCD_perf_t perf;

  QDP_initialize(&argc, &argv);

  if (argc < 3)
    print_usage(argv[0]);

  j = 0;    // Explicitly check to see if source point has been given
  for (i = 1; i < argc; i++) {
    char c, *argp = &argv[i][1];
    c = argv[i][0];
    if (argv[i][1] == '\0') {
      if (i + 1 < argc) {
        argp = argv[i + 1];
        i++;
      }
      else
        argp = NULL;
    }
    switch(c) {
    case 'l' : latfile = argp; break;
    case 'L' : ls = atoi(argp); break;
    case 'm' : m = atof(argp); break;
    case 'M' : M5 = atof(argp); break;
    case 'n' : nsmear = atoi(argp); break;
    case 'o' : outpat = argp; break;
    case 'p' : prof = atoi(argp); break;
    case 'r' : res = atof(argp); break;
    case 's' : rsmear = atof(argp); break;
    case 'x' :
      j = i;
      while (i + 1 < argc
            && ((argv[i + 1][0] == '-') || isdigit(argv[i + 1][0])))
        ++i;
      break;
    default : print_usage(argv[0]);
    }
  }

  QDP_profcontrol(prof);                      // Profiling

  if (latfile == NULL) {
    printf0("Error: lattice file not specified!\n\n");
    print_usage(argv[0]);
  }
  if (outpat == NULL) {
    printf0("Error: output pattern not specified!\n\n");
    print_usage(argv[0]);
  }

  sprintf(outpat_mid, "%s%s", "m", outpat);

  QCD_latticePrecGetFromFile(&latsize, &ndim, &latprec, latfile);
  QCD_latticeSet(latsize, ndim);
  QDP_create_layout();

  pnt = malloc(ndim * sizeof(pnt));
  getPoint(j, ndim, pnt, argc, argv);         // Get source point

  if (QDP_this_node == 0) {
    QCD_latticePrint();
    printf("Ls = %d\n", ls);
    printf("M5 = %g\n", M5);
    printf("m = %g\n", m);
    printf("res = %g\n", res);
    printf("source point =");
    for (i = 0; i < ndim; i++)
      printf(" %d", pnt[i]);
    printf("\n");
    printf("smear factor = %g\n", rsmear);
    printf("smear iterations = %i\n", nsmear);
    printf("profiling = %d\n", prof);
    printf("outpat = %s\n", outpat);
  }

  // Load gauge field
  gauge = malloc(ndim * sizeof(QDP_ColorMatrix *));
  for (mu = 0; mu < ndim; mu++)
    gauge[mu] = QDP_create_M();
  QCD_gaugePrecLoad(gauge, latfile, latprec);

  // Calculate and print plaquette
  plaqs = malloc((ndim - 2) * (ndim - 1) * sizeof(QLA_Real));
  QCD_gaugePlaq(gauge, &plaq, &splaq, &tplaq, plaqs);
  printf0("Plaquette = %.8g\n", plaq);
  free(plaqs);    // Done with this

  // Init inverter
  dw = QCD_dwNew();                     // Calls dwInit
  QOP_verbose(QOP_VERB_OFF);            // Debugging
  dw->ls = ls;
  dw->res_arg.rsqmin = res * res;       // Squared residual!
  dw->float_inner = 1;                  // Turns on float_inner
  dw->preres = PRERES;
  dw->inv_arg.max_iter = MAX_IT;
  dw->inv_arg.restart = RESTART;
  dw->inv_arg.max_restarts = MAX_RESTARTS;
  QCD_gaugeSetBC(gauge, BCs, ndim);     // Should give APBC in all four dirs
  QCD_dwLoadLinks(dw, gauge, 1);        // APBC already accounted for above

  td   = QDP_create_D();
  wall = QDP_create_P();
  mid  = QDP_create_P();

  in  = malloc(ls * sizeof(QDP_DiracFermion *));
  out = malloc(ls * sizeof(QDP_DiracFermion *));
  prop = malloc(ls * sizeof(QDP_DiracPropagator *));
  for (i = 0; i < ls; i++) {
    in[i]  = QDP_create_D();
    out[i] = QDP_create_D();
    prop[i] = QDP_create_P();
  }

  // Build propagators by looping through colors and spins
  for (ic = 0; ic < 3; ic++) {
    for (is = 0; is < 4; is++) {
      // Put source on left or right wall depending on chirality
      QDP_D_eq_zero(td, QDP_all);                       // Delete old point
      QCD_setPoint_D(td, pnt, ndim, ic, is, 1.0, 0.0);  // Add new point

      if (nsmear)
        QCD_fermionGaussianSmear(td, gauge, rsmear, nsmear, 1.0);

      for (i = 0; i < ls; i++) {
        QDP_D_eq_zero(out[i], QDP_all);
        if (i == 0)
          QDP_D_eq_spproj_D(in[i], td, 4, 1, QDP_all);
        else if (i == ls - 1)
          QDP_D_eq_spproj_D(in[i], td, 4, -1, QDP_all);
        else
          QDP_D_eq_zero(in[i], QDP_all);
      }

      printf0("color = %d, spin = %d\n", ic, is);

      // Pass full 5d fermions; residual in data structure
      perf = (QCD_perf_t){0, 0, 0};
//      QCD_dwInvertInner(&perf, dw, M5, m, out, in);
      QCD_dwInvert5d(&perf, dw, M5, m, out, in);
      invtime += perf.secs;

      if (nsmear) {
        for (i = 0; i < ls; i++)
          QCD_fermionGaussianSmear(out[i], gauge, rsmear, nsmear, 1.0);
      }

      // Array of propagators along fifth dimension
      for (i = 0; i < ls; i++)
        QDP_P_eq_diracvec_D(prop[i], out[i], ic, is, QDP_all);

      // Build four dimensional domain wall propagator
      // and midpoint propagator from appropriate chiral components
      QDP_D_eq_spproj_D(td, out[0], 4, -1, QDP_all);
      QDP_D_peq_spproj_D(td, out[ls - 1], 4, 1, QDP_all);
      QDP_P_eq_diracvec_D(wall, td, ic, is, QDP_all);

      QDP_D_eq_spproj_D(td, out[ls / 2], 4, -1, QDP_all);
      QDP_D_peq_spproj_D(td, out[ls / 2 - 1], 4, 1, QDP_all);
      QDP_P_eq_diracvec_D(mid, td, ic, is, QDP_all);
    }
  }
  printf0("Inversion time = %.3f seconds\n", invtime);

  // Done with these
  QCD_dwFree(dw);
  QDP_destroy_D(td);
  for (i = 0; i < ls; i++) {
    QDP_destroy_D(in[i]);
    QDP_destroy_D(out[i]);
  }
  free(in);
  free(out);

  // Print wall and midpoint zero-momentum-projected local--local correlators
  dtime = print_mesons(wall, pnt[3], ndim - 1, outpat);
  printf0("Meson time = %g secs\n", dtime);
  dtime = print_mesons(mid, pnt[3], ndim - 1, outpat_mid);
  printf0("Midpoint time = %g secs\n", dtime);
  QDP_destroy_P(mid);    // Done with this

  // Save 4d propagator
  sprintf(outprop, "prop%s", outpat);
  QCD_propSave(wall, outprop);
  QDP_destroy_P(wall);    // Done with this

  // Print full and zero-momentum-projected conserved--local correlators
  dtime = print_corrs(prop, gauge, pnt[3], ndim - 1, ls, outpat);
  printf0("ZA time = %g secs\n", dtime);

  // Clean up and close down
  for (mu = 0; mu < ndim; mu++)
    QDP_destroy_M(gauge[mu]);
  free(gauge);
  for (i = 0; i < ls; i++)
    QDP_destroy_P(prop[i]);
  free(prop);

  QDP_finalize();
  return 0;
}
// -----------------------------------------------------------------
