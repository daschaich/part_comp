// -----------------------------------------------------------------
// Includes and definitions
#include "utils.h"
#include <complex.h>    // Including fftw3.h AFTER complex.h makes
#include <fftw3.h>      // fftw_complex the native complex double type
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print lengthy usage information if necessary
void print_usage(char *argv0) {
  char *fmt = " %-9s %s\n";

  printf("%s [options]\n", argv0);
  printf("options:\n");
  printf(fmt, "o pat", "input/output file pattern (required)");
  printf(fmt, "Q #", "maximum momentum");
  printf(fmt, "x #", "source point (four components)");
  printf("\n");
  QDP_abort(1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate Qhat_mu, Qhat^2 and Q^2
void getQ(int *pnt, int *size, double *Qhat, double *Qhatsq, double *Qsq) {
  double td;
  Qhatsq[0] = 0;
  Qsq[0] = 0;
  for (int mu = 0; mu < 4; mu++) {
    if (pnt[mu] <= size[mu] / 2) {
      td = 2 * PI * pnt[mu] / size[mu];
      Qhat[mu] = 2 * sin(td / 2);
    }
    else {  // Negative contributions
      td = 2 * PI * (pnt[mu] - size[mu]) / size[mu];
      Qhat[mu] = 2 * sin(td / 2);
    }

    *Qsq += td * td;
    *Qhatsq += Qhat[mu] * Qhat[mu];
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the row-major index for a 4d point
int getIndex(int pnt[4], int *size) {
  return pnt[3] + size[3] * (pnt[2] + size[2] * (pnt[1] + size[1] * pnt[0]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Match a given Qsq to the appropriate component
// in the given array, adding it if necessary
int map(double Qsq, double *moms, int *lenpt) {
  int i, len = lenpt[0];
  int test, key = floor(1e4 * Qsq);     // Avoid roundoff shenanigans

  // If Qsq already in array, return appropriate index
  if (len > 0) {
    for (i = 0; i < len; i++) {
      test = floor(1e4 * moms[i]);
      if (key == test)
        return i;
    }
  }

  // Otherwise add new Qsq to the array
  moms[len] = Qsq;
  lenpt[0] = len + 1;
  return len;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Perform Fourier transform and consolidate results
int decomp(QDP_Prop prop, int *latsize, QLA_Real Qmax,
           double *moms, int *src) {

  QLA_Real *data;
  int ndim, i, j, mu, nu, index = 0, N, label, len = 0, *count;
  int pnt[4] = {0, 0, 0, 0};
  double Qmu = 0, Qhat[4], Qhatsq = 0, Qsq = 0;
  complex double phase, new;

  fftw_complex *FT;                   // In-place transform
  fftw_plan plan;

  ndim = QDP_ndim();
  N = QDP_volume();

  // Arrays don't need to be so large, but to be safe
  count = (int *) malloc(N * sizeof(int));
  for (i = 0; i < N; i++) {
    longi[i] = 0;
    trans[i] = 0;
  }
  FT = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));

  // FFTW_MEASURE good for multiple transforms of the same size
  plan = fftw_plan_dft(ndim, latsize, FT, FT, 1, FFTW_MEASURE);

  // Perform Fourier transform for each (mu, nu)
  for (mu = 0; mu < ndim; mu++) {
    for (nu = 0; nu < ndim; nu++) {
      // Load (real) data, converting QDP_index to row-major index
      data = QDP_expose_R(pix[mu][nu]);
      for (i = 0; i < QDP_sites_on_node; i++) {
        QDP_get_coords(pnt, QDP_this_node, i);

        // Account for non-zero source point
        for (j = 0; j < ndim; j++)
          pnt[j] = (pnt[j] + latsize[j] - src[j]) % latsize[j];

        index = getIndex(pnt, latsize);
        FT[index] = data[i];
      }
      QDP_reset_R(pix[mu][nu]);

      fftw_execute(plan);

      // Calculate transverse--longitudinal decomposition.
      // Checked that both components are purely real
      for (i = 0; i < QDP_sites_on_node; i++) {
        count[i] = 0;             // i >= label, so won't overwrite data
        QDP_get_coords(pnt, QDP_this_node, i);

        // Account for non-zero source point
        for (j = 0; j < ndim; j++)
          pnt[j] = (pnt[j] + latsize[j] - src[j]) % latsize[j];

        getQ(pnt, latsize, Qhat, &Qhatsq, &Qsq);
        if (Qsq <= Qmax * Qmax && Qsq != 0) {
          index = getIndex(pnt, latsize);
          label = map(Qsq, moms, &len);
          count[label]++;               // For averaging

          // Apply phase to conserved--local correlator
          if (pnt[mu] <= latsize[mu] / 2)
            Qmu = PI * pnt[mu] / latsize[mu];
          else  // Negative contributions
            Qmu = PI * (pnt[mu] - latsize[mu]) / latsize[mu];

          phase = cos(Qmu) + I * sin(Qmu);
          new = phase * FT[index];

          longi[label] -= Qhat[mu] * new * Qhat[nu] / Qhatsq;
          if (mu == nu)   // Trace
            trans[label] += new / 3;
        }
      }
    }
  }

  // Normalize output
  for (i = 0; i < len; i++) {
    trans[i] += longi[i] / 3;
    longi[i] /= count[i];
    trans[i] /= count[i];
  }

  // Clean up
  free(count);
  fftw_destroy_plan(plan);
  fftw_free(FT);

  return len;
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Main method -- load propagator and perform Fourier transforms
int main(int argc, char *argv[]) {
  int i, j, mu, nu, ndim, *latsize, N, len = 0, src[] = {0,0,0,0};
  double *moms, Qmax = 10;
  char prec;
  QDP_Real *corr[4][4];
  char *outpat = NULL, fn[80];
  QDP_Reader *qr;
  QDP_String *md;
  FILE *out = NULL;

  if (argc < 2)
    print_usage(argv[0]);

  QDP_initialize(&argc, &argv);
  QDP_profcontrol(0);           // No profiling

  j = 0;    // Check for source point

  for(i = 1; i < argc; i++) {
    char c, *argp = &argv[i][1];
    c = argv[i][0];
    if(argv[i][1] == '\0') {
      if(i + 1 < argc) {
        argp = argv[i + 1];
        i++;
      }
      else argp = NULL;
    }
    switch (c) {
    case 'o': outpat = argp; break;
    case 'Q': Qmax = atof(argp); break;
    case 'x' :
      j = i;
      while (i + 1 < argc
            && ((argv[i + 1][0] == '-') || isdigit(argv[i + 1][0])))
        ++i;
      break;
    default: print_usage(argv[0]);
    }
  }

  if (outpat == NULL) {
    printf("Error: output pattern not specified!\n");
    print_usage(argv[0]);
  }

  sprintf(fn, "prop%s", outpat);
  QCD_latticePrecGetFromFile(&latsize, &ndim, &prec, fn);
  QCD_latticeSet(latsize, ndim);
  QDP_create_layout();

  N = QDP_volume();

  // Set up source point -- default zero
  getPoint(j, ndim, src, argc, argv);

  QCD_latticePrint();
  printf("Qmax = %g\n", Qmax);
  printf("source point =");
  for (i = 0; i < ndim; i++)
    printf(" %d", src[i]);
  printf("\noutpat = %s\n", outpat);

  // Arrays don't need to be so large, but to be safe
  moms = malloc(N * sizeof(double));

  // Load propagator and apply Fourier transform
  corr[mu][nu] = QDP_create_P();
  printf("Loading propagator file %s\n", fn);
  md = QDP_string_create();
  qr = QDP_open_read(md, fn);   // fn was set to axic above
  QDP_read_P(qr, md, corr[mu], 4);

  QDP_close_read(qr);
  len = project(prop, latsize, Qmax, moms, src)

  // Print out results
  // Keep imaginary parts to verify that output is purely real
  sprintf(fn, "project%s", outpat);
  printf("Writing output file %s\n", fn);
  out = fopen(fn, "w");
  for (i = 0; i < len; i++) {     // Won't be sorted
    fprintf(out, "%.4g\t"//TODO\n",
                  moms[i], //TODO);
  }
  fclose(out);

  // Clean up
  free(moms);
  for (mu = 0; mu < ndim; mu++) {
    for (nu = 0; nu < ndim; nu++)
      QDP_destroy_R(corr[mu][nu]);
  }

  QDP_finalize();
  return 0;
}
// -----------------------------------------------------------------
