#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// CFL stability factor: dt = min(dx2,dy2) / STABILITY_FACTOR
#define STABILITY_FACTOR 4.e0

static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline int    imin(int a, int b)       { return a < b ? a : b; }
static inline int    imax_i(int a, int b)     { return a > b ? a : b; }

// Linear index into phi, the full grid including boundaries: (is:ie, ks:ke).
static inline int idx(int i, int k, int is, int ks, int kouter) {
  return (i - is) * kouter + (k - ks);
}

// Linear index into phin, the inner grid excluding boundaries:
// (is+halo:ie-halo, ks+halo:ke-halo).
static inline int idxn(int i, int k, int is, int ks, int kinner) {
  return (i - is - 1) * kinner + (k - ks - 1);
}

// Index convention:
//   i = 1st array coordinate = y
//   k = 2nd array coordinate = x   (Caution: y,x and not x,y)
//   *start / *max  : domain bounds
//   is/ie, ks/ke   : this rank's working bounds (here = domain bounds)
//   *outer         : extent including boundaries
//   *inner         : extent excluding boundaries (outer - 2*halo)
int main(int argc, char *argv[])
{
  int i, k, it;
  int const istart = 0, kstart = 0;
  int const halo = 1;           // boundary/halo width
  double const eps = 1.e-08;

  // Parameters (with defaults):
  //   argv[1] = grid size (imax == kmax)
  //   argv[2] = itmax
  //   argv[3] = print_result (1 or 0 for printing / not printing the result)
  int imax = 80, kmax = 80, itmax = 20000;
  int print_result = 0;
  if (argc > 1) {
    imax = kmax = atoi(argv[1]);
    if (imax <= 0) {
      fprintf(stderr, "invalid grid size: %s\n", argv[1]);
      return 1;
    }
  }
  if (argc > 2) {
    itmax = atoi(argv[2]);
    if (itmax <= 0) {
      fprintf(stderr, "invalid itmax: %s\n", argv[2]);
      return 1;
    }
  }
  if (argc > 3) {
    print_result = atoi(argv[3]);
  }

  int const is = istart, ie = imax, ks = kstart, ke = kmax;
  int const iouter = ie - is + 1, kouter = ke - ks + 1;
  int const iinner = iouter - 2 * halo, kinner = kouter - 2 * halo;

  // index range: (is:ie, ks:ke)
  double *phi = malloc(iouter * kouter * sizeof(double));
  // index range: (is+halo:ie-halo, ks+halo:ke-halo)
  double *phin = malloc(iinner * kinner * sizeof(double));
  if (phi == NULL || phin == NULL) {
    fprintf(stderr, "allocation failed\n");
    free(phi);
    free(phin);
    return 1;
  }

  double dx = 1.e0 / (kmax - kstart);   // grid spacing in x
  double dy = 1.e0 / (imax - istart);   // grid spacing in y
  double dx2 = dx * dx;
  double dy2 = dy * dy;
  double dx2i = 1.e0 / dx2;             // 1/dx² (precomputed)
  double dy2i = 1.e0 / dy2;             // 1/dy² (precomputed)
  double dt = dmin(dx2, dy2) / STABILITY_FACTOR;   // time step

  // start values 0
  for (i = imax_i(1, is); i <= imin(ie, imax - 1); i++) { // do not overwrite the boundary condition
    for (k = ks; k <= imin(ke, kmax - 1); k++) {
      phi[idx(i, k, is, ks, kouter)] = 0.e0;
    }
  }

  // start values 1
  if (ke == kmax) {
    for (i = is; i <= ie; i++) {
      phi[idx(i, kmax, is, ks, kouter)] = 1.e0;
    }
  }

  // start values dx
  if (is == istart) {
    for (k = ks; k <= imin(ke, kmax - 1); k++) {
      phi[idx(istart, k, is, ks, kouter)] = (k - kstart) * dx;
    }
  }
  if (ie == imax) {
    for (k = ks; k <= imin(ke, kmax - 1); k++) {
      phi[idx(imax, k, is, ks, kouter)] = (k - kstart) * dx;
    }
  }

  printf("\nHeat Conduction 2d\n");
  printf("\ndx =%12.4e, dy =%12.4e, dt=%12.4e, eps=%12.4e\n", dx, dy, dt, eps);

  struct timespec t_start, t_end;
  clock_gettime(CLOCK_MONOTONIC, &t_start);

  // iteration
  for (it = 1; it <= itmax; it++) {
    double dphimax = 0.;   // max change this sweep (convergence test)
    for (i = is + halo; i <= ie - halo; i++) {
      for (k = ks + halo; k <= ke - halo; k++) {
        double dphi = (phi[idx(i + 1, k, is, ks, kouter)] + phi[idx(i - 1, k, is, ks, kouter)]
                       - 2. * phi[idx(i, k, is, ks, kouter)]) * dy2i
                    + (phi[idx(i, k + 1, is, ks, kouter)] + phi[idx(i, k - 1, is, ks, kouter)]
                       - 2. * phi[idx(i, k, is, ks, kouter)]) * dx2i;
        dphi = dphi * dt;
        dphimax = dmax(dphimax, dphi);
        phin[idxn(i, k, is, ks, kinner)] = phi[idx(i, k, is, ks, kouter)] + dphi;
      }
    }

    // save values
    for (i = is + halo; i <= ie - halo; i++) {
      for (k = ks + halo; k <= ke - halo; k++) {
        phi[idx(i, k, is, ks, kouter)] = phin[idxn(i, k, is, ks, kinner)];
      }
    }

    if (dphimax < eps) break;
  }

  clock_gettime(CLOCK_MONOTONIC, &t_end);
  double wall_sec = (t_end.tv_sec - t_start.tv_sec)
                  + (t_end.tv_nsec - t_start.tv_nsec) * 1.e-9;

  if (print_result) {
    printf("   i=");
    for (i = is; i <= ie; i++)
      printf(" %4d", i);
    printf("\n");
    for (k = ks; k <= ke; k++) {
      printf("k=%3d", k);
      for (i = is; i <= ie; i++)
        printf(" %4.2f", phi[idx(i, k, is, ks, kouter)]);
      printf("\n");
    }
  }

  printf("\n%d iterations\n", it);
  // machine-readable summary row (identical format in heat-mpi-big.c):
  // grid_size,nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s
  printf("RESULT,%d,%d,%d,%d,%d,%.6g,%.6g,%.6g\n",
         imax, 1, 1, 1, it, wall_sec, 0.0, 0.0);

  free(phi);
  free(phin);
  return 0;
}
