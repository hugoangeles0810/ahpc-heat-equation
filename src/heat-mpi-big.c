#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// CFL stability factor: dt = min(dx2,dy2) / STABILITY_FACTOR
#define STABILITY_FACTOR 4.e0

static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline int    imin(int a, int b)       { return a < b ? a : b; }
static inline int    imax_i(int a, int b)     { return a > b ? a : b; }

// Linear index into phi, the full local grid including halos/boundaries: (is:ie, ks:ke).
static inline int idx(int i, int k, int is, int ks, int kouter) {
  return (i - is) * kouter + (k - ks);
}

// Linear index into phin, the inner grid excluding halos/boundaries:
// (is+halo:ie-halo, ks+halo:ke-halo).
static inline int idxn(int i, int k, int is, int ks, int kinner) {
  return (i - is - 1) * kinner + (k - ks - 1);
}

// Index convention:
//   i = 1st array coordinate = 1st process coordinate = y
//   k = 2nd array coordinate = 2nd process coordinate = x  (Caution: y,x and not x,y)
//   *outer = extent including halos; *inner = extent excluding halos (outer - 2*halo)
int main(int argc, char *argv[])
{
  int i, k, it;
  int const print_halo = 1;      // 1/0: include halos in the print
  int const istart = 0, kstart = 0;
  int const halo = 1;            // boundary/halo width
  double const eps = 1.e-08;
  double *phi, *phin;

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

  // variables for parallelization
  int numprocs, my_rank, right, left, upper, lower;
  MPI_Comm comm;          // Cartesian communicator
  int dims[2], coords[2]; // only for MPI_Cart_create and ..._coords
  int period[2];
  int idim, kdim, icoord, kcoord;
  int isize, iinner0, in1,  ksize, kinner0, kn1; // only for calculation of iinner, ...
  int iinner, iouter, is, ie, kinner, kouter, ks, ke;
  int ic, i_first, i_last,  kc, k_first, k_last; // only for printing
  int const stride = 10; // for reduced calculation of the abort criterion
  MPI_Request req[4];
  MPI_Status statuses[4];
  int gsizes[2], lsizes[2], starts[2]; // only for MPI_Type_create_subarray
  MPI_Datatype horizontal_border, vertical_border;
  double dphimaxpartial; // for local/global calculation of the abort criterion
  double start_time, end_time, comm_time, criterion_time;

  // index ranges:  i = is..ie,  icoord = 0..idim-1
  //                k = ks..ke,  kcoord = 0..kdim-1   (s = start index, e = end index)

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // the 2-dimensional domain decomposition:
  dims[0] = 0; dims[1] = 0;
  period[0] = 0; period[1] = 0;
  MPI_Dims_create(numprocs, 2, dims);
  idim = dims[0];
  kdim = dims[1];
  MPI_Cart_create(MPI_COMM_WORLD,2,dims,period,1,&comm);
  MPI_Comm_rank(comm, &my_rank);
  MPI_Cart_coords(comm, my_rank, 2, coords);
  icoord = coords[0];
  kcoord = coords[1];
  MPI_Cart_shift(comm, 0, 1, &left, &right);
  // ranks (left,right) are coordinates ((icoord-1,kcoord),(icoord+1,kcoord))
  MPI_Cart_shift(comm, 1, 1, &lower,&upper);
  // ranks (lower,upper) are coordinates ((icoord,kcoord-1),(icoord,kcoord+1))

// Step 1: calculating the own subdomain in each process
// --------------------------------------------------------------
//
// whole y indecees   |------------- isize = imax-istart+1 ------------|
//    start/end index  ^-istart                                  imax-^
//
// 1. interval        |--- iouter1 --|
//                    |--|--------|--|
//                   halo iinner1 halo
//    start/end index  ^-is      ie-^
// 2. interval                 |--|--------|--|
// 3. interval                          |--|--------|--|
// 4. interval                                   |--|-------|--|
//                                                   iinner0 = iinner1 - 1
// 5. interval = idim's interval                         |--|-------|--|
//
// In each iteration on each interval, the inner area is computed
// by using the values of the last iteration in the whole outer area.
//
// To fit exactly into isize, we use in1 intervals with iinner1 = iinner0 + 1
// and (idim-in1) intervals with iinner0:
//   isize = imax-istart+1 = 2*halo + in1*iinner1 + (idim-in1)*iinner0

  isize = imax - istart + 1; // total number of elements, including the "2*halo" boundary elements
  iinner0 = (isize - 2*halo)  / idim; // smaller inner size through divide with rounding off
  in1 = isize - 2*halo - idim * iinner0; // number of processes that must have "iinner0+1" unknowns
  if (icoord < in1) {  // the first in1 processes will have "iinner0+1" unknowns
    iinner = iinner0 + 1;
    is = (istart+halo) + icoord * iinner - halo; // "is" is the position of the first halo or
                                                 // boundary element of the subdomain
  } else {             // and all other processes will have iinner0 unknowns
    iinner = iinner0;
    is = istart + in1 * (iinner0+1) + (icoord-in1) * iinner;
  }
  iouter = iinner + 2*halo;
  ie = is + iouter - 1;

  // same for x coordinate:
  ksize = kmax - kstart + 1;
  kinner0 = (ksize - 2*halo)  / kdim;
  kn1 = ksize - 2*halo - kdim * kinner0;
  if (kcoord < kn1) {
    kinner = kinner0 + 1;
    ks = (kstart+halo) + kcoord * kinner - halo;
  } else {
    kinner = kinner0;
    ks = kstart + kn1 * (kinner0+1) + (kcoord-kn1) * kinner;
  }
  kouter = kinner + 2*halo;
  ke = ks + kouter - 1;

  if (my_rank == 0) {
    printf("\n isize=%3d idim=%3d || in1*(iinner0+1)=%3d *%3d + (idim-in1)*iinner0=%3d *%3d + 2*halo=2*%1d || sum = %3d\n",
               isize,    idim,       in1, iinner0+1,             idim-in1, iinner0,             halo,
                                                                            in1*(iinner0+1) + (idim-in1)*iinner0 + 2*halo);
    printf("\n ksize=%3d kdim=%3d || kn1*(kinner0+1)=%3d *%3d + (kdim-kn1)*kinner0=%3d *%3d + 2*halo=2*%1d || sum = %3d\n",
               ksize,    kdim,       kn1, kinner0+1,             kdim-kn1, kinner0,             halo,
                                                                            kn1*(kinner0+1) + (kdim-kn1)*kinner0 + 2*halo);
  }

  // 'phin' must not be an empty array, i.e. iinner must be at least 1 in all processes
  // (total unknowns in direction "i" at least 1*idim, and same for k), otherwise the halo
  // communication would not work. All processes must make the same decision, so they use
  // "isize,ksize" and "idim,kdim" rather than only their own "iinner,kinner".
  if(((isize - 2*halo) < idim) ||
     ((ksize - 2*halo) < kdim)) {
    if(my_rank == 0) {
        printf("phin is in some processes an empty array because isize-2*halo=%3d < idim=%3d or ksize-2*halo=%3d < kdim=%3d\n",
                                                                 isize-2*halo,     idim,       ksize-2*halo,     kdim);
    }
    MPI_Finalize();
    exit(0);
  }

  // index range: (is:ie, ks:ke)
  phi  = malloc(iouter * kouter * sizeof(double));
  // index range: (is+halo:ie-halo, ks+halo:ke-halo)
  phin = malloc(iinner * kinner * sizeof(double));
  if (phi == NULL || phin == NULL) {
    fprintf(stderr, "rank %d: allocation failed\n", my_rank);
    free(phi);
    free(phin);
    MPI_Abort(comm, 1);
  }

  // create and commit the derived datatypes for the halo borders
  gsizes[0]=iouter; gsizes[1]=kouter;
  lsizes[0]=halo;   lsizes[1]=kinner;  starts[0]=0; starts[1]=0;
  MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &vertical_border);
  MPI_Type_commit(&vertical_border);
  gsizes[0]=iouter; gsizes[1]=kouter;
  lsizes[0]=iinner; lsizes[1]=halo;    starts[0]=0; starts[1]=0;
  MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &horizontal_border);
  MPI_Type_commit(&horizontal_border);

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

  if (my_rank == 0) {
    printf("\nHeat Conduction 2d\n");
    printf("\ndx =%12.4e, dy =%12.4e, dt=%12.4e, eps=%12.4e\n", dx, dy, dt, eps);
  }
  start_time = MPI_Wtime();
  comm_time = 0;
  criterion_time = 0;

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

    // Step 2: common abort criterion for all processes.
    // For optimization, only allreduce every "stride"th iteration.
    criterion_time = criterion_time - MPI_Wtime();
    if ((it % stride) == 0) {
      if (numprocs > 1) {
        // each process only calculates a local (not global) dphimax
        dphimaxpartial = dphimax;
        MPI_Allreduce(&dphimaxpartial, &dphimax, 1, MPI_DOUBLE, MPI_MAX, comm);
      }
      if (dphimax < eps) {
        criterion_time = criterion_time + MPI_Wtime();
        break;
      }
    }
    criterion_time = criterion_time + MPI_Wtime();

    // Step 3: the halo communication
    comm_time = comm_time - MPI_Wtime();
    // send and receive to/from upper/lower
    if (kdim > 1) { // otherwise both lower and upper are MPI_PROC_NULL in all processes
      MPI_Irecv(&phi[idx(is+halo, ks,      is, ks, kouter)], 1, horizontal_border, lower, 1, comm, &req[0]);
      MPI_Irecv(&phi[idx(is+halo, ke,      is, ks, kouter)], 1, horizontal_border, upper, 2, comm, &req[1]);
      MPI_Isend(&phi[idx(is+halo, ke-halo, is, ks, kouter)], 1, horizontal_border, upper, 1, comm, &req[2]);
      MPI_Isend(&phi[idx(is+halo, ks+halo, is, ks, kouter)], 1, horizontal_border, lower, 2, comm, &req[3]);
      MPI_Waitall(4, req, statuses);
    }

    // send and receive to/from left/right
    if (idim > 1) { // otherwise both left and right are MPI_PROC_NULL in all processes
      MPI_Irecv(&phi[idx(is,      ks+halo, is, ks, kouter)], 1, vertical_border, left,  3, comm, &req[0]);
      MPI_Irecv(&phi[idx(ie,      ks+halo, is, ks, kouter)], 1, vertical_border, right, 4, comm, &req[1]);
      MPI_Isend(&phi[idx(ie-halo, ks+halo, is, ks, kouter)], 1, vertical_border, right, 3, comm, &req[2]);
      MPI_Isend(&phi[idx(is+halo, ks+halo, is, ks, kouter)], 1, vertical_border, left,  4, comm, &req[3]);
      MPI_Waitall(4, req, statuses);
    }
    comm_time = comm_time + MPI_Wtime();

  }
  end_time = MPI_Wtime();

  if (print_result) {
    for (ic = 0; ic <= idim-1; ic++) {
      for (kc = 0; kc <= kdim-1; kc++) {
        if ((ic == icoord) && (kc == kcoord)) {
          i_first = is ;
          i_last  = ie ;
          k_first = ks ;
          k_last  = ke ;
          if (! print_halo) { // do not print the halos
            if (ic > 0) i_first = is + halo;
            if (ic < idim-1) i_last = ie - halo;
            if (kc > 0) k_first = ks + halo;
            if (kc < kdim-1) k_last = ke - halo;
          }
          if (kc == 0) {
              printf("\nprinting the %3dth horizontal block\n", ic);
              printf("               i=");
              for (i=i_first; i<=i_last; i++) printf(" %4d",i); printf("\n");
          } else {
            if (print_halo)  printf("\n"); // additional empty line between the processes
          }
          for (k=k_first; k<=k_last; k++) {
              printf("ic=%2d kc=%2d k=%3d", ic, kc, k);
              for (i=i_first; i<=i_last; i++) printf(" %4.2f", phi[idx(i, k, is, ks, kouter)]); printf("\n");
          }
        }
        MPI_Barrier(comm); // separate the printing of each block by different processes.
                           // Caution: this works in most cases, but does not guarantee that the
                           //          output lines on the common stdout are in the expected sequence.
      }
    }
  }

  if (my_rank == 0) {
    printf("\n%d iterations\n", it);
    // machine-readable summary row (identical format in heat-big.c):
    // grid_size,nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s
    printf("RESULT,%d,%d,%d,%d,%d,%.6g,%.6g,%.6g\n",
           imax, numprocs, idim, kdim, it,
           end_time - start_time, comm_time, criterion_time);
  }

  MPI_Type_free(&vertical_border);
  MPI_Type_free(&horizontal_border);
  free(phi);
  free(phin);
  MPI_Finalize();
  return 0;
}
