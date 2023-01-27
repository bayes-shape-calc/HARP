#ifndef HARP_RM_C_H_
#define HARP_RM_C_H_

// External calls from python
void render_atoms_gaussian(int natoms, double *dmodel, int *dboundsijk, double *weights, double *origin, double *dxyz, long *nxyz, double *xyz, double *sigmas, int nsigma_cutoff, double offset);
void harp(int natoms, int nadfs, int nblobs, double *ln_out, double *density, double *com, double *origin, double *dxyz, long *nxyz, double *xyz, double *weights, double *blobs, double *adfs, int nsigma_cutoff, double offset);
double ln_evidence(double *x, double *y, long *nxyz, double lnprior_loc, double lnprior_scale);

#endif
