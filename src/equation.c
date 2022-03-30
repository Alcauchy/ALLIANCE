////////////////////////////////////////////////////////////////////////////////
// 26/01/2022 created by Gene Gorbunov
//                                   EQUATION
//
// equation_getLinearTerm
// equation_getNonlinearElectromagnetic
// equation_getNonlinearElectrostatic
// equation_getNonlinearTerm
// equation_getRHS
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "equation.h"
#include "array.h"
#include "fields.h"
#include "distrib.h"

/***************************************
 * equation_getLinearTerm()
 ***************************************/
void equation_getLinearTerm(const COMPLEX *in, const COMPLEX *plus_boundary, const COMPLEX *minus_boundary, COMPLEX *out) {
    size_t ind6D;
    size_t ind6DPlus;
    size_t ind6DMinus;
    size_t indBoundary;
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 1; im < array_local_size.nm - 1; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ind6D = get_flat_c(is, il, im, ix, iy, iz);
                            ind6DPlus = get_flat_c(is, il, im + 1, ix, iy, iz);
                            ind6DMinus = get_flat_c(is, il, im - 1, ix, iy, iz);
                            //printf("in = %f\n", cabs(in[ind6DPlus]));
                            out[ind6D] = space_iKz[iz] *
                                         (space_sqrtM[im + 1] * in[ind6DPlus] + space_sqrtM[im] * in[ind6DMinus]);

                            //printf("out = %f\n", cabs(out[ind6D]));
                        }
                    }
                }
            }
        }
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t il = 0; il < array_local_size.nl; il++) {
                    for (size_t is = 0; is < array_local_size.ns; is++) {
                        ind6D = get_flat_c(is, il, 0, ix, iy, iz);
                        ind6DPlus = get_flat_c(is, il, 1, ix, iy, iz);
                        indBoundary = ix * array_local_size.nky * array_local_size.nkz * array_local_size.nl *
                                      array_local_size.ns
                                      + iy * array_local_size.nkz * array_local_size.nl * array_local_size.ns
                                      + iz * array_local_size.nl * array_local_size.ns
                                      + il * array_local_size.ns
                                      + is;
                        out[ind6D] = space_iKz[iz] *
                                     (space_sqrtM[1] * in[ind6DPlus] + space_sqrtM[0] * minus_boundary[indBoundary]);

                        ind6D = get_flat_c(is, il, array_local_size.nm - 1, ix, iy, iz);
                        ind6DMinus = get_flat_c(is, il, array_local_size.nm - 1, ix, iy, iz);
                        out[ind6D] = space_iKz[iz] * (space_sqrtM[array_local_size.nm - 2] * in[ind6DMinus] +
                                                      space_sqrtM[array_local_size.nm] * plus_boundary[indBoundary]);
                    }
                }
            }
        }
    }

};

/***************************************
 * equation_getNonlinearElectromagnetic()
 ***************************************/
void equation_getNonlinearElectromagnetic() {};

/***************************************
 * equation_getNonlinearElectrostatic()
 ***************************************/
void equation_getNonlinearElectrostatic() {};

/***************************************
 * equation_getNonlinearTerm()
 ***************************************/
void equation_getNonlinearTerm() {};

/***************************************
 * equation_getRHS()
 ***************************************/
void equation_getRHS(const COMPLEX *in, COMPLEX *out) {
    /* allocating memory for buffers */
    COMPLEX *rhs = malloc(array_local_size.total_comp * sizeof(*rhs));
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *minus_boundary = calloc(array_local_size.total_comp, sizeof(*minus_boundary));
    COMPLEX *plus_boundary = calloc(array_local_size.total_comp, sizeof(*plus_boundary));
    /* exchanging boundaries to compute linear term and fields */
    /* computing h */
    fields_sendG(in);
    fields_getFields(g00, g10, g01);
    fields_getChi();
    distrib_getH(h, in);

    /* border exchange */
    mpi_exchangeMBoundaries(h, plus_boundary, minus_boundary);
    /* computing linear term */
    equation_getLinearTerm(h, plus_boundary, minus_boundary, rhs);
    /* computing nonlinear term */
    equation_getNonlinearTerm();
    for (size_t i = 0; i < array_local_size.total_comp; i++){
        out[i] = rhs[i];
    }
    free(rhs);
    free(h);
    free(minus_boundary);
    free(plus_boundary);
};