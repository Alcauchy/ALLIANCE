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

#define CHI_EM 3
#define CHI_EL 1
#define CHI_PHI 0
#define CHI_A 1
#define CHI_B 2
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
                            out[ind6D] = space_iKz[iz] * var_var.vT[is] *
                                         (space_sqrtM[im + 1] * in[ind6DPlus] + space_sqrtM[im] * in[ind6DMinus]);

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
                        /*treating left boundary*/
                        ind6D = get_flat_c(is, il, 0, ix, iy, iz);
                        ind6DPlus = get_flat_c(is, il, 1, ix, iy, iz);
                        indBoundary = ix * array_local_size.nky * array_local_size.nkz * array_local_size.nl *
                                      array_local_size.ns
                                      + iy * array_local_size.nkz * array_local_size.nl * array_local_size.ns
                                      + iz * array_local_size.nl * array_local_size.ns
                                      + il * array_local_size.ns
                                      + is;
                        out[ind6D] = space_iKz[iz] * var_var.vT[is] *
                                (space_sqrtM[1] * in[ind6DPlus] + space_sqrtM[0] * minus_boundary[indBoundary]);
                        /*treating right boundary*/
                        ind6D = get_flat_c(is, il, array_local_size.nm - 1, ix, iy, iz);
                        ind6DMinus = get_flat_c(is, il, array_local_size.nm - 2, ix, iy, iz);
                        out[ind6D] = space_iKz[iz] * var_var.vT[is] * (space_sqrtM[array_local_size.nm - 1] * in[ind6DMinus] +
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
 * equation_computeGradChi(size_t ix, size_t iy, size_t iz, size_t is, COMPLEX *out_dhdhx, COMPLEX *out_dhdy)
 ***************************************/
void equation_computeGradChi(size_t ix, size_t iy, size_t iz, size_t is, COMPLEX *out_dChidx, COMPLEX *out_dChidy){
    size_t ind4D;
    size_t indChiBuf;
    switch(systemType){
        case(ELECTROSTATIC):
            ind4D = getIndChi(ix,iy,iz,is);
            indChiBuf = getIndChiBufEL_c(ix,iy,iz,is);
            out_dChidx[indChiBuf] = space_iKx[ix] * fields_chi.phi[ind4D];
            out_dChidy[indChiBuf] = space_iKy[iy] * fields_chi.phi[ind4D];
            break;
        case(ELECTROMAGNETIC):
            ind4D = getIndChi(ix,iy,iz,is);
            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is, CHI_PHI);
            out_dChidx[indChiBuf] = space_iKx[ix] * fields_chi.phi[ind4D];
            out_dChidy[indChiBuf] = space_iKy[iy] * fields_chi.phi[ind4D];

            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is, CHI_A);
            out_dChidx[indChiBuf] = space_iKx[ix] * fields_chi.A[ind4D];
            out_dChidy[indChiBuf] = space_iKy[iy] * fields_chi.A[ind4D];

            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is, CHI_B);
            out_dChidx[indChiBuf] = space_iKx[ix] * fields_chi.B[ind4D];
            out_dChidy[indChiBuf] = space_iKy[iy] * fields_chi.B[ind4D];
            break;
        default:
            printf("[MPI process %d] Error computing gradient chi! aborting...",mpi_my_rank);
            exit(1);
    }


}


/***************************************
 * equation_getNonlinearTerm()
 ***************************************/
void equation_getNonlinearTerm(COMPLEX *h, COMPLEX *plus_boundary, COMPLEX *minus_boundary, COMPLEX *out) {
    size_t ind6D;
    size_t ind4D;
    size_t indPhi;
    size_t indB;
    size_t indA;
    size_t chi_buf_size;
    size_t chi_buf_size_r;
    switch(systemType){
        case ELECTROSTATIC:
            chi_buf_size = array_local_size.nkx *
                    array_local_size.nky *
                    array_local_size.nkz *
                    array_local_size.ns *
                    CHI_EL;
            chi_buf_size_r = array_local_size.nkx *
                             array_local_size.nky *
                             array_local_size.nz *
                             array_local_size.ns *
                             CHI_EL;
            break;
        case ELECTROMAGNETIC:
            chi_buf_size = array_local_size.nkx *
                    array_local_size.nky *
                    array_local_size.nkz *
                    array_local_size.ns *
                    CHI_EM;
            chi_buf_size_r = array_local_size.nkx *
                           array_local_size.nky *
                           array_local_size.nz *
                           array_local_size.ns *
                           CHI_EM;
            break;
    }
    COMPLEX *dchidx = malloc(chi_buf_size * sizeof(*dchidx));
    COMPLEX *dchidy = malloc(chi_buf_size * sizeof(*dchidy));
    double *dchidx_r = malloc(array_local_size.total_real * sizeof(*dchidx_r));
    double *dchidy_r = malloc(array_local_size.total_real * sizeof(*dchidy_r));

    COMPLEX *dhdx = malloc(array_local_size.total_comp * sizeof(*dhdx));
    COMPLEX *dhdy = malloc(array_local_size.total_comp * sizeof(*dhdy));

    double *dhdx_r = malloc(array_local_size.total_real * sizeof(*dhdx_r));
    double *dhdy_r = malloc(array_local_size.total_real * sizeof(*dhdy_r));
    /*computing gradients simultaneously to save time*/
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            /*computing h gradients*/
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            dhdx[ind6D] = space_iKx[ix] * h[ind6D];
                            dhdy[ind6D] = space_iKy[iy] * h[ind6D];
                            /*computing chi gradients*/
                            equation_computeGradChi(ix, iy, iz, is, dchidx, dchidy);
                        }
                    }
                }
            }
        }
    }
    //
    //make fftw of everything
    //
    fftw_c2r(dhdx, dhdx_r);
    fftw_c2r(dhdy,dhdy_r);
    fftw_c2r_chi(dchidx,dchidx_r);
    fftw_c2r_chi(dchidy,dchidy_r);
    //
    //compute product
    //

    //
    // make ifftw of rhs
    //

    //
    // free memory
    //
    free(dchidx);
    free(dchidy);
    free(dchidx_r);
    free(dchidy_r);
    free(dhdx);
    free(dhdy);
    free(dhdx_r);
    free(dhdy_r);
};

/***************************************
 * equation_getRHS()
 ***************************************/
void equation_getRHS(const COMPLEX *in, COMPLEX *out) {
    /* allocating memory for buffers */
    size_t boundary_size = array_local_size.nkx *
            array_local_size.nky *
            array_local_size.nkz *
            array_local_size.nl *
            array_local_size.ns;
    COMPLEX *rhs = malloc(array_local_size.total_comp * sizeof(*rhs));
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *minus_boundary = calloc(boundary_size, sizeof(*minus_boundary));
    COMPLEX *plus_boundary = calloc(boundary_size, sizeof(*plus_boundary));
    /* exchanging boundaries to compute linear term and fields */
    /* computing h */
    fields_sendG(in);
    fields_getFields(g00, g10, g01);
    fields_getChi();
    distrib_getH(h, in);

    /* boundary exchange */
    mpi_exchangeMBoundaries(h, plus_boundary, minus_boundary);

    /* computing linear term */
    equation_getLinearTerm(h, plus_boundary, minus_boundary, rhs);

    /* computing nonlinear term */
    equation_getNonlinearTerm(h, plus_boundary, minus_boundary, rhs);
    for (size_t i = 0; i < array_local_size.total_comp; i++){
        out[i] = rhs[i];
    }
    free(rhs);
    free(h);
    free(minus_boundary);
    free(plus_boundary);
};