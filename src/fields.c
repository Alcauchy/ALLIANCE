/**************************************
* @file fields.c
*   \brief field computation and manipulation module
*
*   Rerquired to compute
*   \f$A_{||}(\mathbf{k}), B_{||}(\mathbf{k}), \phi(\mathbf{k})\f$
*   potentials, as well as gyrokinetic potentials
*   \f$\chi^{A}_s(\mathbf{k}),\chi^{B}_s(\mathbf{k}),\chi^{\phi}_s(\mathbf{k})\f$
*
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 26/01/2022 created by Gene Gorbunov
//                                   FIELDS
//
// fields_init
// fields_getA
// fields_getB
// fields_getPhi
// fields_getFields
// fields_getChi
// fields_getChiPhi
// fields_getChiB
// fields_getChiA
// fields_sendF
// fields_getFieldsFromH
// fields_getAFromH
// fields_getBFromH
// fields_getPhiFromH
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "fields.h"

#define CHI_PHI 0
#define CHI_A 1
#define CHI_B 2

struct fields_fields fields_fields;
struct fields_chi fields_chi;

/* parameters needed to compute A_parallel:
 * divisor (see formula for A_parallel) and q*n*vT*J_00*/
double *A_denom;
double *qnvTsJ;

/* parameters needed to compute B_parallel:
 * I_B = n * T / B0 * J_10 */
double *I_B;

/* parameters needed to compute phi:
 * I_phi = q * n * J_00 */
double *I_phi;

/* parameters which are common for both phi and B_parallel:
 * a_pot = Sum_s q^2 * n / T * (1-J_00^2)
 * b_pot = 1 + beta * Sum_s n * T / B0^2 * J10^2
 * c_pot = Sum_s q * n / B0 * J10 * J00
 * phiB_denom = 2 * a_pot * b_pot + beta * c^2 */
double *a_pot;
double *b_pot;
double *c_pot;
double *phiB_denom;

int *global_nm_index;
/* buffers to send and store g^0_{s0},g^1_{s0},g^0_{s1} to compute fields*/
COMPLEX *f00;
COMPLEX *f10;
COMPLEX *f01;

/***************************************
 * \fn void fields_init():
 * \brief intializes fields
 *
 * pre-computes some comstants required to compute fields.
 * Called in #init_start function
 ***************************************/
void fields_init() {
    global_nm_index = malloc(array_local_size.nm * sizeof(*global_nm_index));
    for (size_t i = 0; i < array_local_size.nm; i++){
        global_nm_index[i] = array_global_size.nm / mpi_dims[0] * mpi_my_col_rank + i;
    }
    switch (systemType)
    {
        case ELECTROSTATIC:
            fields_fields.phi = malloc(array_local_size.nkx *
                                       array_local_size.nky *
                                       array_local_size.nkz *
                                       sizeof(*fields_fields.phi));
            fields_fields.A = 0;
            fields_fields.B = 0;

            fields_chi.phi = malloc(array_local_size.nkx *
                                    array_local_size.nky *
                                    array_local_size.nkz *
                                    array_local_size.ns *
                                    sizeof(*fields_chi.phi));
            fields_chi.A = 0;
            fields_chi.B = 0;

            f00 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*f00));
            f10 = 0;
            f01 = 0;
            break;

        case ELECTROMAGNETIC:
            fields_fields.phi = malloc(array_local_size.nkx *
                                       array_local_size.nky *
                                       array_local_size.nkz *
                                       sizeof(*fields_fields.phi));
            fields_fields.A = malloc(array_local_size.nkx *
                                     array_local_size.nky *
                                     array_local_size.nkz *
                                     sizeof(*fields_fields.A));
            fields_fields.B = malloc(array_local_size.nkx *
                                     array_local_size.nky *
                                     array_local_size.nkz *
                                     sizeof(*fields_fields.B));

            fields_chi.phi = malloc(array_local_size.nkx *
                                    array_local_size.nky *
                                    array_local_size.nkz *
                                    array_local_size.ns *
                                    sizeof(*fields_chi.phi));
            fields_chi.A = malloc(array_local_size.nkx *
                                    array_local_size.nky *
                                    array_local_size.nkz *
                                    array_local_size.ns *
                                    sizeof(*fields_chi.A));
            fields_chi.B = malloc(array_local_size.nkx *
                                    array_local_size.nky *
                                    array_local_size.nkz *
                                    array_local_size.ns *
                                    sizeof(*fields_chi.B));

            f00 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*f00));
            f10 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*f10));
            f01 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*f01));

            /* preparing auxiliary data to compute A_parallel*/
            A_denom = malloc(array_local_size.nkx *
                             array_local_size.nky * sizeof(*A_denom));
            qnvTsJ = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.ns * sizeof(*qnvTsJ));
            size_t flatInd;
            double *sum = malloc(array_local_size.nkx * array_local_size.nky * sizeof(*sum));
            for(size_t ix  = 0 ; ix < array_local_size.nkx; ix++)
            {
                for(size_t iy  = 0 ; iy < array_local_size.nky; iy++)
                {
                    flatInd = ix * array_local_size.nky + iy;
                    sum[flatInd] = 0.;
                    for (size_t is = 0; is < array_local_size.ns; is++)
                    {
                        sum[flatInd] += var_var.q[is] * var_var.q[is] * var_var.n[is] *
                                        var_var.vT[is] * var_var.vT[is] /
                                        var_var.T[is] * var_J0[var_getJIndex(ix,iy,is)] *
                                                        var_J0[var_getJIndex(ix,iy,is)];
                    }
                    sum[flatInd] *= var_var.beta/4.;
                    if (fabs((space_kPerp2[flatInd] + sum[flatInd])) > 1e-16){
                        A_denom[flatInd] = 1.0/(space_kPerp2[flatInd] + sum[flatInd]);
                    }
                    else{
                        A_denom[flatInd] = 1e16;
                    }
                    //printf("A_denom[flatInd] = %f\n", A_denom[flatInd]);
                }
            }
            for(size_t ix  = 0 ; ix < array_local_size.nkx; ix++)
            {
                for(size_t iy  = 0 ; iy < array_local_size.nky; iy++)
                {
                    for(size_t is  = 0 ; is < array_local_size.ns; is++)
                    {
                        qnvTsJ[var_getJIndex(ix,iy,is)] = 0.5 * sqrt(0.5) * var_var.beta * var_var.q[is] *
                                                          var_var.n[is] * var_var.vT[is] *
                                                          var_J0[var_getJIndex(ix,iy,is)];
                    }
                }
            }


            /* preparing auxiliary data to compute Phi */
            I_phi = malloc(array_local_size.nkx *
                            array_local_size.nky *
                            array_local_size.ns * sizeof(*I_phi));
            size_t ind3D;
            for (size_t ix = 0; ix < array_local_size.nkx; ix++)
            {
                for (size_t iy = 0; iy < array_local_size.nky; iy++)
                {
                    for (size_t is = 0; is < array_local_size.ns; is++)
                    {
                        ind3D = ix * array_local_size.nky  * array_local_size.ns +
                                iy  * array_local_size.ns +
                                is;
                        I_phi[ind3D] = var_var.q[is] * var_var.n[is] * var_J0[var_getJIndex(ix,iy,is)];
                    }
                }
            }

            /* preparing auxiliary data to compute B_parallel */
            I_B = malloc(array_local_size.nkx *
                           array_local_size.nky *
                           array_local_size.ns * sizeof(*I_B));
            for (size_t ix = 0; ix < array_local_size.nkx; ix++)
            {
                for (size_t iy = 0; iy < array_local_size.nky; iy++)
                {
                    for (size_t is = 0; is < array_local_size.ns; is++)
                    {
                        ind3D = ix * array_local_size.nky  * array_local_size.ns +
                                iy  * array_local_size.ns +
                                is;
                        I_B[ind3D] = var_var.n[is] * var_var.T[is] * var_J1[var_getJIndex(ix,iy,is)];
                        I_B[ind3D] /= var_var.B0;
                    }
                }
            }
            /* preparing common auxiliary data for phi and B_parallel */
            a_pot = malloc(array_local_size.nkx *
                           array_local_size.nky * sizeof(*a_pot));
            b_pot = malloc(array_local_size.nkx *
                           array_local_size.nky * sizeof(*b_pot));
            c_pot = malloc(array_local_size.nkx *
                           array_local_size.nky * sizeof(*c_pot));
            phiB_denom = malloc(array_local_size.nkx *
                           array_local_size.nky * sizeof(*phiB_denom));
            size_t ind2D;
            for (size_t ix = 0; ix < array_local_size.nkx; ix++)
            {
                for (size_t iy = 0; iy < array_local_size.nky; iy++)
                {
                    ind2D = ix * array_local_size.nky + iy;
                    a_pot[ind2D] = 0;
                    b_pot[ind2D] = 1.;
                    c_pot[ind2D] = 0;
                    for (size_t is = 0; is < array_local_size.ns; is++)
                    {
                        a_pot[ind2D] += var_var.q[is] * var_var.q[is] *
                                        var_var.n[is] / var_var.T[is] *
                                        (1. - var_J0[var_getJIndex(ix,iy,is)] *
                                             var_J0[var_getJIndex(ix,iy,is)]);
                        b_pot[ind2D] += var_var.beta * var_var.n[is] *
                                        var_var.T[is] / (var_var.B0 * var_var.B0) *
                                        var_J1[var_getJIndex(ix,iy,is)] * var_J1[var_getJIndex(ix,iy,is)];
                        c_pot[ind2D] += var_var.q[is] * var_var.n[is] / var_var.B0 *
                                        var_J0[var_getJIndex(ix,iy,is)] *
                                        var_J1[var_getJIndex(ix,iy,is)];
                    }
                    phiB_denom[ind2D] = 2. * a_pot[ind2D] * b_pot[ind2D] +
                                             var_var.beta * c_pot[ind2D] * c_pot[ind2D];
                    //printf("phib_denom = %f\n", phiB_denom[ind2D]);
                }
            }
            break;

        default:
            printf("[MPI process %d] mistake defining fields! Aborting...",mpi_my_rank);
            exit(1);
    }
};

/***************************************
 * \fn void fields_getA(const COMPLEX *g)
 * \brief compute A field
 * \param g: 4D complex array (kx,ky,kz,s). Must be first Hermite and zeroth Laguerre moment of modified gyrokinetic distribution function g.
 *
 * computes
 * \f$A_{||}(\mathbf{k})\f$
 * potential from
 * \f$g^1_{s0}\f$ (<tt>g</tt> parameter)
 ***************************************/
void fields_getA(const COMPLEX *g) {
    size_t flatInd;
    size_t flatInd2D;
    switch(kinetic)
    {
        case ADIABATIC:
            break;
        case NONADIABATIC:
            for(size_t ix = 0; ix < array_local_size.nkx; ix++)
            {
                for(size_t iy = 0; iy < array_local_size.nky; iy++)
                {
                    for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                    {
                        fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] = 0.;
                        for(size_t is = 0; is < array_local_size.ns; is++)
                        {
                            flatInd = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                      iy * array_local_size.nkz * array_local_size.ns + iz * array_local_size.ns + is;
                            fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] += qnvTsJ[var_getJIndex(ix,iy,is)] *
                                   g[flatInd];
                        }
                        flatInd2D = ix * array_local_size.nky + iy;
                        if (fabs(space_kPerp2[flatInd2D]>1e-16))
                        {
                            fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] *= A_denom[flatInd2D];
                        }
                        else
                        {
                            fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] *= 0.0;
                        }

                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Wrong call of get_A function! Aborting...",mpi_my_rank);
            exit(1);
    }
};

/***************************************
 * \fn void fields_getB(const COMPLEX* g0, const COMPLEX* g1)
 * \brief computes B potential
 * \param g0: 4D complex array (kx,ky,kz,s). Zeroth Hermite and Laguerre moment of modified gyrokinetic distribution function.
 * \param g1: 4D complex array (kx,ky,kz,s). Zeroth Hermite and first Laguerre moment of the modified distribution function.
 *
 * Computes
 * \f$B_{\perp}(\mathbf{k})\f$
 * from
 * \f$g^0_{s0}\f$ (<tt>g0</tt> parameter) and
 * \f$g^1_{s0}\f$ (<tt>g1</tt> parameter).
 ***************************************/
void fields_getB(const COMPLEX* g0, const COMPLEX* g1) {
    size_t ind2D;
    size_t ind3D;
    size_t ind4D;
    switch(kinetic)
    {
        case ADIABATIC:
            break;
        case NONADIABATIC:
            for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                for(size_t iy = 0; iy < array_local_size.nky; iy++){
                    ind2D = ix * array_local_size.nky + iy;
                    if (fabs(phiB_denom[ind2D]) > 1e-16){
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                ind3D = ix * array_local_size.nky  * array_local_size.ns +
                                        iy  * array_local_size.ns +
                                        is;
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] += a_pot[ind2D] * I_B[ind3D] *
                                                                                        (g0[ind4D] + g1[ind4D]);//(g0[ind4D] + g1[ind4D]);
                                fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] += c_pot[ind2D] * I_phi[ind3D] *
                                                                                     g0[ind4D];
                            }
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] *= - var_var.beta / phiB_denom[ind2D];
                        }

                    }
                    else{
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                            /*for(size_t is = 0; is < array_local_size.ns; is++) {
                                ind3D = ix * array_local_size.nky * array_local_size.ns +
                                        iy * array_local_size.ns +
                                        is;
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = - 1.0 * var_var.beta * I_B[ind3D] * // * 0.5 *
                                                                                    (g0[ind4D] + g1[ind4D])/b_pot[ind2D];
                            }*/
                        }
                    }
                }
            }
            break;

        default:
            printf("[MPI process %d] Wrong call of get_B function! Aborting...",mpi_my_rank);
            exit(1);
    }
};

/***************************************
 * \fn void fields_getPhi(const COMPLEX* g0, const COMPLEX* g1)
 * \brief computes phi potential
 * \param g0: 4D complex array (kx,ky,kz,s). Zeroth Hermite and Laguerre moment of modified gyrokinetic distribution function.
 * \param g1: 4D complex array (kx,ky,kz,s). Zeroth Hermite and first Laguerre moment of the modified distribution function.
 *
 * Computes
 * \f$\phi(\mathbf{k})\f$
 * from
 * \f$g^0_{s0}\f$ (<tt>g0</tt> parameter) and
 * \f$g^1_{s0}\f$ (<tt>g1</tt> parameter).
 ***************************************/
void fields_getPhi(const COMPLEX* g0, const COMPLEX* g1) {
    size_t ind2D;
    size_t ind3D;
    size_t ind4D;
    switch(kinetic)
    {
        case ADIABATIC:
            switch(systemType)
            {
                case ELECTROSTATIC:
                    break;
                case ELECTROMAGNETIC:
                    break;
            }
            break;
        case NONADIABATIC:
            switch(systemType)
            {
                case ELECTROSTATIC:
                    break;
                case ELECTROMAGNETIC:
                    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
                    {
                        for(size_t iy = 0; iy < array_local_size.nky; iy++)
                        {
                            ind2D = ix * array_local_size.nky + iy;

                            if (fabs(phiB_denom[ind2D]) > 1e-16)
                            {
                                for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                                {
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                                    for(size_t is = 0; is < array_local_size.ns; is++)
                                    {
                                        ind3D = ix * array_local_size.nky * array_local_size.ns +
                                                iy  * array_local_size.ns +
                                                is;
                                        ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                                iy * array_local_size.nkz * array_local_size.ns +
                                                iz * array_local_size.ns +
                                                is;
                                        fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] += 2. * b_pot[ind2D] *
                                                                                               I_phi[ind3D] *
                                                                                               g0[ind4D];
                                        fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] -= var_var.beta *
                                                                                               c_pot[ind2D] *
                                                                                               I_B[ind3D] *
                                                                                               (g0[ind4D] + g1[ind4D]);
                                    }
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] *= 1.0 / phiB_denom[ind2D];
                                }

                            }
                            else
                            {
                                for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                                {
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                                }
                            }
                        }
                    }
                    break;
            }
            break;
    }
};

/***************************************
 * \fn void fields_getFields(COMPLEX *f00, COMPLEX *f10, COMPLEX *f01)
 * \brief wrapper to get all the fields simultaneously
 * \param g00: 4D complex array (kx,ky,kz,s). Zeroth Hermite and Laguerre moment of modified gyrokinetic distribution function.
 * \param g10: 4D complex array (kx,ky,kz,s). Must be first Hermite and zeroth Laguerre moment of modified gyrokinetic distribution function g.
 * \param g01: 4D complex array (kx,ky,kz,s). Zeroth Hermite and first Laguerre moment of the modified distribution function.
 *
 * Wrapper for functions #fields_getPhi, #fields_getB, #fields_getA.
 ***************************************/
void fields_getFields(COMPLEX *g00, COMPLEX *g10, COMPLEX *g01) {
    fields_getPhi(g00,g01);
    fields_getB(g00,g01);
    fields_getA(g10);
};

/***************************************
 * \fn void fields_getChi():
 * \brief computes gyrokinetic potentials chi
 *
 * Wrapper for functions
 * #fields_getChiPhi,
 * #fields_getChiA,
 * #fields_getChiB
 ***************************************/
void fields_getChi() {
    switch (systemType)
    {
        case ELECTROSTATIC:
            fields_getChiPhi();
            break;
        case ELECTROMAGNETIC:
            fields_getChiPhi();
            fields_getChiA();
            fields_getChiB();
            break;
    }
};

/***************************************
 * \fn fields_getChiPhi()
 * \brief computes chiPhi gyrokinetic potential from phi potential
 *
 * computes \f$chi^{\phi}_s (\mathbf{k})\f$
 ***************************************/
void fields_getChiPhi(){
    size_t ind3D = 0;
    size_t ind4D = 0;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
            {
                for(size_t is = 0; is < array_local_size.ns; is++)
                {
                    ind3D = ix * array_local_size.nky * array_local_size.nkz +
                            iy * array_local_size.nkz +
                            iz;
                    ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                            iy * array_local_size.nkz * array_local_size.ns +
                            iz * array_local_size.ns +
                            is;
                    fields_chi.phi[ind4D] = var_J0[var_getJIndex(ix,iy,is)] * fields_fields.phi[ind3D];
                }
            }
        }
    }
}

/***************************************
 * \fn void fields_getChiB()
 * \brief computes chiB gyrokinetic potential from B potential
 *
 * computes \f$ \chi^{B}_s (\mathbf{k})\f$
 ***************************************/
void fields_getChiB(){
    size_t ind3D;
    size_t ind4D;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
            {
                for(size_t is = 0; is < array_local_size.ns; is++)
                {
                    ind3D = ix * array_local_size.nky * array_local_size.nkz +
                            iy * array_local_size.nkz +
                            iz;
                    ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                            iy * array_local_size.nkz * array_local_size.ns +
                            iz * array_local_size.ns +
                            is;
                    fields_chi.B[ind4D] = var_var.T[is] / (var_var.q[is] * var_var.B0) *
                                          var_J1[var_getJIndex(ix,iy,is)] *
                                          fields_fields.B[ind3D];
                }
            }
        }
    }
}

/***************************************
 * \fn void fields_getChiA()
 * \brief computes chiA gyrokinetic potential from A potential
 *
 * computes \f$chi^{A}_s (\mathbf{k})\f$
 ***************************************/
void fields_getChiA(){
    size_t ind3D = 0;
    size_t ind4D = 0;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
            {
                for(size_t is = 0; is < array_local_size.ns; is++)
                {
                    ind3D = ix * array_local_size.nky * array_local_size.nkz +
                            iy * array_local_size.nkz +
                            iz;
                    ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                            iy * array_local_size.nkz * array_local_size.ns +
                            iz * array_local_size.ns +
                            is;
                    fields_chi.A[ind4D] = - var_var.vT[is] * var_J0[var_getJIndex(ix,iy,is)] * fields_fields.A[ind3D];
                }
            }
        }
    }
}

/***************************************
 * \fn fields_sendF(COMPLEX *f)
 * \param f: complex array. Modified or non-modified gyrokinetic distribution function
 * \brief sends moments of gyrokinetic distribution function which are required to compute fields
 *
 * sends
 * \f$g^{1}_{s0}(\mathbf{k})f\f$,
 * \f$g^{0}_{s1}(\mathbf{k})f\f$,
 * \f$g^{0}_{s0}(\mathbf{k})f\f$
 *
 * to all processes to compute potentials locally.
 ***************************************/
void fields_sendF(COMPLEX *f){
    /* fill buffers with data from processor which has required data */
    size_t ind4D;
    int count;
    int broadcast_root0 = mpi_whereIsM[0];  // processor which works with 0th Hermite moment
    int broadcast_root1 = mpi_whereIsM[2];  // processor which works with 1st Hermite moment
    size_t m0 = mpi_whereIsM[1];            // local index of the 0th Hermite moment
    size_t m1 = mpi_whereIsM[3];            // local index of the 1st Hermite moment
    switch(systemType){
        case ELECTROSTATIC:
            if (mpi_my_col_rank = broadcast_root0){
                for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                    for(size_t iy = 0; iy < array_local_size.nky; iy++){
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                f00[ind4D] = f[get_flat_c(is, 0, m0, ix, iy, iz)];
                            }
                        }
                    }
                }
            }
            MPI_Bcast(f00, count, MPI_DOUBLE_COMPLEX, broadcast_root0, mpi_col_comm);
            break;
        case ELECTROMAGNETIC:
            if (mpi_my_col_rank == broadcast_root0){
                for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                    for(size_t iy = 0; iy < array_local_size.nky; iy++){
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                f00[ind4D] = f[get_flat_c(is, 0, m0, ix, iy, iz)];
                                f01[ind4D] = f[get_flat_c(is, 1, m0, ix, iy, iz)];
                            }
                        }
                    }
                }
            }
            if (mpi_my_col_rank == broadcast_root1){
                for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                    for(size_t iy = 0; iy < array_local_size.nky; iy++){
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                f10[ind4D] = f[get_flat_c(is, 0, m1, ix, iy, iz)];
                            }
                        }
                    }
                }
            }
            count = array_local_size.nkx * array_local_size.nky * array_local_size.nkz * array_local_size.ns;
            MPI_Bcast(f00, count, MPI_DOUBLE_COMPLEX, broadcast_root0, mpi_col_comm);
            MPI_Bcast(f01, count, MPI_DOUBLE_COMPLEX, broadcast_root0, mpi_col_comm);
            MPI_Bcast(f10, count, MPI_DOUBLE_COMPLEX, broadcast_root1, mpi_col_comm);
            break;
    }
};

/***************************************
 * \fn void fields_getFieldsFromH(COMPLEX *h00, COMPLEX *h10, COMPLEX *h01)
 * \brief wrapper to get all the fields simultaneously, computed from gyrokinetic distribution function
 * \param h00: 4D complex array (kx,ky,kz,s). Zeroth Hermite and Laguerre moment of gyrokinetic distribution function.
 * \param h10: 4D complex array (kx,ky,kz,s). Must be first Hermite and zeroth Laguerre moment of gyrokinetic distribution function h.
 * \param h01: 4D complex array (kx,ky,kz,s). Zeroth Hermite and first Laguerre moment of the distribution function.
 *
 * Wrapper for functions #fields_getPhiFromH, #fields_get_BFromH, #fields_getAFromH.
 ***************************************/
void fields_getFieldsFromH(COMPLEX *h00, COMPLEX *h10, COMPLEX *h01){
    fields_getPhiFromH(h00);
    fields_getBFromH(h00, h01);
    fields_getAFromH(h10);
};

/***************************************
 * \fn void fields_getAFromH(const COMPLEX* h)
 * \brief compute A field
 * \param h: 4D complex array (kx,ky,kz,s). Must be first Hermite and zeroth Laguerre moment of gyrokinetic distribution function h.
 *
 * computes
 * \f$A_{||}(\mathbf{k})\f$
 * potential from
 * \f$h^1_{s0}\f$ (<tt>h</tt> parameter)
 ***************************************/
void fields_getAFromH(const COMPLEX* h){
    size_t flatInd;
    size_t flatInd2D;
    switch(kinetic)
    {
        case ADIABATIC:
            break;
        case NONADIABATIC:
            for(size_t ix = 0; ix < array_local_size.nkx; ix++)
            {
                for(size_t iy = 0; iy < array_local_size.nky; iy++)
                {
                    for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                    {
                        fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] = 0.;
                        for(size_t is = 0; is < array_local_size.ns; is++)
                        {
                            flatInd = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                      iy * array_local_size.nkz * array_local_size.ns + iz * array_local_size.ns + is;
                            fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] += var_var.q[is] * var_var.n[is] *
                                                                                 var_var.vT[is] * var_J0[var_getJIndex(ix,iy,is)] *
                                                                                 sqrt(0.5) * h[flatInd];
                        }
                        flatInd2D = ix * array_local_size.nky + iy;
                        if (fabs(space_kPerp2[flatInd2D]>1e-16))
                        {
                            fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] *= 0.5 * var_var.beta / space_kPerp2[flatInd2D];
                        }
                        else
                        {
                            fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] *= 0.0;
                        }

                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Wrong call of get_A function! Aborting...",mpi_my_rank);
            exit(1);
    }
};

/***************************************
 * \fn void fields_getBFromH(const COMPLEX *h0, const COMPLEX *h1)
 * \brief computes B potential
 * \param h0: 4D complex array (kx,ky,kz,s). Zeroth Hermite and Laguerre moment of gyrokinetic distribution function.
 * \param h1: 4D complex array (kx,ky,kz,s). Zeroth Hermite and first Laguerre moment of distribution function.
 *
 * Computes
 * \f$B_{\perp}(\mathbf{k})\f$
 * from
 * \f$h^0_{s0}\f$ (<tt>h0</tt> parameter) and
 * \f$h^1_{s0}\f$ (<tt>h1</tt> parameter).
 ***************************************/
void fields_getBFromH(const COMPLEX *h0, const COMPLEX *h1) {
    size_t ind4D;
    size_t ind2D;
    switch(kinetic)
    {
        case ADIABATIC:
            break;
        case NONADIABATIC:
            for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                for(size_t iy = 0; iy < array_local_size.nky; iy++){
                    ind2D = ix * array_local_size.nky + iy;
                    if (fabs(phiB_denom[ind2D]) > 1e-16){
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                            for(size_t is = 0; is < array_local_size.ns; is++)
                            {
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] += var_var.n[is] * var_var.T[is] / var_var.B0 *
                                                                                     var_J1[var_getJIndex(ix,iy,is)] *
                                                                                     (h0[ind4D] + h1[ind4D]);
                            }
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] *= - var_var.beta * 0.5;
                        }
                    }
                    else{
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                        }
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Wrong call of get_B function! Aborting...",mpi_my_rank);
            exit(1);
    }

};

/***************************************
 * \fn void fields_getPhiFromH(const COMPLEX* h)
 * \brief computes phi potential
 * \param h: 4D complex array (kx,ky,kz,s). Zeroth Hermite and Laguerre moment of gyrokinetic distribution function.
 *
 * Computes
 * \f$\phi(\mathbf{k})\f$
 * from
 * \f$h^0_{s0}\f$ (<tt>h</tt> parameter)
 ***************************************/
void fields_getPhiFromH(const COMPLEX* h){
    size_t ind4D;
    size_t ind2D;
    double q2nT = 0;
    for (size_t is = 0; is < array_local_size.ns; is++)
    {
        q2nT += var_var.q[is] * var_var.q[is] * var_var.n[is] / var_var.T[is];
    }
    switch(kinetic)
    {
        case ADIABATIC:
            switch(systemType)
            {
                case ELECTROSTATIC:
                    break;
                case ELECTROMAGNETIC:
                    break;
            }
            break;
        case NONADIABATIC:
            switch(systemType)
            {
                case ELECTROSTATIC:
                    break;
                case ELECTROMAGNETIC:
                    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                        for(size_t iy = 0; iy < array_local_size.nky; iy++){
                            ind2D = ix * array_local_size.nky + iy;
                            if (fabs(phiB_denom[ind2D]) > 1e-16){
                                for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                                    for(size_t is = 0; is < array_local_size.ns; is++)
                                    {
                                        ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                                iy * array_local_size.nkz * array_local_size.ns +
                                                iz * array_local_size.ns +
                                                is;
                                        fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] += var_var.q[is] * var_var.n[is] *
                                                                                               var_J0[var_getJIndex(ix,iy,is)] *
                                                                                               h[ind4D];
                                    }
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] /= q2nT;
                                }
                            }
                            else{
                                for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                                }
                            }

                        }
                    }
                    break;
            }
            break;
    }

};

/***************************************
 * \fn void fields_getGradX(COMPLEX *out)
 * \brief computes chi gradient in x direction
 * \param out: output complex array of size (kx,ky,kz,s,Nfields).
 *
 * computes gradient in x direction for chi potentials.
 * <tt>Nfields<\tt> can be 1 or 3,
 * and chosen automatically at start of the simulation
 * depending on the simulation type (electrostatic or electromagnetic)
 ***************************************/
void fields_getGradX(COMPLEX *out){
    size_t ind4D;
    size_t indChiBuf;
    switch(systemType){
        case(ELECTROSTATIC):
            for (size_t ix = 0; ix < array_local_size.nkx; ix ++){
                for (size_t iy = 0; iy < array_local_size.nky; iy ++){
                    for (size_t iz = 0; iz < array_local_size.nkz; iz ++){
                        for (size_t is = 0; is < array_local_size.ns; is ++){
                            ind4D = getIndChi(ix,iy,iz,is);
                            indChiBuf = getIndChiBufEL_c(ix,iy,iz,is);
                            out[indChiBuf] = space_iKx[ix] * fields_chi.phi[ind4D];
                        }
                    }
                }
            }
            break;
        case(ELECTROMAGNETIC):
            for (size_t ix = 0; ix < array_local_size.nkx; ix ++){
                for (size_t iy = 0; iy < array_local_size.nky; iy ++){
                    for (size_t iz = 0; iz < array_local_size.nkz; iz ++){
                        for (size_t is = 0; is < array_local_size.ns; is ++){
                            ind4D = getIndChi(ix,iy,iz,is);
                            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is,CHI_PHI);
                            out[indChiBuf] = space_iKx[ix] * fields_chi.phi[ind4D];

                            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is,CHI_A);
                            out[indChiBuf] = space_iKx[ix] * fields_chi.A[ind4D];

                            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is,CHI_B);
                            out[indChiBuf] = space_iKx[ix] * fields_chi.B[ind4D];
                            //printf("kx = %f\n", cabs(space_iKx[ix]));
                        }
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Error computing gradient chi! aborting...",mpi_my_rank);
            exit(1);
    }
}

/***************************************
 * \fn void fields_getGradY(COMPLEX *out)
 * \brief computes chi gradient in y direction
 * \param out: output complex array of size (kx,ky,kz,s,Nfields).
 *
 * computes gradient in y direction for chi potentials.
 * <tt>Nfields<\tt> can be 1 or 3,
 * and chosen automatically at start of the simulation
 * depending on the simulation type (electrostatic or electromagnetic)
 ***************************************/
void fields_getGradY(COMPLEX *out){
    size_t ind4D;
    size_t indChiBuf;
    switch(systemType){
        case(ELECTROSTATIC):
            for (size_t ix = 0; ix < array_local_size.nkx; ix ++){
                for (size_t iy = 0; iy < array_local_size.nky; iy ++){
                    for (size_t iz = 0; iz < array_local_size.nkz; iz ++){
                        for (size_t is = 0; is < array_local_size.ns; is ++){
                            ind4D = getIndChi(ix,iy,iz,is);
                            indChiBuf = getIndChiBufEL_c(ix,iy,iz,is);
                            out[indChiBuf] = space_iKy[iy] * fields_chi.phi[ind4D];
                        }
                    }
                }
            }
            break;
        case(ELECTROMAGNETIC):
            for (size_t ix = 0; ix < array_local_size.nkx; ix ++){
                for (size_t iy = 0; iy < array_local_size.nky; iy ++){
                    for (size_t iz = 0; iz < array_local_size.nkz; iz ++){
                        for (size_t is = 0; is < array_local_size.ns; is ++){
                            ind4D = getIndChi(ix,iy,iz,is);
                            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is,CHI_PHI);
                            out[indChiBuf] = space_iKy[iy] * fields_chi.phi[ind4D];

                            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is,CHI_A);
                            out[indChiBuf] = space_iKy[iy] * fields_chi.A[ind4D];

                            indChiBuf = getIndChiBufEM_c(ix,iy,iz,is,CHI_B);
                            out[indChiBuf] = space_iKy[iy] * fields_chi.B[ind4D];
                            //printf("ky = %f\n", cabs(space_iKy[iy]));
                        }
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Error computing gradient chi! aborting...",mpi_my_rank);
            exit(1);
    }
}