//
// Created by alcauchy on 26/01/2022.
//
#include "fields.h"

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
COMPLEX *g00;
COMPLEX *g10;
COMPLEX *g01;

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

            g00 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*g00));
            g10 = 0;
            g01 = 0;
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

            g00 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*g00));
            g10 = malloc(array_local_size.nkx *
                         array_local_size.nky *
                         array_local_size.nkz *
                         array_local_size.ns *
                         sizeof(*g10));
            g01 = malloc(array_local_size.nkx *
                               array_local_size.nky *
                               array_local_size.nkz *
                               array_local_size.ns *
                               sizeof(*g01));

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
                    A_denom[flatInd] = 1.0/(space_kPerp2[flatInd] + sum[flatInd]);

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
                }
            }
            break;

        default:
            printf("[MPI process %d] mistake defining fields! Aborting...",mpi_my_rank);
            exit(1);
    }
};

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
                        fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)] *= A_denom[flatInd2D];
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Wrong call of get_A function! Aborting...",mpi_my_rank);
            exit(1);
    }
};

void fields_getB(const COMPLEX* g0, const COMPLEX* g1) {
    size_t ind2D;
    size_t ind3D;
    size_t ind4D;
    switch(kinetic)
    {
        case ADIABATIC:
            break;
        case NONADIABATIC:
            for(size_t ix = 0; ix < array_local_size.nkx; ix++)
            {
                for(size_t iy = 0; iy < array_local_size.nky; iy++)
                {
                    ind2D = ix * array_local_size.nky + iy;

                    if (fabs(phiB_denom[ind2D]) > 1e-16)
                    {
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                        {
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                            for(size_t is = 0; is < array_local_size.ns; is++)
                            {
                                ind3D = ix * array_local_size.nky  * array_local_size.ns +
                                        iy  * array_local_size.ns +
                                        is;
                                ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                        iy * array_local_size.nkz * array_local_size.ns +
                                        iz * array_local_size.ns +
                                        is;
                                fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] += a_pot[ind2D] * I_B[ind3D] *
                                                                                     (g0[ind4D] + g1[ind4D]);
                                fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] += c_pot[ind2D] * I_phi[ind3D] *
                                                                                     g0[ind4D];

                            }
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] *= - var_var.beta / phiB_denom[ind2D];
                        }

                    }
                    else
                    {
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                        {
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
                                    fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] *=  1.0 / phiB_denom[ind2D];
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

void fields_getFields(COMPLEX *g00, COMPLEX *g10, COMPLEX *g01) {
    fields_getPhi(g00,g10);
    fields_getB(g00,g01);
    fields_getA(g10);
};

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

void fields_getChiB(){
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
                    fields_chi.B[ind4D] = var_var.T[is] / (var_var.q[is] * var_var.B0) *
                                          var_J1[var_getJIndex(ix,iy,is)] *
                                          fields_fields.B[ind3D];
                }
            }
        }
    }
}

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

void fields_sendG(COMPLEX *g){
    /* fill buffers with data from processor which has required data */
    size_t ind4D;
    int count;
    int broadcast_root0 = 0;
    int broadcast_root1 = 0;
    switch(systemType)
    {
        case ELECTROSTATIC:
            break;
        case ELECTROMAGNETIC:
            for(size_t im = 0; im < array_local_size.nm; im++)
            {
                if (global_nm_index[im] == 0)
                {
                    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
                    {
                        for(size_t iy = 0; iy < array_local_size.nky; iy++)
                        {
                            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                            {
                                for(size_t is = 0; is < array_local_size.ns; is++)
                                {
                                    ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                            iy * array_local_size.nkz * array_local_size.ns +
                                            iz * array_local_size.ns +
                                            is;
                                    g00[ind4D] = g[get_flat_c(is,0,im,ix,iy,iz)];
                                    g01[ind4D] = g[get_flat_c(is,1,im,ix,iy,iz)];
                                }
                            }
                        }
                    }

                }
                if (global_nm_index[im] == 1)
                {
                    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
                    {
                        for(size_t iy = 0; iy < array_local_size.nky; iy++)
                        {
                            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                            {
                                for(size_t is = 0; is < array_local_size.ns; is++)
                                {
                                    ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                            iy * array_local_size.nkz * array_local_size.ns +
                                            iz * array_local_size.ns +
                                            is;
                                    g10[ind4D] = g[get_flat_c(is,0,im,ix,iy,iz)];
                                }
                            }
                        }
                    }
                }
            }
            count = array_local_size.nkx * array_local_size.nky * array_local_size.nkz * array_local_size.ns;
            MPI_Bcast(g00,count,MPI_DOUBLE_COMPLEX,broadcast_root0,mpi_col_comm);
            MPI_Bcast(g01,count,MPI_DOUBLE_COMPLEX,broadcast_root0,mpi_col_comm);
            MPI_Bcast(g10,count,MPI_DOUBLE_COMPLEX,broadcast_root0,mpi_col_comm);

            break;
    }
};

void fields_getFieldsFromH(COMPLEX *h00, COMPLEX *h10, COMPLEX *h01){
    fields_getPhiFromH(h00);
    fields_getBFromH(h00, h01);
    fields_getAFromH(h10);
};

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

void fields_getBFromH(const COMPLEX *h0, const COMPLEX *h1) {
    size_t ind4D;
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
                        fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] = 0.j;
                        for(size_t is = 0; is < array_local_size.ns; is++)
                        {
                            ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                    iy * array_local_size.nkz * array_local_size.ns +
                                    iz * array_local_size.ns +
                                    is;
                            fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] += var_var.n[is] * var_var.T[is] / var_var.B0 *
                                                                                     var_J0[var_getJIndex(ix,iy,is)] *
                                                                                     (h0[ind4D] + h1[ind4D]);
                        }
                        fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)] *= - var_var.beta * 0.5;
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] Wrong call of get_B function! Aborting...",mpi_my_rank);
            exit(1);
    }

};

void fields_getPhiFromH(const COMPLEX* h){
    size_t ind4D;
    double q2nT = 0;
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
                            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                            {
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
                                    q2nT += var_var.q[is] * var_var.q[is] * var_var.n[is] / var_var.T[is];
                                }
                                fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)] /= q2nT;
                            }
                        }
                    }
                    break;
            }
            break;
    }

};

