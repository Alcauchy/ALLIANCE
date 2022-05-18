////////////////////////////////////////////////////////////////////////////////
// 09/02/2022 created by Gene Gorbunov
//                                   TEST UTILITIES
//
// test_freeEnergyComputation
// test_mainFunction
// test_fieldComputation
// test_fieldComparison
// test_kSpecComputations
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////

#include "utils_tests.h"
#include "diagnostics.h"

#define PRINT_A 0
#define PRINT_B 0
#define PRINT_PHI 0
#define TOLERANCE 1e-12

/***************************
 *  vortexInit(COMPLEX *g)
 * *************************/
void vortexInit(COMPLEX *g){
    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
    {
        if(global_nkx_index[ix] == 4)
        {
            g[get_flat_c(0,0,0,ix,3,3)] = 1.j - 0.9;
            g[get_flat_c(0,1,0,ix,8,1)] = 10.j + 0.5;
            g[get_flat_c(0,0,1,ix,2,2)] = 2.+1.j;
        }
    }
}

/***************************
 *  test_fieldComputation()
 * *************************/
void test_fieldComputation(){
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    for(size_t i = 0; i < array_local_size.total_comp; i++)
    {
        g[i] = 0;
    }
    vortexInit(g);
    init_conditions(g);
    fields_sendG(g);
    size_t ind4D;
    size_t ind2D;

    switch(systemType)
    {
        case ELECTROSTATIC:
            break;
        case ELECTROMAGNETIC:
            printf("[MPI process %d] LAUNCHING FIELD COMPUTATION TEST\n",mpi_my_rank);
            fields_getPhi(g00,g01);
            fields_getB(g00,g01);
            fields_getA(g10);
            if(PRINT_A)
            {
                for(size_t ix = 0; ix < array_local_size.nkx; ix++)
                {
                    for(size_t iy = 0; iy < array_local_size.nky; iy++)
                    {
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                        {
                            printf("A[%zu,%zu,%zu] = %f + i %f\n",ix,iy,iz,
                                   creal(fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)]),
                                   cimag(fields_fields.A[get_flatIndexComplex3D(ix,iy,iz)]));
                        }
                    }
                }
            }
            if(PRINT_B){
                for(size_t ix = 0; ix < array_local_size.nkx; ix++)
                {
                    for(size_t iy = 0; iy < array_local_size.nky; iy++)
                    {
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                        {
                            printf("B[%zu,%zu,%zu] = %f + i %f\n",ix,iy,iz,
                                   creal(fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)]),
                                   cimag(fields_fields.B[get_flatIndexComplex3D(ix,iy,iz)]));
                        }
                    }
                }
            }
            if(PRINT_PHI){
                for(size_t ix = 0; ix < array_local_size.nkx; ix++)
                {
                    for(size_t iy = 0; iy < array_local_size.nky; iy++)
                    {
                        for(size_t iz = 0; iz < array_local_size.nkz; iz++)
                        {
                            printf("A[%zu,%zu,%zu] = %f + i %f\n",ix,iy,iz,
                                   creal(fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)]),
                                   cimag(fields_fields.phi[get_flatIndexComplex3D(ix,iy,iz)]));
                        }
                    }
                }
            }
            fields_getChiPhi();
            fields_getChiA();
            fields_getChiB();
            fields_getChiPhi();
            hdf_saveFieldA("A.h5");
            hdf_saveFieldB("B.h5");
            hdf_saveFieldPhi("phi.h5");
            hdf_create_file_c("testFields",g);
            break;
        default:
            printf("[MPI process %d] FIELD COMPUTATION TEST FAILED\n",mpi_my_rank);
            exit(1);
    }
    printf("[MPI process %d] FIELD COMPUTATION TEST PASSED\n",mpi_my_rank);
};

/***************************
 *  test_mainFunction()
 * *************************/
void test_mainFunction(){
    COMPLEX* h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX* g = malloc(array_local_size.total_comp * sizeof(*g));
    init_conditions(h);
    fields_sendG(h);
    fields_getFieldsFromH(g00, g10, g01);
    fields_getChi();
    distrib_getG(g, h);
    for(int it = 0; it < solver.Nt; it++)
    {
        solver_makeStep(g, h);
        //solver_updateDt();
        if(parameters.save_diagnostics && it % parameters.iter_diagnostics == 0)
        {
            diag_compute(g, h, it);
        }
        hdf_saveData(h, it);
    }
};

/***************************
 *  test_fieldComparison()
 * *************************/
void test_fieldComparison(){
    COMPLEX* h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX* g = malloc(array_local_size.total_comp * sizeof(*g));
    COMPLEX* phi_h = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.nkz * sizeof(*phi_h));
    COMPLEX* A_h = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.nkz * sizeof(*A_h));
    COMPLEX* B_h = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.nkz * sizeof(*B_h));

    init_conditions(h);
    fields_sendG(h);
    fields_getFieldsFromH(g00, g10, g01);
    fields_getChi();
    for(size_t i = 0; i < array_local_size.nkx * array_local_size.nky * array_local_size.nkz; i++)
    {
        phi_h[i] = fields_fields.phi[i];
        A_h[i] = fields_fields.A[i];
        B_h[i] = fields_fields.B[i];
    }
    distrib_getG(g, h);
    fields_sendG(g);
    fields_getFields(g00, g10, g01);
    for(size_t i = 0; i < array_local_size.nkx * array_local_size.nky * array_local_size.nkz; i++)
    {
        if(fabs(phi_h[i] - fields_fields.phi[i]) > TOLERANCE)
        {
           printf("[MPI process %d] phi[%zu] not equal !!\n",mpi_my_rank, i);
        }

        if(fabs(A_h[i]-fields_fields.A[i]) > TOLERANCE)
        {
            printf("[MPI process %d] A[%zu] not equal !!\n",mpi_my_rank, i);
        }

        if(fabs(B_h[i] - fields_fields.B[i]) > TOLERANCE)
        {
            printf("[MPI process %d] B[%zu] not equal !!\n",mpi_my_rank, i);
        }
    }

}

/***************************
 *  test_kSpecComputations()
 * *************************/
void test_kSpecComputations(){
    COMPLEX* h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX* g = malloc(array_local_size.total_comp * sizeof(*g));
    init_conditions(h);
    fields_sendG(h);
    fields_getFieldsFromH(g00, g10, g01);
    fields_getChi();
    distrib_getG(g, h);
    diag_computeSpectra(g, h,0);
    for(size_t i = 0; i < parameters.k_shells; i++){
        printf("[MPI process %d] spec[%zu] = %f\n", mpi_my_rank,i,diag_kSpec[i]);
    }
    for(size_t i = 0; i < array_local_size.nm; i++){
        printf("[MPI process %d] spec_m[%zu] = %f\n", mpi_my_rank,i,diag_mSpec[i]);
    }
}

/***************************
 *  test_linearRHS()
 * *************************/
 void test_linearRHS(){
    COMPLEX* h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX* g = malloc(array_local_size.total_comp * sizeof(*g));
    COMPLEX* rhs = malloc(array_local_size.total_comp * sizeof(*rhs));
    init_conditions(h);
    fields_sendG(h);
    fields_getFieldsFromH(g00, g10, g01);
    fields_getChi();
    distrib_getG(g, h);
    diag_computeSpectra(g, h,0);
    equation_getRHS(g, h, rhs);
    hdf_create_file_c("h.h5",h);
    hdf_create_file_c("rhs.h5",rhs);
 }

/***************************
*  test_inplaceFFTW_chi()
* *************************/
void test_inplaceFFTW_chi(){
    //COMPLEX* h = malloc(array_local_size.total_comp * sizeof(*h));
    //COMPLEX* g = malloc(array_local_size.total_comp * sizeof(*g));
    //init_conditions(h);
    //fields_sendG(h);
    //fields_getFieldsFromH(g00, g10, g01);
    //fields_getChi();
    //distrib_getG(g, h);
    size_t indChi;
    COMPLEX *chi = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.nkz * array_local_size.ns * 3 * sizeof(*chi));
    double *chi_r = malloc(array_local_size.nkx * array_local_size.nky * (array_local_size.nz + 2) * array_local_size.ns * 3 * sizeof(*chi_r));
    // initialising chi array
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t is = 0; is < array_local_size.ns; is++){
                    for(size_t ifield = 0 ; ifield < 3; ifield ++){
                        if(global_nkx_index[ix] == 1 && iy == 0 && iz == 0){
                            indChi = getIndChiBufEM_c(ix,iy,iz,is, ifield);
                            chi[indChi] = 1.0;
                        }
                        else{
                            indChi = getIndChiBufEM_c(ix,iy,iz,is, ifield);
                            chi[indChi] = 0.0;
                        }
                    }
                }
            }
        }
    }

    fftw_c2r_chi();

    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nz + 2; iz++){
                for(size_t is = 0; is < array_local_size.ns; is++){
                    for(size_t ifield = 0 ; ifield < 3; ifield ++){
                        if (iz == 0 && iy == 0 && ifield == 0){
                            indChi = getIndChiBufEM_r(ix,iy,iz,is, ifield);
                            printf("[MPI process %d] chi[%zu,%zu,%zu,%zu] = %f\n", mpi_my_rank,global_nkx_index[ix], iy, iz, is, ifield, chi_r[indChi]);
                        }
                    }
                }
            }
        }
    }

    fftw_r2c_chi();

    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t is = 0; is < array_local_size.ns; is++){
                    for(size_t ifield = 0 ; ifield < 3; ifield ++){
                        if(global_nkx_index[ix] == 1 && iy == 0 && iz == 0){
                            indChi = getIndChiBufEM_c(ix,iy,iz,is, ifield);
                            printf("[MPI process %d] chi[%zu,%zu,%zu,%zu] = %f\n", mpi_my_rank,global_nkx_index[ix], iy, iz, is, ifield, creal(chi[indChi]));
                        }
                    }
                }
            }
        }
    }

}

/***************************
*  test_nonlinearTerm()
* *************************/
void test_nonlinearTerm(){
    size_t ind6D;
    size_t ind4D;
    COMPLEX *h = calloc(array_local_size.total_comp,sizeof(*h));
    COMPLEX *buf = calloc(array_local_size.total_comp,sizeof(*buf));
    for(size_t ix = 0; ix < array_local_size.nkx; ix ++){
        for(size_t iy = 0; iy < array_local_size.nky; iy ++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz ++){
                for(size_t im = 0; im < array_local_size.nm; im ++){
                    for(size_t il = 0; il < array_local_size.nl; il ++){
                        for(size_t is = 0; is < array_local_size.ns; is ++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            if(global_nkx_index[ix] == 1 && iy == 0 && iz == 0 && im == 0 && is == 0 && il == 0){
                                h[ind6D] = 1.;

                            }
                            else{
                                h[ind6D] = 0;
                            }
                            ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                    iy * array_local_size.nkz * array_local_size.ns +
                                    iz * array_local_size.ns +
                                    is;
                            if(global_nkx_index[ix] == 0 && iy == 1 && iz == 0){
                                fields_chi.phi[ind4D] = 1.;
                                fields_chi.A[ind4D] = 0.;
                                fields_chi.B[ind4D] = 0.;
                            }
                            else{
                                fields_chi.phi[ind4D] = 0.;
                                fields_chi.A[ind4D] = 0.;
                                fields_chi.B[ind4D] = 0.;
                            }
                        }
                    }
                }
            }
        }
    }
    fftw_copy_buffer_c(fftw_hBuf, h);
    fftw_c2r();
    fftw_r2c();
    ind6D = get_flat_c(0,0,0,1,0,0);
    if(mpi_my_rank == 0) printf("%f\n",cabs(fftw_hBuf[ind6D]));
    equation_getNonlinearTerm(h, buf);
    for(size_t ix = 0; ix < array_local_size.nkx; ix ++){
        for(size_t iy = 0; iy < array_local_size.nky; iy ++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                ind6D = get_flat_c(0,0,0,ix,iy,iz);
                printf("[MPI process %d] %f\n", cabs(space_iKy[iy]));
                printf("[MPI process %d] buf[%d,%d,%d] = %.8f\n",mpi_my_rank,ix,iy,iz,cabs(buf[ind6D]));
            }
        }
    }


}

/***************************
*  test_Poisson()
* *************************/
void test_Poisson(){
    COMPLEX *h = calloc(array_local_size.total_comp, sizeof(*h));
    COMPLEX *g = calloc(array_local_size.total_comp, sizeof(*g));
    size_t total_size = array_local_size.nkx * array_local_size.nky * (array_local_size.nz + 2) * array_local_size.nm * array_local_size.nl * array_local_size.ns;
    double *h_r = calloc(total_size, sizeof(*h_r));
    double *g_r= calloc(total_size, sizeof(*g_r));
    double *buff_r= calloc(total_size, sizeof(*buff_r));
    printf( "%zu, %zu\n", total_size, array_local_size.total_real);
    init_conditions(h);
    init_conditions(g);
    //for (size_t ii = 0; ii< array_local_size.total_comp; ii++) g[ii] = 1;
    size_t ind6D;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t il = 0; il < array_local_size.nl; il++){
                    for(size_t im = 0; im < array_local_size.nm; im++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            //if(global_nkx_index[ix] == 1 && iz != 0 ) h[ind6D] = (1.+1.j) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                            //if(global_nkx_index[ix] == 4 && iz != 0 ) g[ind6D] = (2.+1.j) * array_local_size.nkx * array_local_size.nky * (array_local_size.nz);

                        }
                    }
                }
            }
        }
    }
    for (size_t ii = 0; ii< array_local_size.total_comp; ii++) printf("%f, %f, %f \n",cabs(h[ii]),cabs(g[ii]),cabs(g[ii]-h[ii]));
    //dealiasing23(h);
    //dealiasing23(g);
    // get dhdx and transform it
    distrib_getXGrad(h, fftw_hBuf);
    fftw_c2r();
    double *hBuf = fftw_hBuf;
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = hBuf[ii];

    // get dgdy and transform it
    distrib_getYGrad(g, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) g_r[ii] = hBuf[ii];

    // get product dhdx * dgdy
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] += h_r[ii] * g_r[ii];

    // dhdy
    distrib_getYGrad(h, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = hBuf[ii];

    //dgdx
    distrib_getXGrad(g, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) g_r[ii] = hBuf[ii];

    // dhdx * dgdy - dhdy*dgdx
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] -= h_r[ii] * g_r[ii];

    //dealiasing
    for(size_t ii = 0; ii < total_size; ii++) hBuf[ii] = buff_r[ii];
    fftw_r2c();
    //dealiasing23(fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] = hBuf[ii];

    // transform h to real space
    for(size_t ii = 0; ii < array_local_size.total_comp; ii++) fftw_hBuf[ii] = g[ii];
    fftw_c2r();
    double total = 0;
    for(size_t ii = 0; ii < total_size; ii++) total += buff_r[ii];//* hBuf[ii]

    double total_reduced = 0;
    MPI_Allreduce(&total,
                  &total_reduced,
                      1,
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);

    printf("[MPI process %d] POISSON BRACKET = %16.5e\n",mpi_my_rank, total_reduced);

    //
    // computing flux as {g, h^2}
    //
    for  (size_t ii = 0 ; ii < total_size; ii++) buff_r[ii] = 0;
    // compute h^2 first, then we can continue as usual
    fftw_copy_buffer_c(fftw_hBuf, h);
    fftw_c2r();
    fftw_copy_buffer_r(h_r, fftw_hBuf);
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = h_r[ii] * h_r[ii];
    fftw_copy_buffer_r(fftw_hBuf, h_r);
    fftw_r2c();
    dealiasing23(fftw_hBuf);
    fftw_copy_buffer_c(h, fftw_hBuf);

    // get dhdx and transform it
    distrib_getXGrad(h, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = hBuf[ii];

    // get dgdy and transform it
    distrib_getYGrad(g, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) g_r[ii] = hBuf[ii];

    // get product dhdx * dgdy
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] += h_r[ii] * g_r[ii];

    // dhdy
    distrib_getYGrad(h, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = hBuf[ii];

    //dgdx
    distrib_getXGrad(g, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) g_r[ii] = hBuf[ii];

    // dhdx * dgdy - dhdy*dgdx
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] -= h_r[ii] * g_r[ii];

    //dealiasing
    for(size_t ii = 0; ii < total_size; ii++) hBuf[ii] = buff_r[ii];
    fftw_r2c();
    dealiasing23(fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] = hBuf[ii];

    // transform h to real space
    for(size_t ii = 0; ii < array_local_size.total_comp; ii++) fftw_hBuf[ii] = g[ii];
    fftw_c2r();
    double total_1 = 0;
    for(size_t ii = 0; ii < total_size; ii++) total_1 += buff_r[ii];
    MPI_Allreduce(&total_1,
                  &total_reduced,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);
    printf("[MPI process %d] POISSON BRACKET 1  = %16.5e\n",mpi_my_rank, total_reduced);




}

/***************************
*  test_xGrad()
* *************************/
void test_xGrad(){
    size_t w = 2;
    COMPLEX *h = calloc(array_local_size.total_comp, sizeof(*h));
    COMPLEX *dhdx = calloc(array_local_size.total_comp, sizeof(*dhdx));
    size_t ind6D;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t il = 0; il < array_local_size.nl; il++){
                    for(size_t im = 0; im < array_local_size.nm; im++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            if(global_nkx_index[ix] == w && iz == 0 && iy == 0) h[ind6D] = (0.5) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                            if(global_nkx_index[ix] == array_local_size.nkx - w && iz == 0 && iy == 0) h[ind6D] = (0.5) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                        }
                    }
                }
            }
        }
    }
    double *real_buf = fftw_hBuf;
    fftw_copy_buffer_c(fftw_hBuf, h);
    fftw_c2r();
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        ind6D = get_flat_r(0,0,0,ix,0,0);
        printf("h[%zu] = %f\n",ix, real_buf[ind6D]);
    }
    distrib_getXGrad(h,dhdx);
    fftw_copy_buffer_c(fftw_hBuf, dhdx);
    fftw_c2r();
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        ind6D = get_flat_r(0,0,0,ix,0,0);
        printf("dhdx[%zu] = %f\n",ix, real_buf[ind6D]);
    }
}

/***************************
*  test_Poisson1()
* *************************/
void test_Poisson1(){
    size_t wh_x = 2;
    size_t wh_y = 3;
    size_t wg_x = 5;
    size_t wg_y = 2;
    COMPLEX *h = calloc(array_local_size.total_comp, sizeof(*h));
    COMPLEX *g = calloc(array_local_size.total_comp, sizeof(*g));
    double *h_r = calloc(array_local_size.total_real, sizeof(*h_r));
    double *g_r= calloc(array_local_size.total_real, sizeof(*g_r));
    double *buff_r= calloc(array_local_size.total_real, sizeof(*buff_r));
    size_t ind6D;
    // initialize arrays with some cosine
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t il = 0; il < array_local_size.nl; il++){
                    for(size_t im = 0; im < array_local_size.nm; im++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            if(global_nkx_index[ix] == wh_x && iz == 0 && iy == wh_y){
                                h[ind6D] = (0.5) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                            }
                            if(global_nkx_index[ix] == array_local_size.nkx - wh_x && iz == 0 && iy == array_local_size.nky - wh_y){
                                h[ind6D] = (0.5) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                            }
                            if(global_nkx_index[ix] == wg_x && iz == 0 && iy == wg_y){
                                g[ind6D] = (0.5) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                            }
                            if(global_nkx_index[ix] == array_local_size.nkx - wg_x && iz == 0 && iy == array_local_size.nky - wg_y){
                                g[ind6D] = (0.5) *  array_local_size.nkx * array_local_size.nky * (array_local_size.nz);
                            }
                        }
                    }
                }
            }
        }
    }
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                ind6D = get_flat_r(0,0,0,ix,iy,iz);
                h_r[ind6D] = cos(2. * M_PI * wh_x / array_global_size.nkx * ix ) * sin(2. * M_PI * wh_y / array_global_size.nky * iy );
                g_r[ind6D] = sin(2. * M_PI * wg_x / array_global_size.nkx * ix ) * sin(2. * M_PI * wg_y / array_global_size.nky * iy );
            }
        }
    }
    // fft(h_r)
    fftw_copy_buffer_r(fftw_hBuf,h_r);
    fftw_r2c();
    fftw_copy_buffer_c(h,fftw_hBuf);
    // fft(g_r)
    fftw_copy_buffer_r(fftw_hBuf,g_r);
    fftw_r2c();
    fftw_copy_buffer_c(g,fftw_hBuf);
    //dhdx
    distrib_getXGrad(h,fftw_hBuf);
    fftw_c2r();
    fftw_copy_buffer_r(h_r,fftw_hBuf);
    //dgdy
    distrib_getYGrad(g,fftw_hBuf);
    fftw_c2r();
    fftw_copy_buffer_r(g_r,fftw_hBuf);
    // dhdx * dgdy
    for(size_t ii = 0; ii<array_local_size.total_real; ii++) buff_r[ii] += h_r[ii] * g_r[ii];
    //dhdy
    distrib_getYGrad(h,fftw_hBuf);
    fftw_c2r();
    fftw_copy_buffer_r(h_r,fftw_hBuf);
    //dgdx
    distrib_getXGrad(g,fftw_hBuf);
    fftw_c2r();
    fftw_copy_buffer_r(g_r,fftw_hBuf);
    // dhdxy *dgdx
    for(size_t ii = 0; ii<array_local_size.total_real; ii++) buff_r[ii] -= h_r[ii] * g_r[ii];
    double total = 0;
    for(size_t ii = 0; ii <array_local_size.total_real; ii++) total += buff_r[ii];
    printf("POISSON BRACKET = %16.5e\n", total);
    free(h);
    free(h_r);
    free(g);
    free(g_r);
    free(buff_r);
}

/***************************
*  test_enforceRealityConditions()
* *************************/
void test_enforceRealityConditions(){
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    init_conditions(h);
    distrib_enforceReality(h);
    hdf_create_file_c("test_real_cond.h5",h);
}
/***************************
*  test_enforceRealityConditions()
* *************************/
void test_enforceZero(){
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    init_conditions(h);
    distrib_setZeroNHalf(h);
    hdf_create_file_c("test_enforceZero.h5",h);
}