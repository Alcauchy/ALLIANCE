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
            printf("[MPI process %d] B[%zu] not equal !! B_h = (%f,%f), B_g = (%f, %f)\n",mpi_my_rank, i,creal(B_h[i]),cimag(B_h[i]), creal(fields_fields.B[i]),cimag(fields_fields.B[i]));
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
/*void test_nonlinearTerm(){
    size_t ind6D;
    size_t ind4D;
    COMPLEX *h = calloc(array_local_size.total_comp,sizeof(*h));
    double *buf = calloc(array_local_size.total_real,sizeof(*buf));
    init_conditions(h);
    //random initialization of chi function
    for(size_t ii = 0; ii < array_local_size.nkx * array_local_size.nky * array_local_size.nkz * array_local_size.ns; ii++){
        fields_chi.A[ii] =0;//((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
        fields_chi.B[ii] =0;//((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
        fields_chi.phi[ii] = ((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
    }
    //equation_getNonlinearTerm(h,buf);
    //
    // dhdx and dchidy
    //
    distrib_getXGrad(h,fftw_hBuf);
    fields_getGradY(fftw_chiBuf);
    fftw_c2r();
    fftw_c2r_chi();
    double *h_r = fftw_hBuf;
    double *chi_r = fftw_chiBuf;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nz + 2; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_r(is,il,im,ix,iy,iz);
                            ind4D = getIndChiBufEM_r(ix,iy,iz,is,0);
                            buf[ind6D] += h_r[ind6D] * chi_r[ind4D];
                        }
                    }
                }
            }
        }
    }
    //
    // dhdy and dchidx
    //
    distrib_getYGrad(h,fftw_hBuf);
    fields_getGradX(fftw_chiBuf);
    fftw_c2r();
    fftw_c2r_chi();
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nz + 2; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_r(is,il,im,ix,iy,iz);
                            ind4D = getIndChiBufEM_r(ix,iy,iz,is,0);
                            buf[ind6D] -= h_r[ind6D] * chi_r[ind4D];
                        }
                    }
                }
            }
        }
    }
    double buf_total = 0;
    for (size_t ii = 0; ii < array_local_size.total_real; ii ++) buf_total += buf[ii];
    printf("total = %f\n", buf_total);

}*/

void test_nonlinearTerm(){
    size_t ind6D;
    size_t ind4D;
    COMPLEX *h = calloc(array_local_size.total_comp,sizeof(*h));
    COMPLEX *buf = calloc(array_local_size.total_comp,sizeof(*buf));
    init_conditions(h);
    //random initialization of chi function
    for(size_t ii = 0; ii < array_local_size.nkx * array_local_size.nky * array_local_size.nkz * array_local_size.ns; ii++){
        fields_chi.A[ii] =((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
        fields_chi.B[ii] =((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
        fields_chi.phi[ii] = ((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
    }
    equation_getNonlinearTerm(h,buf);
}


/***************************
*  test_Poisson()
* *************************/
void test_Poisson(){
    COMPLEX *h = calloc(array_local_size.total_comp, sizeof(*h));
    COMPLEX *g = calloc(array_local_size.total_comp, sizeof(*g));
    size_t total_size = array_local_size.nkx * array_local_size.nky * (array_local_size.nz + 2) * array_local_size.nm * array_local_size.nl * array_local_size.ns;
    double *h_r = calloc(total_size, sizeof(*h_r));
    double *h_r1 = calloc(total_size, sizeof(*h_r1));
    double *g_r= calloc(total_size, sizeof(*g_r));
    double *buff_r= calloc(total_size, sizeof(*buff_r));
   // printf( "%zu, %zu\n", total_size, array_local_size.total_real);
    int my_proc = mpi_whereIsX[2*2];
    int my_local_x = mpi_whereIsX[2*2 + 1];
    if (mpi_my_row_rank == my_proc){
        size_t ind6D = get_flat_c(0,0,0,my_local_x,0,0);
        h[ind6D] = array_global_size.nkx * array_global_size.nky * array_global_size.nz;
    }
    my_proc = mpi_whereIsX[2*(array_global_size.nkx - 2)];
    my_local_x = mpi_whereIsX[2*(array_global_size.nkx - 2) + 1];
    if (mpi_my_row_rank == my_proc){
        size_t ind6D = get_flat_c(0,0,0,my_local_x,0,0);
        h[ind6D] = array_global_size.nkx * array_global_size.nky * array_global_size.nz;
    }
    //init_conditions(h);
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
    //for (size_t ii = 0; ii< array_local_size.total_comp; ii++) printf("%f, %f, %f \n",cabs(h[ii]),cabs(g[ii]),cabs(g[ii]-h[ii]));
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

    ///////////////////////////////
    // computing flux as {g, h^2}
    //////////////////////////////
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
    for (size_t ii = 0; ii < array_local_size.total_comp; ii++){h[ii] = 0;}
    if (mpi_my_row_rank == my_proc){
        size_t ind6D = get_flat_c(0,0,0,my_local_x,0,0);
        h[ind6D] = array_global_size.nkx * array_global_size.nky * array_global_size.nz;
    }
    my_proc = mpi_whereIsX[2*(array_global_size.nkx - 2)];
    my_local_x = mpi_whereIsX[2*(array_global_size.nkx - 2) + 1];
    if (mpi_my_row_rank == my_proc){
        size_t ind6D = get_flat_c(0,0,0,my_local_x,0,0);
        h[ind6D] = array_global_size.nkx * array_global_size.nky * array_global_size.nz;
    }
    init_conditions(h);
    ///////////////////////////////
    // computing flux as {g, h}h
    //////////////////////////////
    for  (size_t ii = 0 ; ii < total_size; ii++) buff_r[ii] = 0;
    // compute h first, then we can continue as usual
    fftw_copy_buffer_c(fftw_hBuf, h);
    fftw_c2r();
    fftw_copy_buffer_r(h_r1, fftw_hBuf);
    // get dhdx and transform it
    distrib_getXGrad(h, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = hBuf[ii];

    // get dgdy and transform it
    distrib_getYGrad(g, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) g_r[ii] = hBuf[ii];

    // get product dhdx * dgdy
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] += h_r[ii] * g_r[ii] ;

    // dhdy
    distrib_getYGrad(h, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) h_r[ii] = hBuf[ii];

    //dgdx
    distrib_getXGrad(g, fftw_hBuf);
    fftw_c2r();
    for(size_t ii = 0; ii < total_size; ii++) g_r[ii] = hBuf[ii];

    // dhdx * dgdy - dhdy*dgdx
    for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] -= h_r[ii] * g_r[ii] ;

    //dealiasing
    //for(size_t ii = 0; ii < total_size; ii++) hBuf[ii] = buff_r[ii];
    //fftw_r2c();
    //dealiasing23(fftw_hBuf);
    //fftw_c2r();
    //for(size_t ii = 0; ii < total_size; ii++) buff_r[ii] = hBuf[ii];

    for(size_t ii = 0; ii < total_size; ii++) total_1 += buff_r[ii] * h_r1[ii];
    MPI_Allreduce(&total_1,
                  &total_reduced,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);
    printf("[MPI process %d] POISSON BRACKET 1  = %16.5e\n",mpi_my_rank, total_reduced);
    hdf_create_file_r("buf_r.h5",buff_r);
    hdf_create_file_r("h_r1.h5",h_r1);



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
    hdf_create_file_c("test_real_cond_0.h5",h);
    distrib_enforceReality(h);
    hdf_create_file_c("test_real_cond_1.h5",h);
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

/***************************
*  test_fieldsFFT()
* *************************/
void test_fieldsFFT(){
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    init_conditions(h);
    fftw_copy_buffer_c(fftw_hBuf,h);
    fftw_c2r();
    hdf_create_file_r("h_field.h5",fftw_hBuf);
    fields_sendG(h);
    fields_getFieldsFromH(g00,g10,g01);
    fftw_copyFieldBuf_c(fftw_field,fields_fields.phi);
    double *field_r = fftw_field;
    fftw_c2r_field();
    for(size_t ix = 0; ix < array_local_size.nkx; ix ++){
        for(size_t iy = 0; iy < array_local_size.nky; iy ++){
            for(size_t iz = array_local_size.nz; iz < array_local_size.nz + 2; iz ++){
                size_t ind3D = ix * array_local_size.nky * (array_local_size.nz + 2) +
                               iy * (array_local_size.nz + 2) + iz;
                printf("%f\n",field_r[ind3D]);
            }
        }
    }
    hdf_saveField_r(fftw_field,"phi.h5");

    fftw_copyFieldBuf_c(fftw_field,fields_fields.A);
    fftw_c2r_field();
    hdf_saveField_r(fftw_field,"A.h5");

    fftw_copyFieldBuf_c(fftw_field,fields_fields.B);
    fftw_c2r_field();
    hdf_saveField_r(fftw_field,"B.h5");
}

/***************************
*  test_everything()
* *************************/
void test_everything(){
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    //initializing g
    init_conditions(h);
    fields_sendG(h);
    //computing fields
    fields_getFieldsFromH(g00,g10,g01);
    fields_getChi();
    // computing h function
    distrib_getG(g, h);
    //computing diagnostics
    diag_compute(g, h, 0);
    //saving distribution functions, both real and complex space
    char name[64];
    sprintf(name, "%s%s%s%d%s", ".","/","g_",0,".h5");
    hdf_create_file_c(name,g);
    sprintf(name, "%s%s%s%d%s", ".","/","h_",0,".h5");
    hdf_create_file_c(name,h);
    //fftw transform to save real distribution functions
    fftw_copy_buffer_c(fftw_hBuf, g);
    fftw_c2r();
    sprintf(name, "%s%s%s%d%s", ".","/","gr_",0,".h5");
    hdf_create_file_r(name,(double *) fftw_hBuf);

    fftw_copy_buffer_c(fftw_hBuf, h);
    fftw_c2r();
    sprintf(name, "%s%s%s%d%s", ".","/","hr_",0,".h5");
    hdf_create_file_r(name,(double *) fftw_hBuf);
    //fftw transform of the fields
    fftw_copyFieldBuf_c(fftw_field,fields_fields.phi);
    fftw_c2r_field();
    hdf_saveField_r(fftw_field,"phi.h5");

    fftw_copyFieldBuf_c(fftw_field,fields_fields.A);
    fftw_c2r_field();
    hdf_saveField_r(fftw_field,"A.h5");

    fftw_copyFieldBuf_c(fftw_field,fields_fields.B);
    fftw_c2r_field();
    hdf_saveField_r(fftw_field,"B.h5");

    hdf_saveData(g,0);
    free(g);
    free(h);
}

/***************************
*  test_RHS()
* *************************/
void test_RHS(){
    //
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *out = malloc(array_local_size.total_comp * sizeof(*out));
    //initializing g
    init_conditions(h);
    fields_sendG(h);
    //computing fields
    fields_getFieldsFromH(g00,g10,g01);
    /*for (size_t ii = 0; ii < array_local_size.nkx * array_local_size.nky * array_local_size.nkz; ii++){
        fields_fields.A[ii] = 0;
    }*/
    fields_getChi();
    // computing h function
    distrib_getG(g, h);
    //computing diagnostics
    diag_compute(g, h, 0);
    //compute dh/dx
    distrib_getXGrad(h,fftw_hBuf);
    fftw_c2r();
    char name[64];
    sprintf(name, "%s%s%s%d%s", ".","/","dhdx_",0,".h5");
    hdf_create_file_r(name,fftw_hBuf);
    //compute and save dh/dy
    distrib_getYGrad(h,fftw_hBuf);
    fftw_c2r();
    sprintf(name, "%s%s%s%d%s", ".","/","dhdy_",0,".h5");
    hdf_create_file_r(name,fftw_hBuf);
    //compute and save dchidx
    fields_getGradX(fftw_chiBuf);
    fftw_c2r_chi();
    sprintf(name, "%s%s%s%d%s", ".","/","dchidx_",0,".h5");
    hdf_createChiFile_r(name, fftw_chiBuf);
    //compute and save dchidy
    fields_getGradY(fftw_chiBuf);
    fftw_c2r_chi();
    sprintf(name, "%s%s%s%d%s", ".","/","dchidy_",0,".h5");
    hdf_createChiFile_r(name, fftw_chiBuf);

    //compute rhs and save it
    equation_getNonlinearTerm(h,out);
    fftw_copy_buffer_c(fftw_hBuf,out);
    fftw_c2r();
    sprintf(name, "%s%s%s%d%s", ".","/","rhs_",0,".h5");
    //size_t ind6D = get_flat_r(0,0,15,0,0,0);
    //printf(ind6D = get_flat_r(0,0,,0,0,0));

    hdf_create_file_r(name, fftw_hBuf);
    //free memory
    free(out);
    free(g);
    free(h);
}

/***************************
*  test_CFL()
* *************************/
void test_CFL(){
    //
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *out = malloc(array_local_size.total_comp * sizeof(*out));
    //initializing g
    init_conditions(h);
    fields_sendG(h);
    //computing fields
    fields_getFieldsFromH(g00,g10,g01);
    /*for (size_t ii = 0; ii < array_local_size.nkx * array_local_size.nky * array_local_size.nkz; ii++){
        fields_fields.A[ii] = 0;
    }*/
    fields_getChi();
    // computing h function
    distrib_getG(g, h);
    double C_max = 0.5;
    //computing cfl from complex data
    //kmax = nx / L * 2pi
    double kxMax = 2. * M_PI / 100. * array_global_size.nkx;
    double kyMax = 2. * M_PI / 100. * array_global_size.nky;
    //double *dt_ar =malloc();
    free(out);
    free(g);
    free(h);
}

/***************************
*  test_transposedFFTW()
* *************************/
void test_transposedFFTW(){
    COMPLEX *h = calloc(array_local_size.total_comp , sizeof(*h));
    int kx_local,kx_localNeg;
    int kx = 0;
    int kxNeg = array_global_size.nkx - kx;
    int proc, procNeg;
    proc = mpi_whereIsX[2 * kx];
    kx_local = mpi_whereIsX[2 * kx + 1];
    procNeg = mpi_whereIsX[2*kxNeg];
    kx_localNeg = mpi_whereIsX[2*kxNeg + 1];
    if(mpi_my_row_rank == proc){
        size_t iy = 2;
        size_t ind6D = get_flat_c(0,0,0,kx_local,iy,0);
        h[ind6D] = array_global_size.nkx*array_global_size.nky*array_global_size.nz/2.;
        ind6D = get_flat_c(0,0,0,kx_local,array_global_size.nky - iy,0);
        h[ind6D] = array_global_size.nkx*array_global_size.nky*array_global_size.nz/2.;

        size_t indChi = getIndChiBufEM_c(kx_local,iy,0,0,0);
        fftw_chiBuf[indChi] =  array_global_size.nkx*array_global_size.nky*array_global_size.nz/2.;
        indChi = getIndChiBufEM_c(kx_local,array_global_size.nky - iy,0,0,0);
        fftw_chiBuf[indChi] =  array_global_size.nkx*array_global_size.nky*array_global_size.nz/2.;

        size_t ind3D = get_flatIndexComplex3D(kx_local,iy,0);//kx_local * array_local_size.nky * array_local_size.nkz + iy * array_local_size.nkz;
        fftw_field[ind3D] = (array_global_size.nkx*array_global_size.nky*array_global_size.nz)/2.;
        ind3D = kx_local * array_local_size.nky * array_local_size.nkz + (array_global_size.nky - iy) * array_local_size.nkz;
        fftw_field[ind3D] = array_global_size.nx*array_global_size.ny*array_global_size.nz/2.;
    }
    double *hr = h;
    fftw_copy_buffer_c(fftw_hBuf,h);
    fftw_c2r();
    fftw_copy_buffer_r(h,fftw_hBuf);
    double *hbuf = calloc(array_local_size.total_real , sizeof(*hbuf));
    fftw_copy_buffer_r(hbuf,fftw_hBuf);
    if(mpi_my_coords[0] == 0){
        for (size_t iy = 0; iy < array_local_size.ny;iy++){
            size_t ind6D = get_flat_r(0,0,0,0,iy,0);
            //printf("[MPI process %d] f[%zu] = %f\n",mpi_my_rank,iy,hr[ind6D]);
            //printf("[MPI process %d] b[%zu] = %f\n",mpi_my_rank,ind6D,hbuf[ind6D]);
        }
    }
    hdf_create_file_r("h.h5", hr);
    fftw_transposeToXY();
    fftw_transposeToYX();
    if(mpi_my_coords[0] == 0){
        for (size_t iy = 0; iy < array_local_size.ny;iy++){
            size_t ind6D = get_flat_r(0,0,0,0,iy,0);
           // printf("[MPI process %d] f[%zu] = %f\n",mpi_my_rank,iy,hr[ind6D]);
        }
    }
    fftw_c2r_chi();
    double *chi_r = fftw_chiBuf;
    if(mpi_my_coords[0] == 0){
        for (size_t iy = 0; iy < array_local_size.ny;iy++){
            size_t indChi = getIndChiBufEM_r(0,iy,0,0,0);
            printf("[MPI process %d] chi[%zu] = %f\n",mpi_my_rank,iy,chi_r[indChi]);
            //printf("[MPI process %d] b[%zu] = %f\n",mpi_my_rank,ind6D,hbuf[ind6D]);
        }
    }
    fftw_transposeToXY_chi();
    fftw_transposeToYX_chi();
    if(mpi_my_coords[0] == 0){
        for (size_t iy = 0; iy < array_local_size.ny;iy++){
            size_t indChi = getIndChiBufEM_r(0,iy,0,0,0);
            printf("[MPI process %d] chi[%zu] = %f\n",mpi_my_rank,iy,chi_r[indChi]);
            //printf("[MPI process %d] b[%zu] = %f\n",mpi_my_rank,ind6D,hbuf[ind6D]);
        }
    }
    fftw_r2c_chi();

    // fftw transposed of fields
    fftw_c2r_field();
    double *f_r = fftw_field;
    if(mpi_my_coords[0] == 0){
        for (size_t iy = 0; iy < array_local_size.ny;iy++){
            size_t ind3D = iy * array_local_size.nx * (array_local_size.nz + 2);
            printf("[MPI process %d] field[%zu] = %f\n",mpi_my_rank,iy,f_r[ind3D]);
            //printf("[MPI process %d] b[%zu] = %f\n",mpi_my_rank,ind6D,hbuf[ind6D]);
        }
    }
    fftw_transposeToXY_field();
    fftw_transposeToYX_field();
    fftw_r2c_field();
}