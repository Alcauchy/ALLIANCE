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
        solver_makeStep(g);
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