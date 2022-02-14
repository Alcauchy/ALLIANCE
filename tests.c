//
// Created by alcauchy on 09/02/2022.
//

#include "tests.h"
#define PRINT_A 0
#define PRINT_B 1
#define PRINT_PHI 0
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

void test_fieldComputation(){
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    for(size_t i = 0; i < array_local_size.total_comp; i++)
    {
        g[i] = 0;
    }
    vortexInit(g);
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

