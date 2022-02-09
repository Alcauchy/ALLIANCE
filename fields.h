//
// Created by alcauchy on 26/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_FIELDS_H
#define ALLIANCE_ALPHA_1_0_FIELDS_H
#include <stdio.h>
#include <stdlib.h>
#include "mpi_utils.h"
#include "array.h"
#include "init.h"
void fields_init();
void fields_getA(const COMPLEX* g);
void fields_getB(const COMPLEX* g0, const COMPLEX* g1);
void fields_getPhi(const COMPLEX* g0, const COMPLEX* g1);
void fields_getFields(COMPLEX *g00, COMPLEX *g10, COMPLEX *g01);
void fields_getChi();
void fields_getChiPhi();
void fields_getChiB();
void fields_getChiA();


struct fields_fields{
    COMPLEX *phi;
    COMPLEX *A;
    COMPLEX *B;
};

struct fields_chi{
    COMPLEX *phi;
    COMPLEX *A;
    COMPLEX *B;
};
extern struct fields_fields fields_fields;

#endif //ALLIANCE_ALPHA_1_0_FIELDS_H
