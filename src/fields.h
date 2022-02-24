//
// Created by alcauchy on 26/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_FIELDS_H
#define ALLIANCE_ALPHA_1_0_FIELDS_H
#include <stdio.h>
#include <stdlib.h>
#include "utils_mpi.h"
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
void fields_sendG(COMPLEX *g);   // broadcasts g^0_{s0},g^1_{s0},g^0_{s1} to all the processes to compute fields locally
void fields_getFieldsFromH(COMPLEX *h00, COMPLEX *h10, COMPLEX *h01);
void fields_getAFromH(const COMPLEX* h);
void fields_getBFromH(const COMPLEX *h0, const COMPLEX *h1);
void fields_getPhiFromH(const COMPLEX* h);

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
extern struct fields_chi fields_chi;
extern int *global_nm_index;
extern COMPLEX *g00;
extern COMPLEX *g10;
extern COMPLEX *g01;
extern double *A_denom;
#endif //ALLIANCE_ALPHA_1_0_FIELDS_H
