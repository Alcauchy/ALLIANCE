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
void fields_sendF(COMPLEX *f);   // broadcasts g^0_{s0},g^1_{s0},g^0_{s1} to all the processes to compute fields locally
void fields_getFieldsFromH(COMPLEX *h00, COMPLEX *h10, COMPLEX *h01);
void fields_getAFromH(const COMPLEX* h);
void fields_getBFromH(const COMPLEX *h0, const COMPLEX *h1);
void fields_getPhiFromH(const COMPLEX* h);
void fields_getGradY(COMPLEX *out);
void fields_getGradX(COMPLEX *out);

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
extern COMPLEX *f00;
extern COMPLEX *f10;
extern COMPLEX *f01;
extern double *A_denom;
#endif //ALLIANCE_ALPHA_1_0_FIELDS_H
