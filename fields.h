//
// Created by alcauchy on 26/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_FIELDS_H
#define ALLIANCE_ALPHA_1_0_FIELDS_H
#include <stdio.h>
#include <stdlib.h>
#include "mpi_utils.h"
#include "array.h"
void fields_init();
void fields_getA();
void fields_getB();
void fields_getPhi();
void fields_getFields();
void fields_getFieldsElectrostatic();
void fields_getFieldsElectromagnetic();
void fields_getChiElectrostatic();
void fields_getChiElectromagnetic();

extern void (*fields_getChi)(void);
#endif //ALLIANCE_ALPHA_1_0_FIELDS_H
