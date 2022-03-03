//
// Created by alcauchy on 09/02/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_UTILS_TESTS_H
#define ALLIANCE_ALPHA_1_0_UTILS_TESTS_H
#include <stdio.h>
#include <stdlib.h>
#include "utils_mpi.h"
#include "array.h"
#include "init.h"
#include "utils_hdf.h"
#include "fields.h"
#include "distrib.h"
void test_fieldComputation();
void test_freeEnergyComputation();
void test_mainFunction();
void test_fieldComparison();
void test_kSpecComputations();
#endif //ALLIANCE_ALPHA_1_0_UTILS_TESTS_H