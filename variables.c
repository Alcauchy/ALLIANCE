//
// Created by alcauchy on 08/02/2022.
//

#include "variables.h"


struct phys_params var_var;
double *var_J0;
double *var_J1;

void var_init(){
    var_varInit();
    var_getJ0();
    var_getJ1();
};

void var_getJ0(){
    var_J0 = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.ns * sizeof(*var_J0));
    for (size_t i = 0; i < array_local_size.nkx * array_local_size.nky * array_local_size.ns; i++)
    {
        var_J0[i] = exp(- 0.5 * var_var.b[i]);
    }
}

void var_getJ1(){
    var_J1 = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.ns * sizeof(*var_J1));
    for (size_t i = 0; i < array_local_size.nkx * array_local_size.nky * array_local_size.ns; i++)
    {
        if (fabs(var_var.b[i]) > 1e-16)
        {
            var_J1[i] = (1. - var_J0[i]) * 2./ var_var.b[i];
        }
        else
        {
            var_J1[i] = 1.0;
        }

    }
}

size_t var_getJIndex(size_t ikx, size_t iky, size_t is){
    return ikx * array_local_size.nky * array_local_size.ns +
           iky * array_local_size.ns +
           is;
};

void var_varInit(){
    var_var.m = malloc(parameters.ns * sizeof(*var_var.m));
    var_var.q = malloc(parameters.ns * sizeof(*var_var.q));
    var_var.n = malloc(parameters.ns * sizeof(*var_var.n));
    var_var.T = malloc(parameters.ns * sizeof(*var_var.T));
    var_var.vT = malloc(parameters.ns * sizeof(*var_var.vT));
    var_var.rho = malloc(parameters.ns * sizeof(*var_var.rho));
    var_var.b = malloc(parameters.ns * array_local_size.nky * array_local_size.nkx * sizeof(*var_var.b));

    var_var.beta = parameters.beta;
    var_var.B0 = 1.0;
    for(size_t i = 0; i < parameters.ns; i++)
    {
        var_var.m[i] = parameters.mass[i];
        var_var.q[i] = parameters.charge[i];
        var_var.T[i] = parameters.temperature[i];
        var_var.n[i] = parameters.density[i];
        var_var.vT[i] = sqrt(2.* var_var.T[i]/var_var.m[i]);
        var_var.rho[i] = sqrt(var_var.m[i]/var_var.m[0]);
    }

    for(size_t ix = 0; ix < array_local_size.nkx; ix++)
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t is = 0; is < array_local_size.ns; is++)
            {
                var_var.b[var_getJIndex(ix,iy,is)] = 0.5 * space_kPerp2[ix * array_local_size.nky + iy] *
                                                     var_var.rho[is] * var_var.rho[is];
            }
        }
    }
}