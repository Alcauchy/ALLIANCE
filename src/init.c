/**************************************
* @file init.c
*   \brief initialization module for alliance.
*
*   all the inititalization routines are here.
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 03/02/2022 created by Gene Gorbunov
//                                   INITIALIZATION
//
// init_start
// init_printParameters
// init_initFlags
// init_conditions
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "init.h"
#include "distrib.h"

/**************************************
* \def RANK_IO
*   defines rank of the processor
*   used to output information to console
***************************************/
#define RANK_IO 0

enum adiabatic kinetic;
enum electromagnetic systemType;
enum initial initialConditions;
enum spectrum spectrumType;
enum dealiasing dealiasingType;

/****************************************
* \fn void init_start(char *filename)
* \brief initialization of ALLIANCE
* \param filename: specifies parameter filename
*
* initializes all the modules required for ALLIANCE to work.
 ***************************************/
void init_start(char *filename){
    mpi_init();
    read_parameters(filename);
    init_computation();
    init_physicalSystem();
    hdf_createFiles();
    init_printParameters();
};

/****************************************
* \fn void init_physicalSystem()
* \brief initialization of physical system and parameters
 ***************************************/
void init_physicalSystem(){
    space_init();
    equation_init();
    diag_initSpec();
    var_init();
    fields_init();
    solver_init();
};

/****************************************
* \fn void init_computation()
* \brief initialize hdf, fftw and mpi
 ***************************************/
void init_computation(){
    srand(0);
    init_initFlags();
    mpi_generateTopology();
    mpi_initMExchange();
    fftw_init(mpi_kx_comm);
    hdf_init();
};

/***************************************
 * \fn void init_printParameters()
 * \brief parameter output
 *
 * prints parameters of the simulation
 ***************************************/
void init_printParameters(){
    if (mpi_my_rank == RANK_IO)
    {
        printf("===========SYSTEM PARAMETERS=========\n");
        switch(kinetic)
        {
            case NONADIABATIC:
                printf("system is fully kinetic\n");
                break;
            case ADIABATIC:
                printf("system is adiabatic\n");
                break;
        }
        switch(systemType)
        {
            case ELECTROMAGNETIC:
                printf("system is electromagnetic\n");
                break;
            case ELECTROSTATIC:
                printf("system is electrostatic\n");
                break;
        }
        switch(spectrumType)
        {
            case UNIT:
                printf("spectrum computed using unit k\n");
                break;
            case SHELL:
                printf("spectrum computed using shells\n");
                break;
        }
        switch(dealiasingType)
        {
            case ALIASED:
                printf("no dealiasing is performed\n");
                break;
            case TWOTHIRDS:
                printf("2/3 rule used for dealiasing\n");
                break;
        }
    }

};

/***************************************
 * \fn void init_initFlags()
 * \brief flag initialization
 *
 * initializes flags, which are then used
 * to define if system is adiabatic or kinetic,
 * electromagnetic or electrostatic,
 * and type of initial conditions
 ***************************************/
void init_initFlags(){
   if (parameters.postprocess){
       parameters.initial = 1;
       parameters.compute_k = 1;
       parameters.compute_m = 1;
       parameters.save_diagnostics = 1;
       parameters.save_EMfield = 1;
       parameters.compute_nonlinear = 1;
   }
   kinetic = parameters.adiabatic;
   systemType = parameters.electromagnetic;
   initialConditions = parameters.initial;
   spectrumType = parameters.spectrum;
   dealiasingType = parameters.dealiasing;
};

/***************************************
 * \fn void fill_rand(COMPLEX *data)
 * \brief fills the inital conditions randomly
 * \param data: complex 6d array to fill
 * initializes distribution with spectrum defined in
 * #init_energySpec
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_rand(COMPLEX *ar1) {
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0;
    }
    for (size_t ix = 0; ix < array_global_size.nkx; ix++){
        for (size_t iy = 0; iy < array_global_size.nky; iy++){
            for (size_t iz = 0; iz < array_global_size.nkz; iz++){
                for (size_t im = 0; im < array_global_size.nm; im++){
                    for (size_t il = 0; il < array_global_size.nl; il++){
                        for (size_t is = 0; is < array_global_size.ns; is++){
                            if (mpi_whereIsX[2*ix] == mpi_my_kx_rank && mpi_whereIsM[2 * im] == mpi_my_m_rank && ((im == 0) || (im == 1))){//mpi_whereIsX[2*ix] == mpi_my_kx_rank && mpi_whereIsM[2*im] == mpi_my_m_rank && (im == 0)
                                size_t ix_local = mpi_whereIsX[2 * ix + 1];
                                size_t im_local = mpi_whereIsM[2 * im + 1];
                                size_t ind6D = get_flat_c(is,il,im_local,ix_local,iy,iz);
                                size_t ind3D = ix_local * array_local_size.nky +
                                               iy;
                                double theta = 2. * M_PI * (double) rand() / (double) (RAND_MAX);
                                ar1[ind6D] = cexp(1.j * theta) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                                double amplitude = sqrt(
                                        init_energySpec(space_kPerp[ind3D], im_local, 1e-9, 100., space_kz[iz], 3.) / 2.0 / M_PI);
                                ar1[ind6D] *=amplitude;
                                if(global_nkx_index[ix_local] == 0 && iy == 0 && iz == 0) ar1[ind6D] = 0;
                            }
                            else{
                                rand();
                            }
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn void fill_rand(COMPLEX *data)
 * \brief fills the initial conditions randomly
 * \param data: complex 6d array to fill
 * generates random spectrum in band defined by forcing parameters.
 * Only initializes ions. Amplitude is taken to be W = eps * t_nonl
 *
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_forceDependant(COMPLEX *ar1) {
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0;
    }
    double kf = 0.5 * (parameters.forceKmin + parameters.forceKmax);
    double tau_nonl = 2.0 * M_PI * M_PI * var_var.m[1] / (kf * kf * parameters.forcePower);
    tau_nonl = pow(tau_nonl, 1.0 / 3.0);
    //printf("tau_nonl = %f\n",tau_nonl);
    double Wi = tau_nonl * parameters.forcePower;
    //printf("W_i = %f\n",Wi);
    //amount of z modes to be initialized
    int zN = 2;
    double amplitude = sqrt(2 * Wi / equation_forceNorm / zN);
    //printf("amp = %f\n",amplitude);
    double tot = 0;
    for (size_t ix = 0; ix < array_global_size.nkx; ix++){
        for (size_t iy = 0; iy < array_global_size.nky; iy++){
            for (size_t iz = 0; iz < array_global_size.nkz; iz++){
                for (size_t im = 0; im < array_global_size.nm; im++){
                    for (size_t il = 0; il < array_global_size.nl; il++){
                        for (size_t is = 0; is < array_global_size.ns; is++){
                            if (mpi_whereIsX[2 * ix] == mpi_my_kx_rank && mpi_whereIsM[2 * im] == mpi_my_m_rank && (im == (size_t) parameters.forceM) && il == 0 && iz < zN){
                                size_t ix_local = mpi_whereIsX[2 * ix + 1];
                                size_t im_local = mpi_whereIsM[2 * im + 1];
                                size_t ind6D = get_flat_c(is,il,im_local,ix_local,iy,iz);
                                size_t ind2D = ix_local * array_local_size.nky + iy;
                                double theta = 2. * M_PI * (double) rand() / (double) (RAND_MAX);
                                if (is == 1 && (space_kPerp2[ind2D] > parameters.forceKmin * parameters.forceKmin && space_kPerp2[ind2D] < parameters.forceKmax * parameters.forceKmax)){
                                    ar1[ind6D] = cexp(1.j * theta);
                                    ar1[ind6D] *= amplitude;
                                    if (iz != 0) ar1[ind6D] *= 0.5;
                                }
                                else{
                                    ar1[ind6D] = 0.j;
                                }
                                if(global_nkx_index[ix_local] == 0 && iy == 0 && iz == 0) ar1[ind6D] = 0;
                                tot+=ar1[ind6D]*conj(ar1[ind6D]) * 0.5;

                            }
                            else{
                                rand();
                            }

                        }
                    }
                }
            }
        }
    }
    //printf("tot = %f\n",tot);
}


/***************************************
 * \fn void fill_randM0(COMPLEX *data)
 * \brief fill zeroth Hermite moment with random values
 * \param data: complex 6D array to fill
 *
 * fills 0-th Hermite moment of a distribution function ar1 with random values
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_randM0(COMPLEX *ar1) {
    size_t ind6D;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0.;
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t il = 0; il < array_local_size.nl; il++) {
                    for (size_t im = 0; im < array_local_size.nm; im++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            if (space_globalMIndex[im] == 0 || space_globalMIndex[im] == 1) {
                                //printf("%d\n",im);
                                ind6D = get_flat_c(is, il, im, ix, iy, iz);
                                ar1[ind6D] = ((0.5 - (double) rand() / (double) (RAND_MAX)) +
                                             (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                            }

                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn void fill_randSingleKM(COMPLEX *ar1)
 * \brief fill single chosen wavevector and Hermite moment
 * \param data: complex 6D array
 *
 * initializes single wavevector and Hermite moment
 * of a distribution function with random variable.
 * This function is only for in-module use
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_randSingleKM(COMPLEX *ar1) {
    size_t ind6D;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0.;
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
       // printf("[MPI process %d] ix = %zu, %d\n",mpi_my_rank, ix, global_nkx_index[ix]);
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    //printf("im = %d\n", space_globalMIndex[im]);
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            if (space_globalMIndex[im] == 1 || space_globalMIndex[im] == 0 && global_nkx_index[ix] == 0) {
                                ind6D = get_flat_c(is, il, im, ix, 0, iz);
                                ar1[ind6D] = ((0.5 - (double) rand() / (double) (RAND_MAX)) +
                                                   (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j)
                                                           * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                            };
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn void init_conditions(COMPLEX *data)
 * \brief distribution function initialization
 * \param data: complex 6D array
 *
 * initializes distribution function with chosen method
 * (see #fill_rand,
 * #fill_randM0,
 * #fill_randSingleKM)
 ***************************************/
void init_conditions(COMPLEX *data){
    switch(initialConditions)
    {
        case RANDOM:
            //fill_rand(data);
            //init_localPertirbation(data);
            fill_forceDependant(data);
            //init_reconnection(data);
            //init_fillSinc(data);
            //fill_randM0(data);
            //init_fillCosZ(data);
            //fill_randSingleKM(data);
            break;

        case FROMFILE:
            hdf_readData(parameters.from_simulationName,data);
            break;

        default:
            printf("[MPI process %d] error with initial conditions! Aborting...",mpi_my_rank);
            exit(1);
    }
    distrib_dealiasing(data);
    distrib_enforceReality(data);
    distrib_setZeroNHalf(data);
};

/***************************************
 * \fn double init_energySpec(double k, double m, double amp, double disp)
 * \brief returns energy spectrum
 * \param k: a wavenumber at which spectrum is computed
 * \param m: Hermite moment at which amplitude is computed
 * \param amp: amplitude of the spectrum
 * \param disp: dispersion of the spectrum
 *
 * computes spectrum of form
 * \f$ A \cdot k^2 exp(-2 k^2/\sigma^2) \f$, where
 * \f$\sigma = disp\f$, and
 * \f$A = amp \f$
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
double init_energySpec(double k, double m, double amp, double disp, double k_z, double disp_z) {
    return amp * k * k * k * k * exp( - 2.0 * (k * k / disp / disp )) * k_z * k_z * k_z * k_z * exp( - 2.0 * (k_z * k_z / disp_z / disp_z ));
};

/***************************************
 * \fn double init_sinc(double k, double m, double amp, double disp)
 * \brief returns energy spectrum
 * \param k: a wavenumber at which spectrum is computed
 * \param m: Hermite moment at which amplitude is computed
 * \param amp: amplitude of the spectrum
 * \param disp: dispersion of the spectrum
 *
 * computes sinc function
 * \f$ A \cdot k^2 exp(-2 k^2/\sigma^2) \f$, where
 * \f$\sigma = disp\f$, and
 * \f$A = amp \f$
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
double init_sinc(double amp, double f, double x, double y, double z, double x0, double y0, double z0) {
    double valueX, valueY, valueZ;
    valueX = (fabs((x - x0)) > 1e-13) ? sin(f * (x - x0)) / f/(x - x0) : 1.;
    valueY = (fabs((y - y0)) > 1e-13) ? sin(f * (y - y0)) / f /(y - y0) : 1.;
    valueZ = (fabs((z - z0)) > 1e-13) ? sin(f * (z - z0)) / f /(z - z0) : 1.;

    return amp * valueX * valueY * valueZ;
}

/***************************************
 * \fn void init_fillSinc(COMPLEX *out)
 * \brief returns energy spectrum
 * \param k: a wavenumber at which spectrum is computed
 * \param m: Hermite moment at which amplitude is computed
 * \param amp: amplitude of the spectrum
 * \param disp: dispersion of the spectrum
 *
 * computes sinc function
 * \f$ A \cdot k^2 exp(-2 k^2/\sigma^2) \f$, where
 * \f$\sigma = disp\f$, and
 * \f$A = amp \f$
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void init_fillSinc(COMPLEX *out) {
    size_t proc_id = mpi_whereIsM[0];
    size_t local_m = mpi_whereIsM[1];
    size_t ind6D;
    double *real = fftw_hBuf;
    double *real_check = calloc(array_local_size.total_real, sizeof(*real_check));
    double *real_check1 = calloc(array_local_size.total_real, sizeof(*real_check1));
    double x0,y0,z0;
    double x,y,z;
    double disp[2] = {0.05,0.08};
    x0 = space_dx * array_global_size.nx/2;
    y0 = space_dy * array_global_size.ny/2;
    z0 = space_dz * array_global_size.nz/2;
    if(mpi_my_m_rank == proc_id){
        for(size_t iy = 0; iy < array_local_size.ny; iy++){
            for(size_t ix = 0; ix< array_local_size.nx; ix++){
                for(size_t iz = 0; iz <  array_local_size.nz + 2; iz++){
                    for(size_t is = 0; is< array_local_size.ns; is++){
                        ind6D = get_flat_r(is,0,local_m,ix,iy,iz);
                        x = space_dx * (ix);
                        y = space_dy * (iy + array_global_size.ny * mpi_my_kx_rank / mpi_dims[1]);
                        z = space_dz * (iz);
                        //real[ind6D] = init_sinc(0.01,4,x,y,z,x0,y0,z0);
                        real[ind6D] = init_exp2(1.,disp[is],x,y,z,x0,y0,z0);
                       //real[ind6D] = init_sinX(10.0, 3.0 / space_Lx,x);//*init_sinX(10.0, 3.0 / space_Ly,y);;
                       //real[ind6D] = init_sinX(10.0, 4.0 *disp[is]/ space_Lx,x) * init_cosX(1.0, 6.0/ space_Ly,y);;
                    }
                }
            }
        }
    }
    hdf_create_file_r("init.h5",real);
    fftw_r2c();
    fftw_copy_buffer_c(out,fftw_hBuf);
    //for(size_t ii = 0; ii < array_local_size.total_comp; ii++) out[ii] = 1.;
}

/***************************************
 * \fn double init_exp2(double k, double m, double amp, double disp)
 * \brief returns energy spectrum
 * \param k: a wavenumber at which spectrum is computed
 * \param m: Hermite moment at which amplitude is computed
 * \param amp: amplitude of the spectrum
 * \param disp: dispersion of the spectrum
 *
 * computes exp2 function
 * \f$ A \cdot k^2 exp(-2 k^2/\sigma^2) \f$, where
 * \f$\sigma = disp\f$, and
 * \f$A = amp \f$
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
double init_exp2(double amp, double f, double x, double y, double z, double x0, double y0, double z0) {
    double valueX, valueY, valueZ;
    valueX = exp(-f * ((x - x0) * (x - x0)));
    valueY = exp(-f * ((y - y0) * (y - y0)));
    valueZ = exp(-f * ((z - z0) * (z - z0)));

    return amp * valueX * valueY * valueZ;
}

void init_fillCosZ(COMPLEX *out){
    int mProc;
    int mLocal;
    mProc = mpi_whereIsM[0];
    mLocal = mpi_whereIsM[1];
    double *cosZ = calloc(array_local_size.total_real, sizeof(*cosZ));
    if (mpi_my_m_rank == mProc){
        for (size_t iy = 0; iy < array_local_size.ny; iy ++){
            for (size_t ix = 0; ix < array_local_size.nx; ix ++){
                for (size_t iz = 0; iz < array_local_size.nz; iz ++){
                    for (size_t il = 0; il < array_local_size.nl; il ++){
                        for (size_t is = 0; is < array_local_size.ns; is ++){
                            size_t ind6D = get_flat_r(is,il,mLocal,ix,iy,iz);
                            size_t iyInd = global_ny_index[iy];
                            cosZ[ind6D] = init_cosZ(var_var.q[is] * 1.0,1.0,iz * space_dz)
                                    * init_sinY(1.0,1.0,iyInd * space_dy)
                                    * init_cosX(1.0,1.0,ix * space_dx);

                        }
                    }
                }
            }
        }
    }
    fftw_copy_buffer_r(fftw_hBuf,cosZ);
    fftw_r2c();
    fftw_copy_buffer_c(out,fftw_hBuf);
    free(cosZ);

}

void init_localPertirbation(COMPLEX *out){
    int mProc;
    int mLocal;
    mProc = mpi_whereIsM[0];
    mLocal = mpi_whereIsM[1];
    double *realInit = calloc(array_local_size.total_real, sizeof(*realInit));

    if(mpi_my_m_rank == mProc){
        for (size_t iy = 0; iy < array_local_size.ny; iy ++){
            for (size_t ix = 0; ix < array_local_size.nx; ix ++){
                for (size_t iz = 0; iz < array_local_size.nz; iz ++){
                    for (size_t il = 0; il < array_local_size.nl; il ++){
                        for (size_t is = 0; is < array_local_size.ns; is ++){
                            size_t ind6D = get_flat_r(is,il,mLocal,ix,iy,iz);
                            size_t iyInd = global_ny_index[iy];
                            realInit[ind6D] = init_cosZ(var_var.q[is] * 1.0, 1.0, iz * space_dz)
                                    * init_exp2(60.0, 2.0, ix * space_dx, iyInd * space_dy, 0.0, M_PI, M_PI, 0);
                        }
                    }
                }
            }
        }
    }
    fftw_copy_buffer_r(fftw_hBuf, realInit);
    fftw_r2c();
    fftw_copy_buffer_c(out,fftw_hBuf);
    //hdf_create_file_r("init.h5",realInit);
    free(realInit);
}

/***************************************
 * \fn double init_exp2(double k, double m, double amp, double disp)
 * \brief returns energy spectrum
 * \param x: a wavenumber at which spectrum is computed
 * \param amp: amplitude of the spectrum
 *
 * computes sinus function along chosen axis
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
double init_sinY(double amp, double f, double y){
    return amp * cos( 2 * M_PI * f / space_Ly * y);
}

double init_cosX(double amp, double f, double x){
    return amp * cos( 2 * M_PI * f / space_Lx * x);
}

double init_cosZ(double amp, double f, double x){
    return amp * cos( 2 * M_PI * f / space_Lz * x);
}

void init_reconnection(COMPLEX *out){
    size_t ind6D;
    double *Apar0 = calloc(array_local_size.total_real, sizeof(*Apar0));
    double A0 = 1.299;
    double Aper = 0e-8;
    int local_m = mpi_whereIsM[3];
    int proc_m = mpi_whereIsM[2];
    COMPLEX *Apar0Comp = Apar0;
    // obtain A_par^0 in real space
    if (mpi_my_m_rank == proc_m){
        for (size_t iy = 0; iy < array_local_size.ny; iy ++){
            for (size_t ix = 0; ix < array_local_size.nx; ix ++){
                for (size_t iz = 0; iz < array_local_size.nz; iz ++){
                    // only l=0 m=1 moments are initialized with non-zeros to initialize only A_parallel field. ions are considered to be 0
                    ind6D = get_flat_r(0, 0, local_m, ix, iy, iz);
                    double y = space_dy * (iy + array_global_size.ny * mpi_my_kx_rank / mpi_dims[1]);
                    double cosh = init_cosh(ix, 200.0);
                    double shape = init_shapeFunction(ix);
                    Apar0[ind6D] = 0*(A0 / cosh / cosh * shape + Aper * init_exp2(1.0,60.,ix * space_dx,0,0,M_PI,0,0) *
                                                                      init_sinY(1.0, 1.0, y));
                }
            }
        }
    }
    double *fftw_d_field = fftw_field;
    for (size_t iy = 0; iy < array_local_size.ny; iy ++){
        for (size_t ix = 0; ix < array_local_size.nx; ix ++){
            for (size_t iz = 0; iz < array_local_size.nz; iz ++){
                // only l=0 m=1 moments are initialized with non-zeros to initialize only A_parallel field. ions are considered to be 0
                ind6D = iy * array_local_size.nx * (array_local_size.nz+2) + ix * (array_local_size.nz+2) + iz;
                ind6D = get_flat_r(0, 0, local_m, ix, iy, iz);
                double y = space_dy * (iy + array_global_size.ny * mpi_my_kx_rank / mpi_dims[1]);
                double cosh = init_cosh(ix, 1.0);
                double shape = init_shapeFunction(ix);
                fftw_d_field[ind6D] = A0 / cosh / cosh * shape + Aper * init_exp2(1.0,1.,ix * space_dx,0,0,M_PI,0,0) *
                                                          init_sinY(1.0, 1.0, y);
            }
        }
    }
    //hdf_create_file_r("init.h5",Apar0);

    //FFT and filter at kx[30]:
    //fftw_copy_buffer_r(fftw_hBuf,Apar0);
    fftw_r2c_field();
    fftw_copyFieldBuf_c(fields_fields.A,fftw_field);
    //hdf_saveFieldA("initF.h5");
    //double cutoff = space_dKx * 5.;
    //printf("cutoff = %f\n",cutoff);
    //diag_filterK(fftw_hBuf, Apar0Comp, cutoff);
    //fftw_copy_buffer_c(fftw_hBuf,Apar0Comp);
    //fftw_c2r();
    //hdf_create_file_r("init.h5",fftw_hBuf);
    free(Apar0);
    //now initialize h
    if (mpi_my_m_rank == proc_m){
        for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
            for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                    ind6D = get_flat_c(0, 0, local_m, ix, iy, iz);
                    size_t ind2D = ix * array_local_size.nky + iy;
                    size_t ind3D = ix * array_local_size.nky * array_local_size.ns +
                                   iy * array_local_size.ns;
                    size_t ind3Df = get_flatIndexComplex3D(ix,iy,iz);
                    double beta_coef = 2. * space_kPerp2[ind2D] / var_var.beta;
                    double var_coef = sqrt(2.0) / (var_var.q[0] * var_var.n[0] * var_var.vT[0] * var_J0[ind3D]);//sqrt(2.0) / (var_var.q[0] * var_var.n[0] * var_var.vT[0] * var_J0[ind3D]);

                    fftw_hBuf[ind6D] = (beta_coef * var_coef) *fields_fields.A[ind3Df];//(beta_coef * var_coef) * fields_fields.A[ind3Df];
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    fftw_copy_buffer_c(out, fftw_hBuf);
    fftw_c2r(fftw_hBuf);
    hdf_create_file_r("init.h5",fftw_hBuf);
}

double init_cosh(double x, double a) {
    double xN = x * space_dx - space_Lx/2.;
    return 0.5 * (exp(xN/a) + exp(-xN/a));
}

double init_tanh(double x){
    return (exp(x) - exp(-x)) / (exp(-x) + exp(x));
}

double init_shapeFunction(double x){
    double xN = x * space_dx;
    double Sh;
    Sh = init_tanh(xN * 2. * M_PI / space_Lx) * init_tanh(xN * 2. * M_PI / space_Lx);
    Sh += init_tanh(xN * 2. * M_PI / space_Lx - 2 * M_PI) * init_tanh(xN * 2. * M_PI / space_Lx - 2 * M_PI);
    Sh -= init_tanh(2*M_PI) * init_tanh(2*M_PI);
    Sh /= (2 * init_tanh(M_PI) * init_tanh(M_PI) - init_tanh(2 * M_PI) * init_tanh(2 * M_PI));
    return Sh;
}