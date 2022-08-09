/**************************************
* @file equation.c
*   \brief equation module
*
*   module required to compute RHS of the equation.
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 26/01/2022 created by Gene Gorbunov
//                                   EQUATION
//
// equation_getLinearTerm
// equation_getNonlinearElectromagnetic
// equation_getNonlinearElectrostatic
// equation_getNonlinearTerm
// equation_getRHS
// equation_init
// equation_getForcing
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "equation.h"
#include "array.h"
#include "fields.h"
#include "distrib.h"

#define CHI_EM 3
#define CHI_EL 1
#define CHI_PHI 0
#define CHI_A 1
#define CHI_B 2

size_t *equation_forceKxInd;
size_t *equation_forceKyInd;
size_t *equation_forceKzInd;
int equation_forceKn;
int equation_forceNorm;
double equation_forcingCoef;

/***************************************
 * \fn void equation_getLinearTerm(const COMPLEX *in, const COMPLEX *plus_boundary, const COMPLEX *minus_boundary, COMPLEX *out)
 * \brief computes linear term
 * \param in: complex array
 * \param out: complex array
 * \param plus_boundary: complex array
 * \param minus_boundary: complex array
 *
 * computes linear term <tt>out</tt> from distribution function <tt> in </tt>.
 ***************************************/
void equation_getLinearTerm(const COMPLEX *in, const COMPLEX *plus_boundary, const COMPLEX *minus_boundary, COMPLEX *out) {
    size_t ind6D;
    size_t ind6DPlus;
    size_t ind6DMinus;
    size_t indBoundary;
    //case if there is only one Hermite moment per process
    if (array_local_size.nm == 1){
        for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
            for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            indBoundary = ix * array_local_size.nky * array_local_size.nkz * array_local_size.nl *
                                          array_local_size.ns
                                          + iy * array_local_size.nkz * array_local_size.nl * array_local_size.ns
                                          + iz * array_local_size.nl * array_local_size.ns
                                          + il * array_local_size.ns
                                          + is;
                            ind6D = get_flat_c(is, il, 0, ix, iy, iz);
                            out[ind6D] += space_iKz[iz] * var_var.vT[is] * (space_sqrtM[0] * minus_boundary[indBoundary] +
                                                                            space_sqrtM[array_local_size.nm] * plus_boundary[indBoundary]);
                        }
                    }
                }
            }
        }
    }
    //case if there is only more than one Hermite moment per process
    else{
        for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
            for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                    for (size_t im = 1; im < array_local_size.nm - 1; im++) {
                        for (size_t il = 0; il < array_local_size.nl; il++) {
                            for (size_t is = 0; is < array_local_size.ns; is++) {
                                ind6D = get_flat_c(is, il, im, ix, iy, iz);
                                ind6DPlus = get_flat_c(is, il, im + 1, ix, iy, iz);
                                ind6DMinus = get_flat_c(is, il, im - 1, ix, iy, iz);
                                out[ind6D] += space_iKz[iz] * var_var.vT[is] *
                                              (space_sqrtM[im + 1] * in[ind6DPlus] + space_sqrtM[im] * in[ind6DMinus]);

                            }
                        }
                    }
                }
            }
        }
        for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
            for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            /*treating left boundary*/
                            ind6D = get_flat_c(is, il, 0, ix, iy, iz);
                            ind6DPlus = get_flat_c(is, il, 1, ix, iy, iz);
                            indBoundary = ix * array_local_size.nky * array_local_size.nkz * array_local_size.nl *
                                          array_local_size.ns
                                          + iy * array_local_size.nkz * array_local_size.nl * array_local_size.ns
                                          + iz * array_local_size.nl * array_local_size.ns
                                          + il * array_local_size.ns
                                          + is;
                            out[ind6D] += space_iKz[iz] * var_var.vT[is] *
                                          (space_sqrtM[1] * in[ind6DPlus] + space_sqrtM[0] * minus_boundary[indBoundary]);
                            /*treating right boundary*/
                            ind6D = get_flat_c(is, il, array_local_size.nm - 1, ix, iy, iz);
                            ind6DMinus = get_flat_c(is, il, array_local_size.nm - 2, ix, iy, iz);
                            out[ind6D] += space_iKz[iz] * var_var.vT[is] * (space_sqrtM[array_local_size.nm - 1] * in[ind6DMinus] +
                                                                            space_sqrtM[array_local_size.nm] * plus_boundary[indBoundary]);
                        }
                    }
                }
            }
        }
    }

};

/***************************************
 * \fn void equation_getNonlinearElectromagnetic(double *in, double *chiAr, double *out, double sign)
 * \brief returns nonlinear electromagnetic term
 * \param in: input double array
 * \param chiAr: input double array
 * \param out: output double array
 * \param sign: should be 1 or -1
 *
 * performs multiplication between input 6D complex array <tt>in</tt> and
 * gyrokinetic potential array <tt>chiAr</tt>, in such way that the structure of the product
 * is the same as nonlinear term of drift-kinetic equations.
 * Used by #equation_getNonlinearProduct. <tt>sign</tt> is used to determine the sign of the resulting product.
 * See #equation_getNonlinearTerm for explanation.
 ***************************************/
void equation_getNonlinearElectromagnetic(double *in, double *chiAr, double *out, double sign) {
    size_t indH;
    size_t indHPlus;
    size_t indHMinus;
    size_t indHBound;
    size_t indH_l1;
    size_t indH_l0;
    size_t indChi_phi;
    size_t indChi_A;
    size_t indChi_B;
    double phibh = 0;
    double ah = 0;
    double bh = 0;
    double *plus_boundary = calloc(array_local_size.nx *
                                    array_local_size.ny *
                                    (array_local_size.nz + 2)*
                                    array_local_size.nl *
                                    array_local_size.ns,
                                    sizeof(*plus_boundary));
    double *minus_boundary = calloc(array_local_size.nx *
                                    array_local_size.ny *
                                    (array_local_size.nz + 2)*
                                    array_local_size.nl *
                                    array_local_size.ns,
                                    sizeof(*plus_boundary));
    mpi_exchangeMBoundaries_r(in, plus_boundary, minus_boundary);
    for(size_t iy = 0; iy < array_local_size.ny; iy++){
        for(size_t ix = 0; ix < array_local_size.nx; ix++){
            for(size_t iz = 0; iz < array_local_size.nz + 2; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            //getting flat indices
                            indH = get_flat_r(is,il,im,ix,iy,iz);
                            indH_l1 = get_flat_r(is,1,im,ix,iy,iz);
                            indH_l0 = get_flat_r(is,0,im,ix,iy,iz);
                            indChi_phi = getIndChiBufEM_r(ix,iy,iz,is,CHI_PHI);
                            indChi_A = getIndChiBufEM_r(ix,iy,iz,is,CHI_A);
                            indChi_B = getIndChiBufEM_r(ix,iy,iz,is,CHI_B);
                            // computing (chi_phi + chi_b * h)
                            phibh = (chiAr[indChi_phi] + chiAr[indChi_B]) * in[indH];
                            // computing chi_b * h_l=1
                            bh = (il==0)?  chiAr[indChi_B] * in[indH_l1] : chiAr[indChi_B] * (in[indH_l0] + 2. * in[indH_l1]);
                            if(array_local_size.nm == 1){
                                indHBound = iy * array_local_size.nx * (array_local_size.nz + 2) * array_local_size.nl * array_local_size.ns
                                            + ix * (array_local_size.nz + 2) * array_local_size.nl * array_local_size.ns
                                            + iz * array_local_size.nl * array_local_size.ns
                                            + il * array_local_size.ns
                                            + is;
                                ah = chiAr[indChi_A] * (space_sqrtM[1] * plus_boundary[indHBound] + space_sqrtM[0] * minus_boundary[indHBound]);
                            }
                            else{
                                if(im!=0 && im!=array_local_size.nm - 1){
                                    indHPlus = get_flat_r(is,il,im + 1,ix,iy,iz);
                                    indHMinus = get_flat_r(is,il,im - 1,ix,iy,iz);
                                    ah = chiAr[indChi_A] * (space_sqrtM[im + 1] * in[indHPlus] + space_sqrtM[im] * in[indHMinus]);
                                }
                                    //left boundary
                                else if(im == 0){
                                    indHPlus = get_flat_r(is,il,im + 1,ix,iy,iz);
                                    indHBound = iy * array_local_size.nx * (array_local_size.nz + 2) * array_local_size.nl * array_local_size.ns
                                                + ix * (array_local_size.nz + 2) * array_local_size.nl * array_local_size.ns
                                                + iz * array_local_size.nl * array_local_size.ns
                                                + il * array_local_size.ns
                                                + is;
                                    ah = chiAr[indChi_A] * (space_sqrtM[0] * minus_boundary[indHBound] + space_sqrtM[1] * in[indHPlus]);
                                }
                                    //right boundary
                                else{
                                    indHMinus = get_flat_r(is,il,im - 1,ix,iy,iz);
                                    indHBound = iy * array_local_size.nx * (array_local_size.nz + 2) * array_local_size.nl * array_local_size.ns
                                                + ix * (array_local_size.nz + 2) * array_local_size.nl * array_local_size.ns
                                                + iz * array_local_size.nl * array_local_size.ns
                                                + il * array_local_size.ns
                                                + is;
                                    ah = chiAr[indChi_A] * (space_sqrtM[array_local_size.nm] * plus_boundary[indHBound] + space_sqrtM[array_local_size.nm - 1] * in[indHMinus]);
                                }
                            }
                            out[indH] += sign * (phibh + bh + ah);
                            out[indH] /= var_var.B0;
                        }
                    }
                }
            }
        }
    }
    free(minus_boundary);
    free(plus_boundary);
};

/***************************************
 * \fn void equation_getNonlinearElectrostatic(double *in, double *chiAr, double *out, double sign)
 * \brief returns nonlinear electrostatic term
 * \param in: input double array
 * \param chiAr: input double array
 * \param out: output double array
 * \param sign: should be 1 or -1
 *
 * see #equation_getNonlinearElectromagnetic for explanation
 ***************************************/
void equation_getNonlinearElectrostatic(double *in, double *chiAr, double *out, double sign) {
    size_t indH;
    size_t indChi;
    for(size_t iy = 0; iy < array_local_size.ny; iy++){
        for(size_t ix = 0; ix < array_local_size.nx; ix++){
            for(size_t iz = 0; iz < array_local_size.nz + 2; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            indH = get_flat_r(is,il,im,ix,iy,iz);
                            indChi = getIndChiBufEL_r(ix,iy,iz,is);
                            out[indH] += sign * in[indH] * chiAr[indChi];
                            out[indH] /= var_var.B0;
                        }
                    }
                }
            }
        }
    }
};

/***************************************
 * \fn equation_getNonlinearProduct(double *in, double *chiAr, double *out, double sign)
 * \brief chooses between computing electrostatic or electromagnetic term
 * \param in: input double array
 * \param chiAr: input double array
 * \param out: output double array
 * \param sign: should be 1 or -1
 *
 * depending on flag <tt>systemType</tt> provided by user in parameter file, chooses between
 * #equation_getNonlinearElectrostatic and #equation_getNonlinearElectromagnetic
 ***************************************/
void equation_getNonlinearProduct(double *in, double *chiAr, double *out, double sign){
    switch(systemType)
    {
        case ELECTROSTATIC:
            equation_getNonlinearElectrostatic(in, chiAr, out, sign);
            break;
        case ELECTROMAGNETIC:
            equation_getNonlinearElectromagnetic(in, chiAr, out, sign);
            break;
        default:
            printf("[MPI process %d] ERROR COMPUTING NONLINEAR TERM; ABORTING",mpi_my_rank);
            exit(1);
    }
}

/***************************************
 * \fn void equation_getNonlinearTerm(const COMPLEX *h, COMPLEX *out)
 * \brief computes nonlinear term
 * \param h: input complex array
 * \param out: output complex array
 *
 * function returns nonlinear term.  First it takes  y gradient of distribution function
 * <tt>h</tt>, and x gradient of gyrokinetic potentials chi, and transforms them to real space:\n
 * <tt>
 * distrib_getYGrad(h, fftw_hBuf); \n
 * fields_getGradX(fftw_chiBuf); \n
 * fftw_c2r(); \n
 * fftw_c2r_chi(); \n
 * </tt>
 * after that, it computes
 * \f$\frac{\partial h}{\partial y} \frac{\partial \chi}{\partial x}\f$
 * part of the poisson brackets:
 * <tt> \n
 * equation_getNonlinearProduct((double *)fftw_hBuf, (double *)fftw_chiBuf, buffer, 1.);
 * </tt> \n
 * with the result stored in <tt>buffer</tt>
 * after that, it computes x gradient of <tt>h</tt> and y gradient of gyrokinetic potential chi,
 * and transforms results to real space:\n
 * <tt>
 * distrib_getXGrad(h, fftw_hBuf);\n
 * fields_getGradY(fftw_chiBuf);\n
 * fftw_c2r();\n
 * fftw_c2r_chi();\n
 * </tt>
 * and computes second part of the poisson brackets
 * \f$-\frac{\partial h}{\partial x} \frac{\partial \chi}{\partial y}\f$
 * and adds the result to <tt>buffer</tt>.
 * <tt>buffer</tt> is then transformed back to Fourier space, and dealiasing is performed.
 ***************************************/
void equation_getNonlinearTerm(const COMPLEX *h, COMPLEX *out) {
    size_t ind6D;
    size_t ind4D;
    size_t indPhi;
    size_t indB;
    size_t indA;
    size_t chi_buf_size;
    size_t chi_buf_size_r;
    double *h_r = calloc(array_local_size.total_real, sizeof(*h_r));
    //
    // CHECKS DELETE LATER WHEN NOT NEEDED
    //
    fftw_copy_buffer_c(fftw_hBuf,h);
    fftw_c2r();
    fftw_copy_buffer_r(h_r,fftw_hBuf);
    double *buffer = calloc(array_local_size.total_real, sizeof(*buffer));
    //
    // take gradient: dh/dx and dchi/dy
    //
    distrib_getYGrad(h, fftw_hBuf);
    fields_getGradX(fftw_chiBuf);
    //
    // make fftw of everything
    //
    fftw_c2r();
    fftw_c2r_chi();
    //printf("[MPI process %d] fftw of dh/dx and dchi/dy is done\n",mpi_my_rank);
    //
    //compute product, first part
    //
    equation_getNonlinearProduct((double *)fftw_hBuf, (double *)fftw_chiBuf, buffer, 1.);
    double tot = 0;
    double totb = 0;
    for (size_t ii = 0; ii < array_local_size.total_real; ii++)
    {
        tot += h_r[ii] * buffer[ii];
        totb += buffer[ii];
    }
    //printf("[MPI process %d] total before minus = %.10e\n", mpi_my_rank, totb);
    //printf("[MPI process %d] total flux before minus= %.10e\n", mpi_my_rank, tot);
    //
    // take gradient: dh/dy and dchi/dx
    //
    distrib_getXGrad(h, fftw_hBuf);
    //for (size_t ii = 0; ii < array_local_size.total_comp; ii++) printf("[MPI process %d] buf = %f\n",mpi_my_rank, fftw_hBuf[ii]);
    fields_getGradY(fftw_chiBuf);
    //
    // transform those gradients
    //
    fftw_c2r();
    fftw_c2r_chi();
    //printf("[MPI process %d] fftw of dh/dy and dchi/dx is done\n",mpi_my_rank);
    //
    //compute product, second part
    //
    equation_getNonlinearProduct((double *)fftw_hBuf, (double *)fftw_chiBuf, buffer, -1.);
    //
    // make ifftw of the computed product
    //
    fftw_copy_buffer_r((double *)fftw_hBuf, buffer);
    //tot = 0;
    totb = 0;
    for (size_t ii = 0; ii < array_local_size.total_real; ii++)
    {
        tot += h_r[ii] * buffer[ii];
        totb += buffer[ii];
    }
    //printf("[MPI process %d] total after minus = %.10e\n", mpi_my_rank,totb);
    double total_reduced = 0;
    MPI_Allreduce(&tot,
                  &total_reduced,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);
    //if (mpi_my_rank == 0 ) printf("[MPI process %d] total flux after minus = %.10e\n", mpi_my_rank, total_reduced);
    fftw_r2c();

    //printf("[MPI process %d] inverse fftw of the nonlinear product is done\n",mpi_my_rank);
    //
    // add it to RHS
    //
    for (size_t ii = 0; ii < array_local_size.total_comp; ii++) out[ii] -= fftw_hBuf[ii];
    //
    // perform dealiasing
    //
    dealiasing23(out);
    free(buffer);
    free(h_r);
};


/***************************************
 * \fn void equation_getRHS()
 * \brief computes right hand side of the equation.
 ***************************************/
void equation_getRHS(const COMPLEX *in_g, COMPLEX *in_h, COMPLEX *out) {
    /* allocating memory for buffers */
    size_t boundary_size = array_local_size.nkx *
            array_local_size.nky *
            array_local_size.nkz *
            array_local_size.nl *
            array_local_size.ns;
    COMPLEX *minus_boundary = calloc(boundary_size, sizeof(*minus_boundary));
    COMPLEX *plus_boundary = calloc(boundary_size, sizeof(*plus_boundary));
    /* exchanging boundaries to compute linear term and fields */
    /* computing h */
    fields_sendF(in_g);
    fields_getFields(f00, f10, f01);
    fields_getChi();
    distrib_getH(in_h, in_g);

    /* boundary exchange */
    mpi_exchangeMBoundaries(in_h, plus_boundary, minus_boundary);
    /* computing linear term */
    equation_getLinearTerm(in_h, plus_boundary, minus_boundary, out);
    /* computing dissipation and forcing*/
    equation_getForcing(in_h, out);
    equation_getDissipation(in_h, out);
    /* computing nonlinear term */
    equation_getNonlinearTerm(in_h, out);
    free(minus_boundary);
    free(plus_boundary);
};

/***************************************
 * \fn void equation_getDissipation()
 * \brief computes dissipation.
 ***************************************/
void equation_getDissipation(const COMPLEX *h, COMPLEX *rhs) {
    size_t ind6D;
    size_t ind2D;
    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for (size_t iy = 0; iy < array_local_size.nky; iy++){
            for (size_t iz = 0; iz < array_local_size.nkz; iz++){
                for (size_t im = 0; im < array_local_size.nm; im++){
                    for (size_t il = 0; il < array_local_size.nl; il++){
                        for (size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            ind2D = ix * array_local_size.nky + iy;
                            rhs[ind6D] -= (var_var.mu_k * space_kPerp2[ind2D] + var_var.mu_m * global_nm_index[im])  * h[ind6D] ;
                        }
                    }
                }
            }
        }
    }
};

/***************************************
 * \fn void equation_init()
 * \brief initializes forcing.
 ***************************************/
void equation_init() {
    double k_min = parameters.forceKmin;
    double k_max = parameters.forceKmax;
    size_t ind3D;
    equation_forceKn = 0;
    equation_forceNorm = 0;

    // find out how many wave numbers will be forced
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                ind3D = ix * array_local_size.nky * array_local_size.nkz +
                        iy * array_local_size.nkz +
                        iz;
                if(space_kSq[ind3D] > k_min * k_min && space_kSq[ind3D] < k_max * k_max){
                    equation_forceKn += 1;
                }
            }
        }
    }
    printf("[MPI process %d] forcing %d modes\n",mpi_my_rank,equation_forceKn);
    MPI_Allreduce(&equation_forceKn,&equation_forceNorm,1,MPI_INT, MPI_SUM, mpi_row_comm);
    if (mpi_my_rank == 0) printf("%d\n",equation_forceNorm);

    // putting list of indices of wave numbers which are excited
    equation_forceKxInd = malloc(equation_forceKn * sizeof(*equation_forceKxInd));
    equation_forceKyInd = malloc(equation_forceKn * sizeof(*equation_forceKyInd));
    equation_forceKzInd = malloc(equation_forceKn * sizeof(*equation_forceKzInd));
    size_t count = 0;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                ind3D = ix * array_local_size.nky * array_local_size.nkz +
                        iy * array_local_size.nkz +
                        iz;
                if(space_kSq[ind3D] > k_min * k_min && space_kSq[ind3D] < k_max * k_max){
                    equation_forceKxInd[count] = ix;
                    equation_forceKyInd[count] = iy;
                    equation_forceKzInd[count] = iz;
                    count++;
                }
            }
        }
    }
    //forcing coefficient
    equation_forcingCoef = parameters.forcePower /
                           equation_forceNorm *
                           (array_global_size.nx * array_global_size.ny * array_global_size.nz);
};

/***************************************
 * \fn void equation_getForcing()
 * \brief adds forcing to the RHS.
 ***************************************/
void equation_getForcing(const COMPLEX *h, COMPLEX *rhs) {
    size_t local_m;
    size_t forced_proc;
    size_t indKx;
    size_t indKy;
    size_t indKz;
    size_t ind6D;
    forced_proc = mpi_whereIsM[2];
    local_m = mpi_whereIsM[3];
    // find out how many wave numbers will be forced
    if (mpi_my_col_rank == forced_proc){
        for(size_t i = 0; i < equation_forceKn; i++){
            indKx = equation_forceKxInd[i];
            indKy = equation_forceKyInd[i];
            indKz = equation_forceKzInd[i];
            for(size_t is = 0; is < array_local_size.ns; is++){
                ind6D = get_flat_c(is,0,local_m,indKx,indKy,indKz);
                if (cabs(h[ind6D]) > 1e-16){
                    rhs[ind6D] += equation_forcingCoef/conj(h[ind6D]);
                }
            }
        }
    }
};