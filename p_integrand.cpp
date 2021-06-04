/* p_integrand.c */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <boost/math/special_functions/gamma.hpp>
extern "C"
{
/* Survival function for R^phi*W */
int RW_marginal_C(double *xval, double phi, double gamma, int n_xval, double *result){
    double tmp2 = pow(gamma/2, phi)/boost::math::tgamma(0.5);
    double tmp1, tmp0, a;
    a = 0.5-phi;
    
    for(int i=0; i<n_xval; i++){
        tmp1 = gamma/(2*pow(xval[i],1/phi));
        tmp0 = tmp2/(a*xval[i]);
        result[i] = boost::math::gamma_p(0.5L,tmp1) + boost::math::tgamma((long double)(a+1),tmp1)*tmp0-pow(tmp1,a)*exp(-tmp1)*tmp0;
    }
    return 1;
}

/* Marginal distribution function for R^phi*W + epsilon */
int pRW_me_interp_C(double *xval, double *xp, double *surv_p, double tau_sqd, double phi, double gamma, int n_xval, int n_grid, double *result){
    bool tau_bool = (tau_sqd > 0.05);
    double tp[n_grid];
    double integrand_p[n_grid];
    double tmp, tmp_res; /* temporary constant */
    double tmp_sum = 0; /* temporary trapesoid sum */
    double sd = sqrt(tau_sqd);
    double sd_const = sqrt(2)*sd;
    double sd_const_pi =sqrt(2*M_PI)*sd;
    int i,j, tmp_int; /* iterative constants */

    for (i = 0; i < n_xval; i++) {
        if(tau_bool & (xval[i]<820)){
            /* Calculate integrand on a grid */
            for(j=0; j<n_grid;j++){
                tmp = xval[i]-xp[j];
                tp[j] = tmp;
                integrand_p[j] = exp(-tmp*tmp/(2*tau_sqd)) * surv_p[j];
            }
            
            /* Numerical integral using the trapesoid method */
            for(j=0; j<(n_grid-1);j++){
                tmp_sum+= (tp[j+1]-tp[j])*(integrand_p[j] + integrand_p[j+1])/2;
            }
            tmp_res = 0.5*erfc(-xval[i]/sd_const)-tmp_sum/sd_const_pi;
            tmp_sum = 0;
            
            /* CDF value must be greater than 0 */
            if(tmp_res < 0){
                tmp_res = 0;
            }
            result[i] = tmp_res;
        }
        else{
            tmp_int = RW_marginal_C(&xval[i], phi, gamma, 1, &tmp_res);
            result[i] = 1-tmp_res;
        }
    }
    
    return 1;
}

/* Get the quantile range for certain probability levels */
int find_xrange_pRW_me_C(double min_p, double max_p, double min_x, double max_x, double *xp, double *surv_p, double tau_sqd, double phi, double gamma, int n_grid, double *x_range){
    if (min_x >= max_x){
        printf("Initial value of mix_x must be smaller than max_x.\n");
        exit(EXIT_FAILURE);
    }
    
    /* First the min */
    double p_min_x;
    int tmp_int;
    tmp_int = pRW_me_interp_C(&min_x, xp, surv_p, tau_sqd, phi, gamma, 1, n_grid, &p_min_x);
    while (p_min_x > min_p){
        min_x = min_x-40/phi;
        tmp_int = pRW_me_interp_C(&min_x, xp, surv_p, tau_sqd, phi, gamma, 1, n_grid, &p_min_x);
    }
        
    x_range[0] = min_x;
    
    /* Now the max */
    double p_max_x;
    tmp_int = pRW_me_interp_C(&max_x, xp, surv_p, tau_sqd, phi, gamma, 1, n_grid, &p_max_x);
    while (p_max_x < max_p){
        max_x = max_x*2; /* Upper will set to 20 initially */
        tmp_int = pRW_me_interp_C(&max_x, xp, surv_p, tau_sqd, phi, gamma, 1, n_grid, &p_max_x);
    }
        
    x_range[1] = max_x;
    return 1;
}



/* Density function for R^phi*W */
int RW_density_C(double *xval, double phi, double gamma, int n_xval, double *result){
    double tmp2 = pow(gamma/2, phi)/boost::math::tgamma(0.5);
    double tmp1, tmp0, a;
    a = 0.5-phi;
    
    for(int i=0; i<n_xval; i++){
        tmp1 = gamma/(2*pow(xval[i],1/phi));
        tmp0 = tmp2/(a*pow(xval[i],2));
        result[i] = (boost::math::tgamma((long double)(a+1),tmp1)-pow(tmp1,a)*exp(-tmp1))*tmp0;
    }
    return 1;
}

/* Marginal density function for R^phi*W + epsilon */
int dRW_me_interp_C(double *xval, double *xp, double *den_p, double tau_sqd, double phi, double gamma, int n_xval, int n_grid, double *result){
    double thresh_large = 820;
    if(tau_sqd < 1) {
        thresh_large = 50;
    }
    bool tau_bool = (tau_sqd > 0.05);
    
    double tp[n_grid];
    double integrand_p[n_grid];
    double tmp, tmp_res; /* temporary constant */
    double tmp_sum = 0; /* temporary trapesoid sum */
    double sd = sqrt(tau_sqd);
    double sd_const_pi =sqrt(2*M_PI)*sd;
    int i,j, tmp_int; /* iterative constants */

    for (i = 0; i < n_xval; i++) {
        if(tau_bool & (xval[i]<thresh_large)){
            /* Calculate integrand on a grid */
            for(j=0; j<n_grid;j++){
                tmp = xval[i]-xp[j];
                tp[j] = tmp;
                integrand_p[j] = exp(-tmp*tmp/(2*tau_sqd)) * den_p[j];
            }
            
            /* Numerical integral using the trapesoid method */
            for(j=0; j<(n_grid-1);j++){
                tmp_sum+= (tp[j+1]-tp[j])*(integrand_p[j] + integrand_p[j+1])/2;
            }
            tmp_res = tmp_sum/sd_const_pi;
            tmp_sum = 0;
            result[i] = tmp_res;
        }else if((tau_bool & (xval[i]>=thresh_large))|(!tau_bool & (xval[i]>0))){
            tmp_int = RW_density_C(&xval[i], phi, gamma, 1, &tmp_res);
            result[i] = tmp_res;
        }else{
            result[i] = 0;
        }
    }
    
    return 1;
}



}

