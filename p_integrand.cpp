/* p_integrand.cpp */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <algorithm>

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


/* Transform to uniform scales from RW mixtures */
int RW_me_2_unifs(double *X, double *xp, double *Surv, double tau_sqd, double *phi, double gamma,
                    int n_s, int n_grid, int n_t, double *unifs){
    int tmp_int, X_lookup, Surv_lookup;
    for (int i = 0; i<n_s;i++){
        X_lookup = i*n_t;
        Surv_lookup = i*n_grid;
        tmp_int = pRW_me_interp_C(X+X_lookup, xp, Surv+Surv_lookup, tau_sqd, phi[i], gamma, n_t, n_grid, unifs+X_lookup);
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

/* PCHIP 1-D monotonic cubic interpolation */
using boost::math::interpolators::cubic_hermite;
int pchip(double *x_input, double *y_input, double *x_pred, double *y_pred_res, int n, int n_pred){
    std::vector<double> x(n);
    memcpy(&x[0], x_input, n*sizeof(double));
    std::vector<double> y(n);
    memcpy(&y[0], y_input, n*sizeof(double));
   
    double left_endpoint_derivative = std::numeric_limits<double>::quiet_NaN();
    double right_endpoint_derivative = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> s(x.size(), std::numeric_limits<double>::quiet_NaN());
    if (isnan(left_endpoint_derivative)){
        // O(h) finite difference derivative:
        // This, I believe, is the only derivative guaranteed to be monotonic:
        s[0] = (y[1]-y[0])/(x[1]-x[0]);
    }else{
        s[0] = left_endpoint_derivative;
    }

    for (int k = 1; k < n-1; ++k) {
        double hkm1 = x[k] - x[k-1];
        double dkm1 = (y[k] - y[k-1])/hkm1;

        double hk = x[k+1] - x[k];
        double dk = (y[k+1] - y[k])/hk;
        double w1 = 2*hk + hkm1;
        double w2 = hk + 2*hkm1;
        if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0){
            s[k] = 0;
        }else{
            s[k] = (w1+w2)/(w1/dkm1 + w2/dk);
        }

    }
    // Quadratic extrapolation at the other end:
    if (isnan(right_endpoint_derivative)){
                s[n-1] = (y[n-1]-y[n-2])/(x[n-1] - x[n-2]);
    }else{
                s[n-1] = right_endpoint_derivative;
    }
    
    auto spline = cubic_hermite<std::vector<double>>(std::move(x), std::move(y), std::move(s));

    for(int j=0; j<n_pred; j++){
        y_pred_res[j] = spline(x_pred[j]);
    }
    return 1;
}


/* Quantile inverse function for R^phi*W + epsilon*/
int qRW_me_interp(double *p, double *xp, double *surv_p, double tau_sqd, double phi, double gamma,
                  int n_p, int n_grid, double *result,
                  double *cdf_vals, double *x_vals, int n_x=400, double lower=5, double upper=20
                ){
    int tmp_int, i, zeros;
    double min_p, max_p, delta, tmp_double;
    bool large_delta_large_x = false;
    
    // (1) When phi is varying over space, we need to calculate one quantile for each phi value.
    // Given more accuarte [lower,upper] values for the single p, we can decrease n_x.
    if(n_p<2){ n_x=floor(100*(p[0]+0.1));}
    double x_vals1[n_x+150];
    double cdf_vals1[n_x+150];
    
    // (2) Generate x_vals and cdf_vals to interpolate
    if(isnan(x_vals[0])){
        double x_range[2];
        max_p = *std::max_element(p, p+n_p);
        min_p = *std::min_element(p, p+n_p);
        tmp_int = find_xrange_pRW_me_C(min_p, max_p, lower, upper, xp, surv_p, tau_sqd,
                                       phi, gamma, n_grid, x_range);
        if (isinf(x_range[1])){   // Upper is set to 20 initially
            x_range[1] = std::pow(10,20);
            large_delta_large_x = true;
        }
        if (x_range[0]<=0){
            delta = (0.0001 - x_range[0]) / 149;
            for(i=0; i<150; i++){
                x_vals1[i] = x_range[0] + delta * i;
            }
            delta = (log(x_range[1])-log(0.0001001)) / (n_x-1);
            for(i=0; i<n_x; i++){
                x_vals1[i+150] = exp(log(0.0001001) + delta * i);
            }
            n_x=150+n_x;
        }else{
            delta = (log(x_range[1])-log(x_range[0])) / (n_x-1);
            for(i=0; i<n_x; i++){
                x_vals1[i] = exp(log(x_range[0]) + delta * i);
            }
        }
        tmp_int = pRW_me_interp_C(x_vals1, xp, surv_p, tau_sqd, phi, gamma, n_x, n_grid, cdf_vals1);
    }else{
        double x_vals1[n_x];
        double cdf_vals1[n_x];
        memcpy(x_vals1, x_vals, n_x*sizeof(double));
        if(isnan(cdf_vals[0])){
            tmp_int = pRW_me_interp_C(x_vals1, xp, surv_p, tau_sqd, phi, gamma, n_x, n_grid, cdf_vals1);
        }else{
            memcpy(cdf_vals1, cdf_vals, n_x*sizeof(double));
        }
    }
    
    // (3) Obtain the quantile level using the interpolated function
    if(!large_delta_large_x){
        zeros = 0;
        tmp_double = cdf_vals1[zeros];
        while (tmp_double==0) {
            tmp_double = cdf_vals1[++zeros];
        }
        tmp_int = pchip(cdf_vals1+zeros, x_vals1+zeros, p, result, n_x-zeros, n_p);
    }else{
        int exceeds = 0;
        for (i=0;i<n_p;i++){
            if(p[i]>cdf_vals1[n_x-1]) exceeds+=1;
        }
        tmp_int = pchip(cdf_vals1, x_vals1, p, result, n_x, n_p-exceeds);
        for (i=n_p-exceeds;i<n_p;i++){
            result[i] = x_vals[n_x-1];
        }
    }
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

int density_interp_grid(double *xp, double *phi, double gamma, int n_phi, int n_grid, double *Den, double *Surv){
    int counter = 0;
    int i,j, tmp_int;
    double tmp_surv, tmp_den;
    double tmp2, tmp1, tmp0, a, tmp_incomp, tmp_phi, tmp_phi_inv, tmp_xp;
    double gamma_half = gamma/2;
    
    for(i=0; i<n_phi; i++){
        tmp_phi = phi[i];
        a = 0.5-tmp_phi;
        tmp2 = std::pow(gamma_half, tmp_phi)/(a*sqrt(M_PI));
        tmp_phi_inv = 1/tmp_phi;
        for(j=0; j<n_grid; j++){
            tmp_xp = xp[j];
            tmp1 = gamma_half/std::pow(tmp_xp,tmp_phi_inv);
            tmp0 = tmp2/tmp_xp;
            tmp_incomp = (boost::math::tgamma((long double)(a+1),tmp1)-std::pow(tmp1,a)*exp(-tmp1))*tmp0;
            Surv[counter] = boost::math::gamma_p(0.5L,tmp1) + tmp_incomp;
            Den[counter++] = tmp_incomp/tmp_xp;
        }
    }
    return 1;
}

double dgev_C(double y, double loc, double scale, double shape, bool log_out){
    double t = std::pow(1+shape*((y-loc)/scale), -1/shape);
    double result;
    if(log_out){
        result = -log(scale)+(shape+1)*log(t)-t;
    }else{
        result = std::pow(t, shape+1)*exp(-t)/scale;
    }
    return result;
}

double dnorm_C(double y, double mean, double sd, bool log_out){
    double t=(y-mean)/sd;
    double result;
    if(log_out){
        result = -0.5*log(2*M_PI)-log(sd)-0.5*t*t;
    }else{
        result = exp(-0.5*t*t)/(sqrt(2*M_PI)*sd);
    }
    return result;
}

/* Thresh_X and Thresh_X_above are required */
/* xp, den_p and surv_p are required */
/* Calculate column_wise in C order OR one time  */
double marg_transform_data_mixture_me_likelihood_C(double *Y, double *X, double *X_s, bool *cen, bool *cen_above,
                                                double *Loc, double *Scale, double *Shape,
                                                double tau_sqd, double *phi, double gamma,
                                                double *xp, double *Den, int n_s, int n_grid){
    double sd = sqrt(tau_sqd);
    double sd_const = sqrt(2)*sd;
    double ll=0;
    double RW_den;
    int i, tmp_int, Den_lookup;
    
    for (i=0; i<n_s; i++){
        if(cen[i]){
            ll += log(0.5*erfc(-(X[i]-X_s[i])/sd_const));
        }else if(cen_above[i]){
            ll += log(1-0.5*erfc(-(X[i]-X_s[i])/sd_const));
        }else{
            Den_lookup = i*n_grid;
            tmp_int = dRW_me_interp_C(&X[i], xp, &Den[Den_lookup], tau_sqd, phi[i], gamma, 1, n_grid, &RW_den);
            ll += dnorm_C(X[i], X_s[i], sd, true)+dgev_C(Y[i], Loc[i], Scale[i], Shape[i], true)-log(RW_den);
        }
    }
    return ll;
}

/* Calculate row_wise in F order OR one location */
double marg_transform_data_mixture_me_likelihood_F(double *Y, double *X, double *X_s, bool *cen, bool *cen_above,
                                                double *Loc, double *Scale, double *Shape,
                                                double tau_sqd, double phi, double gamma,
                                                double *xp, double *den_p, int n_t, int n_grid){
    double sd = sqrt(tau_sqd);
    double sd_const = sqrt(2)*sd;
    double ll=0;
    double RW_den;
    int i, tmp_int;
    
    for (i=0; i<n_t; i++){
        if(cen[i]){
            ll += log(0.5*erfc(-(X[i]-X_s[i])/sd_const));
        }else if(cen_above[i]){
            ll += log(1-0.5*erfc(-(X[i]-X_s[i])/sd_const));
        }else{
            tmp_int = dRW_me_interp_C(&X[i], xp, den_p, tau_sqd, phi, gamma, 1, n_grid, &RW_den);
            ll += dnorm_C(X[i], X_s[i], sd, true)+dgev_C(Y[i], Loc[i], Scale[i], Shape[i], true)-log(RW_den);
        }
    }
    return ll;
}

/* Calculate all locations and all times */
double marg_transform_data_mixture_me_likelihood_global(double *Y, double *X, double *X_s, bool *cen, bool *cen_above,
                                                double *Loc, double *Scale, double *Shape,
                                                double tau_sqd, double *phi, double gamma,
                                                double *xp, double *Den, int n_s, int n_t, int n_grid){
    
    double ll=0;
    int site, tmp_int, Den_lookup, X_lookup;
    
    for (site=0; site<n_s; site++){
        X_lookup = site*n_t;
        Den_lookup = site*n_grid;
        ll += marg_transform_data_mixture_me_likelihood_F(&Y[X_lookup], &X[X_lookup], &X_s[X_lookup], &cen[X_lookup], &cen_above[X_lookup], &Loc[X_lookup], &Scale[X_lookup], &Shape[X_lookup], tau_sqd, phi[site], gamma, xp, &Den[Den_lookup], n_t, n_grid);
    }
    return ll;
}

/* Calculate thresh_X and thresh_X_above for an array of phi values */
// First way: iterate through each prob at each location (Not efficient)
int Thresh_X_try(double *phi, double *xp, double *Surv, double prob_below, double prob_above, double tau_sqd, double gamma,
                    double below_approx, double above_approx,
                    int n_phi, int n_grid, double *Thresh_X, double *Thresh_X_above){
    int Surv_lookup, tmp_int;
    int n_x=100;
    double cdf_vals[n_x], x_vals[n_x];
    double tmp_res;
    cdf_vals[0] = std::numeric_limits<double>::quiet_NaN(); x_vals[0] = std::numeric_limits<double>::quiet_NaN();
    
    for (int i = 0; i<n_phi;i++){
        Surv_lookup = i*n_grid;
        tmp_int = qRW_me_interp(&prob_below, xp, &Surv[Surv_lookup], tau_sqd, phi[i], gamma,
                          1, n_grid, &tmp_res,
                          cdf_vals, x_vals, n_x, below_approx-10, below_approx+10);
        Thresh_X[i] = tmp_res;
        tmp_int = qRW_me_interp(&prob_above, xp, &Surv[Surv_lookup], tau_sqd, phi[i], gamma,
                          1, n_grid, &tmp_res,
                          cdf_vals, x_vals, n_x, above_approx-100, above_approx+200);
        Thresh_X_above[i] = tmp_res;
    }
    
    return 1;
}

// Second way: calculate probs altogether at each location
int X_update(double *phi, double *xp, double *Surv, double *probs, double tau_sqd, double gamma,
                    double below_approx, double above_approx,
                    int n_phi, int n_grid, int n_probs, double *X_res){
    int Surv_lookup, tmp_int, X_lookup;
    int n_x=100;
    double cdf_vals[n_x], x_vals[n_x];
    cdf_vals[0] = std::numeric_limits<double>::quiet_NaN(); x_vals[0] = std::numeric_limits<double>::quiet_NaN();
    
    for (int i = 0; i<n_phi;i++){
        Surv_lookup = i*n_grid;
        X_lookup = i*n_probs;
        tmp_int = qRW_me_interp(probs, xp, &Surv[Surv_lookup], tau_sqd, phi[i], gamma,
                          n_probs, n_grid, X_res + X_lookup,
                          cdf_vals, x_vals, n_x, below_approx-10, above_approx+200);
    }
    
    return 1;
}

int unifs_2_RW_me(double *unifs, double *xp, double *Surv, double tau_sqd, double *phi, double gamma,
                    double below_approx, double above_approx,
                    int n_phi, int n_grid, int n_probs_each_loc, double *X_res){
    int Surv_lookup, tmp_int, X_lookup;
    int n_x=100;
    double cdf_vals[n_x], x_vals[n_x];
    cdf_vals[0] = std::numeric_limits<double>::quiet_NaN(); x_vals[0] = std::numeric_limits<double>::quiet_NaN();
    
    for (int i = 0; i<n_phi;i++){
        Surv_lookup = i*n_grid;
        X_lookup = i*n_probs_each_loc;
        tmp_int = qRW_me_interp(unifs + X_lookup, xp, &Surv[Surv_lookup], tau_sqd, phi[i], gamma,
                          n_probs_each_loc, n_grid, X_res + X_lookup,
                          cdf_vals, x_vals, n_x, below_approx-10, above_approx+200);
    }
    
    return 1;
}

void print_c(double *Y, int n_grid){
    printf("%4.2f %4.2f\n",*Y,*(Y+n_grid-1));
}

double print_Vec(double *Y, int n_grid, int n_s){
    int Den_lookup;
    for(int i=0; i<n_s; i++){
        Den_lookup = i*n_grid;
        print_c(&Y[Den_lookup], n_grid);
    }
    
    return Y[0];
}


}
