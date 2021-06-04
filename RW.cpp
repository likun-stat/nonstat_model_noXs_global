/* p_integrand.c */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <boost/math/special_functions/gamma.hpp>
extern "C"

/* Survival function for R^phi*W */
int RW_marginal_C(double *xval, double phi, double gamma, int n_xval, double *result){
    double tmp2 = pow(gamma/2, phi)/boost::math::tgamma(0.5);
    double tmp1, tmp0, a;
    a = 0.5-phi;
    
    for(int i=0; i<n_xval; i++){
        tmp1 = gamma/(2*pow(xval[i],1/phi));
        tmp0 = tmp2/(a*xval[i]);
        result[i] = boost::math::gamma_p(0.5L,tmp1) + boost::math::tgamma((long double)(a+1),tmp1)*tmp0-pow(tmp1,a)*exp(-tmp1)*tmp0; /*boost::math::gamma_p(0.5L,tmp1)+ tmp2*boost::math::tgamma(0.5-phi,tmp1)/xval[i];*/
    }
    return 1;
}
