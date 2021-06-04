#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 13:22:31 2021

Model Simulation & Grid Interpolation

@authors: Likun Zhang & Mark Risser
"""
import numpy as np
import sys
from scipy.stats import norm
from scipy.stats import uniform
import scipy.special as sc
import mpmath
import scipy.integrate as si
import scipy.interpolate as interp
import scipy.optimize as optim
from scipy.stats import genextreme




## integration.cpp
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                        Import C++ function library
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
## i.e., RW_marginal, pRW_me_interp, find_xrange_pRW_me
##
import os, ctypes

# g++ -std=c++11 -shared -fPIC -o p_integrand.so p_integrand.cpp
lib = ctypes.CDLL(os.path.abspath('./nonstat_model_noXs/p_integrand.so'))
i_and_o_type = np.ctypeslib.ndpointer(ndim=1, dtype=np.float64)
grid_type  = np.ctypeslib.ndpointer(ndim=1, dtype=np.float64)

lib.pRW_me_interp_C.restype = ctypes.c_int
lib.pRW_me_interp_C.argtypes = (i_and_o_type, grid_type, grid_type,
                      ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int,
                      i_and_o_type)

lib.RW_marginal_C.restype = ctypes.c_int
lib.RW_marginal_C.argtypes = (i_and_o_type, 
                      ctypes.c_double, ctypes.c_double, ctypes.c_int,
                      i_and_o_type)

lib.find_xrange_pRW_me_C.restype = ctypes.c_int
lib.find_xrange_pRW_me_C.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                     grid_type, grid_type, ctypes.c_double, ctypes.c_double, ctypes.c_double, 
                                     ctypes.c_int, i_and_o_type)

lib.RW_density_C.restype = ctypes.c_int
lib.RW_density_C.argtypes = (i_and_o_type, 
                      ctypes.c_double, ctypes.c_double, ctypes.c_int,
                      i_and_o_type)

lib.dRW_me_interp_C.restype = ctypes.c_int
lib.dRW_me_interp_C.argtypes = (i_and_o_type, grid_type, grid_type,
                      ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int,
                      i_and_o_type)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                       Generate Levy random samples
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
## i.e., Stable variables with alpha=1/2
##
def rlevy(n, m = 0, s = 1):
  if np.any(s < 0):
    sys.exit("s must be positive")
  return s/norm.ppf(uniform.rvs(0,1,n)/2)**2 + m

## The density for R^phi in which R is levy distributed
def dR_power_phi(x, phi, m=0, s=1, log=False):
    x_phi = x**(1/phi)
    if np.any(x_phi <= m):
        sys.exit("some x**phi <= m")
    if np.any(s <= 0): 
        sys.exit("s must be positive")
    tmp = np.sum(np.log(s/(2 * np.pi))/2 - 3 * np.log(x_phi - m)/2 - s/(2 * (x_phi - 
        m)) + (1/phi-1)*np.log(x)-np.log(phi))
    if np.invert(log): 
        tmp = np.exp(tmp)
    return tmp
        
## -------------------------------------------------------------------------- ##






## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##        Calculate unregularized upper incomplete gamma function
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
## The negative a values are allowed
##
def gammaincc_unregulized(a,x):
    if(isinstance(x, (int, np.int64, float))): x=np.array([x])
    if x.any()<0: sys.exit("x must be positive")
    if a>0: 
        return sc.gamma(a)*sc.gammaincc(a,x)
    elif a<0: 
        return gammaincc_unregulized(a+1,x)/a-(x**a)*np.exp(-x)/a
    else:
        return mpmath.gammainc(0,x)

## Compare with mpmath.gammainc
## gammaincc_unregulized is more efficient
# import time
# 
# start_time = time.time()
# gammaincc_unregulized(-3.62,5)
# time.time() - start_time

# start_time = time.time()
# mpmath.gammainc(-3.62,5)
# time.time() - start_time
## -------------------------------------------------------------------------- ##





## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##         Calculate the exact marginal survival function for R^phi*W
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
##
def RW_marginal_uni(x,phi,gamma,survival = True):
    tmp1 = gamma/(2*(x**(1/phi)))
    tmp2 = (gamma/2)**phi/sc.gamma(0.5)
    res = sc.gammainc(0.5,tmp1) + tmp2*gammaincc_unregulized(0.5-phi,tmp1)/x
    if survival:
        return res
    else:
        return 1-res

RW_marginal = np.vectorize(RW_marginal_uni)

    
def RW_marginal_asymp(x,phi,gamma):
    if phi<0.5:
        moment = ((2*gamma)**phi)*sc.gamma(1-2*phi)/sc.gamma(1-phi)
        return moment/x
    elif phi>0.5:
        return np.sqrt(2*gamma/np.pi)*(x**(-1/(2*phi)))/(1-1/(2*phi))
    else:
        return np.sqrt(2*gamma/np.pi)*np.log(x)/x

def RW_quantile_asymp(p,phi,gamma):
    if phi<0.5:
        moment = ((2*gamma)**phi)*sc.gamma(1-2*phi)/sc.gamma(1-phi)
        return moment/(1-p)
    elif phi>0.5:
        return (np.sqrt(2*gamma/np.pi)/(1-1/(2*phi))/(1-p))**(2*phi)
    else:
        tmp = (1-p)/np.sqrt(2*gamma/np.pi)
        return tmp/sc.lambertw(tmp)

# # Compare the exact and asymptotic CDF
# gamma = 1.2; x =10; phi=0.3

# import matplotlib.pyplot as plt
# axes = plt.gca()
# axes.set_ylim([0,0.125])
# X_vals = np.linspace(100,1500,num=200)
# P_vals = RW_marginal(X_vals,phi,gamma)
# P_asym = RW_marginal_asymp(X_vals,phi,gamma)
# plt.plot(X_vals, P_vals, 'b')
# plt.plot(X_vals, P_asym, 'r',linestyle='--');
## -------------------------------------------------------------------------- ##






## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##       Calculate the marginal survival function for R^phi*W + epsilon
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
# ----------------  1. Define integrand in Python: exact form ---------------- #
def mix_distn_integrand(t, xval, phi, tmp1, tmp2, tau_sqd):
    diff = xval - t
    tmp3 = tmp1/(diff**(1/phi))
    res = sc.gammainc(0.5,tmp3) + tmp2*gammaincc_unregulized(0.5-phi,tmp3)/diff
    result = res * np.exp(-t**2/(2*tau_sqd))
    return result

def pRW_me_uni(xval, phi, gamma, tau_sqd):
    tmp1 = gamma/2
    tmp2 = ((gamma/2)**phi)/sc.gamma(0.5)
    sd = np.sqrt(tau_sqd)
    
    I_1 = si.quad(mix_distn_integrand, -np.inf, xval, args=(xval, phi, tmp1, tmp2, tau_sqd)) # 0.00147
    tmp = norm.cdf(xval, loc=0.0, scale=sd)-I_1[0]/np.sqrt(2*np.pi*tau_sqd)
    if tmp<0.999:
        return tmp
    else:
        return RW_marginal_uni(xval,phi,gamma,survival = False)

pRW_me = np.vectorize(pRW_me_uni)


# -----------  2. Define integrand in Python: linear interpolation ----------- #
# Actually BETTER than numerical integration because there are no singular values.
# We use the Trapezoidal rule.
## **** (0). Generate a GRIDDED set of values for P(RW>x) ****
def survival_interp_grid(phi, gamma, grid_size=800):
    xp_1 = np.linspace(0.000001, 200, grid_size, endpoint = False)
    xp_2 = np.linspace(200.5, 900, int(grid_size/4), endpoint = False)
    xp_3 = np.linspace(900.5, 100000, int(grid_size/10), endpoint = False)
    xp = np.concatenate((xp_1, xp_2, xp_3))
    
    xp = xp[::-1] # reverse order
    xp = np.ascontiguousarray(xp, np.float64) #C contiguous order: xp.flags['C_CONTIGUOUS']=True?
    n_xval = len(xp); surv_p = np.empty(n_xval)
    tmp_int = lib.RW_marginal_C(xp, phi, gamma, n_xval, surv_p)
    if tmp_int!=1: sys.exit('C implementaion failed.')
    # surv_p = RW_marginal(xp, phi, gamma)
    return (xp, surv_p)


## **** (1). Vectorize univariate function ****
def pRW_me_uni_interp(xval, xp, surv_p, tau_sqd):
    tp = xval-xp
    integrand_p = np.exp(-tp**2/(2*tau_sqd)) * surv_p
    sd = np.sqrt(tau_sqd)
    
    I_1 = sum(np.diff(tp)*(integrand_p[:-1] + integrand_p[1:])/2)  # 0.00036
    tmp = norm.cdf(xval, loc=0.0, scale=sd)-I_1/np.sqrt(2*np.pi*tau_sqd)
    return tmp
   

def pRW_me_interp_slower(xval, xp, surv_p, tau_sqd):
    return np.array([pRW_me_uni_interp(xval_i, xp, surv_p, tau_sqd) for xval_i in xval])


## **** (2). Broadcast matrices and vectorize columns ****
def pRW_me_interp_py(xval, xp, surv_p, tau_sqd, phi, gamma):
    if(isinstance(xval, (int, np.int64, float))): xval=np.array([xval], dtype='float64')
    tmp = np.zeros(xval.shape) # Store the results
    # Use the smooth process CDF values if tau_sqd<0.05
    if tau_sqd>0.05: 
        which = (xval<820) 
    else: 
        which = np.repeat(False, xval.shape)
    
    # Calculate for values that are less than 820
    if(np.sum(which)>0):
        xval_less = xval[which]
        tp = xval_less-xp[:,np.newaxis]
        integrand_p = np.exp(-tp**2/(2*tau_sqd)) * surv_p[:,np.newaxis]
        sd = np.sqrt(tau_sqd) 
        
        ncol = integrand_p.shape[1]
        I_1 = np.array([np.sum(np.diff(tp[:,index])*(integrand_p[:-1,index] + integrand_p[1:,index])/2) 
              for index in np.arange(ncol)])
        tmp_res = norm.cdf(xval_less, loc=0.0, scale=sd)-I_1/np.sqrt(2*np.pi*tau_sqd)
        # Numerical negative when xval is very small
        if(np.any(tmp_res<0)): tmp_res[tmp_res<0] = 0
        tmp[which] = tmp_res 
        
    
    # Calculate for values that are greater than 820
    if(xval.size-np.sum(which)>0):
        tmp[np.invert(which)] = RW_marginal(xval[np.invert(which)],phi,gamma,survival = False)
    return tmp


## **** (3). Use the C implementation ****
def pRW_me_interp(xval, xp, surv_p, tau_sqd, phi, gamma):
    if(isinstance(xval, (int, np.int64, float))): xval=np.array([xval], dtype='float64')
    n_xval = len(xval); n_grid = len(xp)
    result = np.zeros(n_xval) # Store the results
    tmp_int = lib.pRW_me_interp_C(xval, xp, surv_p, tau_sqd, phi, gamma, n_xval, n_grid, result)
    if tmp_int!=1: sys.exit('C implementaion failed.')
    return result
    

# -----------  3. Define integrand in Python: linear interpolation ----------- #
# The grid in the previous version depends on gamma. It's not ideal.
## **** (0). Generate a GRIDDED set of values for the integrand ****
## When phi > 1 and gamma < 0.01 (fixed?) or gamma > 120, the abnormality kicks in quicker.
## When phi=0.7 and gamma=1e-05, it works fine.
def survival_interp_grid1(phi, grid_size=1000):
    sp_1 = np.linspace(0.000001, 400, grid_size, endpoint = False)
    sp_2 = np.linspace(400.5, 1100, int(grid_size/4), endpoint = False)
    sp_3 = np.linspace(1100.5, 100000, int(grid_size/10), endpoint = False)
    sp = np.concatenate((sp_1, sp_2, sp_3))
    
    tmp = 1/(sp**(1/phi))
    surv_p = sc.gammainc(0.5,tmp) + gammaincc_unregulized(0.5-phi,tmp)/(sp*sc.gamma(0.5))
    return (sp, surv_p)


def pRW_me_interp1(xval, sp, surv_p, tau_sqd, phi, gamma):
    if(isinstance(xval, (int, np.int64, float))): xval=np.array([xval], dtype='float64')
    res = np.zeros(xval.size) # Store the results
    tmp1 = (gamma/2)**phi
    # If the asymp quantile level reaches 0.98, use the smooth distribution func.
    thresh = max(RW_quantile_asymp(0.98,phi,gamma),7.5)  # 7.5 is for gamma<0.0001
    # Use the smooth process CDF values if tau_sqd<0.05
    if tau_sqd>0.05: 
        which = (xval<thresh) 
    else: 
        which = np.repeat(False, xval.shape)
    
    # Calculate for values that are less than 820
    if(np.sum(which)>0):
        xval_less = xval[which]
        tp = xval_less-tmp1*sp[:,np.newaxis]
        integrand_p = np.exp(-tp**2/(2*tau_sqd)) * surv_p[:,np.newaxis]
        sd = np.sqrt(tau_sqd) 
        
        ncol = integrand_p.shape[1]
        I_1 = np.array([np.sum(np.diff(sp)*(integrand_p[:-1,index] + integrand_p[1:,index])/2) 
              for index in np.arange(ncol)])
        tmp_res = norm.cdf(xval_less, loc=0.0, scale=sd)-tmp1*I_1/np.sqrt(2*np.pi*tau_sqd)
        # Numerical negative when xval is very small
        if(np.any(tmp_res<0)): tmp_res[tmp_res<0] = 0
        res[which] = tmp_res 
        
    
    # Calculate for values that are greater than 820
    if(xval.size-np.sum(which)>0):
        res[np.invert(which)] = RW_marginal(xval[np.invert(which)],phi,gamma,survival = False)
    return res


# import matplotlib.pyplot as plt
# axes = plt.gca()
# # axes.set_ylim([0,0.125])
# X_vals = np.linspace(0.001,120,num=300)

# import time
# phi=0.7; gamma=1.2; tau_sqd = 10
# start_time = time.time()
# P_vals = RW_marginal(X_vals,phi, gamma, survival = False)
# time.time() - start_time

# start_time = time.time()
# P_mix = pRW_me(X_vals,phi,gamma,tau_sqd)
# time.time() - start_time

# grid = survival_interp_grid(phi, gamma)
# xp = grid[0]; surv_p = grid[1]
# start_time = time.time()
# P_interp_slower = pRW_me_interp_slower(X_vals, xp, surv_p, tau_sqd)
# time.time() - start_time

# start_time = time.time()
# P_interp_py = pRW_me_interp_py(X_vals, xp, surv_p, tau_sqd, phi, gamma)
# time.time() - start_time

# start_time = time.time()
# P_interp = pRW_me_interp(X_vals, xp, surv_p, tau_sqd, phi, gamma)
# time.time() - start_time

# grid = survival_interp_grid1(phi)
# sp = grid[0]; surv_p1 = grid[1]
# start_time = time.time()
# P_interp1 = pRW_me_interp1(X_vals, sp, surv_p1, tau_sqd, phi, gamma)
# time.time() - start_time


# fig, ax = plt.subplots()
# ax.plot(X_vals, P_vals, 'b', label="Smooth R^phi*W")
# ax.plot(X_vals, P_mix, 'r',linestyle='--', label="With nugget: numerical int")
# ax.plot(X_vals, P_interp_slower, 'g',linestyle=':', label="With nugget: Linear interp 1")
# ax.plot(X_vals, P_interp_py, 'm',linestyle=':', label="With nugget: Linear interp 2")
# ax.plot(X_vals, P_interp, 'y',linestyle='-.', label="With nugget: Linear interp 2 in C++")
# ax.plot(X_vals, P_interp1, 'c',linestyle='-.', label="With nugget: Linear interp w/o gamma")
# # ax.scatter(86.50499743, 0.9, c='red')
# # ax.scatter(1.11750005, 0.19, c='red')
# legend = ax.legend(loc = "lower right",shadow=True)
# plt.show()
## ----- Compared to 0.02579 secs for 1000 values when using pmixture_me() --- ##
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
#|  RW_marginal |    pRW_me    | pRW_me_interp_slower  |   pRW_me_interp_py   | pRW_me_interp
# ---------------------------------------------------------------------------------------------
#|   smooth RW  |  exact marg  |  interp w/ vectorize  |  interp w/ broadcast | interp in C++
# ---------------------------------------------------------------------------------------------
#| 0.02799 secs | 4.85512 secs |      0.25337 secs     |      0.05140 secs    |  0.01389 secs
# ---------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------- ##







## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##       Calculate the quantile inverse function for R^phi*W + epsilon
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
# --------------------  1. Use Shaby's interpolated method ------------------- #
# Improvement: Now we can calculate lower quantile levels that are negative.
def find_xrange_pRW_me(min_p, max_p, x_init, xp, surv_p, tau_sqd, phi, gamma):
    x_range = np.zeros(2)
    min_x = x_init[0]
    max_x = x_init[1]
    # if min_x <= 0 or min_p <= 0.15:
    #     sys.exit('This will only work for x > 0, which corresponds to p > 0.15.')
    if min_x >= max_x:
        sys.exit('x_init[0] must be smaller than x_init[1].')
    
    ## First the min
    p_min_x = pRW_me_interp_py(min_x, xp, surv_p, tau_sqd, phi, gamma)
    while p_min_x > min_p:
        min_x = min_x-40/phi
        p_min_x = pRW_me_interp_py(min_x, xp, surv_p, tau_sqd, phi, gamma)
    x_range[0] = min_x
    
    ## Now the max
    p_max_x = pRW_me_interp_py(max_x, xp, surv_p, tau_sqd, phi, gamma)
    while p_max_x < max_p:
        max_x = max_x*2 # Upper will set to 20 initially
        p_max_x = pRW_me_interp_py(max_x, xp, surv_p, tau_sqd, phi, gamma)
    x_range[1] = max_x
    return x_range

import sklearn
import sklearn.isotonic
# from sklearn.experimental import enable_hist_gradient_boosting  # noqa
# from sklearn.ensemble import HistGradientBoostingRegressor
# x_vals_tmp = x_vals.reshape(-1,1)
# cdf_gbdt = HistGradientBoostingRegressor(monotonic_cst=[1]).fit(x_vals_tmp, cdf_vals)
# cdf_vals_1 = cdf_gbdt.predict(x_vals_tmp)

def qRW_me_interp_py(p, xp, surv_p, tau_sqd, phi, gamma, 
           cdf_vals = np.nan, x_vals = np.nan, n_x=400, lower=5, upper=20):
    if(isinstance(p, (int, np.int64, float))): p=np.array([p])
    large_delta_large_x = False
    # Generate x_vals and cdf_vals to interpolate
    if np.any(np.isnan(x_vals)):
        x_range = find_xrange_pRW_me(np.min(p),np.max(p), np.array([lower,upper]), 
                                     xp, surv_p, tau_sqd, phi, gamma)
        if np.isinf(x_range[1]): # Upper is set to 20 initially
            x_range[1] = 10^20; large_delta_large_x = True
        if np.any(x_range<=0):
            x_vals = np.concatenate((np.linspace(x_range[0], 0.0001, num=150),
                           np.exp(np.linspace(np.log(0.0001001), np.log(x_range[1]), num=n_x))))
        else:
            x_vals = np.exp(np.linspace(np.log(x_range[0]), np.log(x_range[1]), num=n_x))
        cdf_vals = pRW_me_interp_py(x_vals, xp, surv_p, tau_sqd, phi, gamma)
    else:
        if np.any(np.isnan(cdf_vals)):
            cdf_vals = pRW_me_interp_py(x_vals, xp, surv_p, tau_sqd, phi, gamma)
    
    # Obtain the quantile level using the interpolated function
    if not large_delta_large_x:
        zeros = sum(cdf_vals==0)
        try:
            tck = interp.pchip(cdf_vals[zeros:], x_vals[zeros:]) # 1-D monotonic cubic interpolation.
        except ValueError:
            ir = sklearn.isotonic.IsotonicRegression(increasing=True)
            ir.fit(x_vals[zeros:], cdf_vals[zeros:])
            cdf_vals_1 = ir.predict(x_vals[zeros:])
            indices = np.where(np.diff(cdf_vals_1)==0)[0]+1
            tck = interp.pchip(np.delete(cdf_vals_1,indices), np.delete(x_vals[zeros:],indices))
        q_vals = tck(p)
    else:
        which = p>cdf_vals[-1]
        q_vals = np.repeat(np.nan, np.shape(p)[0])
        q_vals[which] = x_range[1]
        if np.any(~which):
            # tck = interp.interp1d(cdf_vals, x_vals, kind = 'cubic')
            tck = interp.pchip(cdf_vals, x_vals)
            q_vals[~which] = tck(p[~which])  
    return q_vals


# Same function implemented in C++
def qRW_me_interp(p, xp, surv_p, tau_sqd, phi, gamma, 
           cdf_vals = np.nan, x_vals = np.nan, n_x=400, lower=5, upper=20):
    if(isinstance(p, (int, np.int64, float))): p=np.array([p])
    large_delta_large_x = False
    # (1) When phi is varying over space, we need to calculate one quantile for each phi value.
    # Given more accuarte [lower,upper] values for the single p, we can decrease n_x.
    if(len(p)==1): n_x=np.int(100*(p+0.1))
    
    # (2) Generate x_vals and cdf_vals to interpolate
    if np.any(np.isnan(x_vals)):
        x_range = np.empty(2); n_grid = len(xp)
        tmp_int = lib.find_xrange_pRW_me_C(np.min(p), np.max(p), lower, upper, 
                                     xp, surv_p, tau_sqd, phi, gamma, n_grid, x_range)
        if tmp_int!=1: sys.exit('C implementaion failed.')
        if np.isinf(x_range[1]): # Upper is set to 20 initially
            x_range[1] = 10^20; large_delta_large_x = True
        if np.any(x_range<=0):
            x_vals = np.concatenate((np.linspace(x_range[0], 0.0001, num=150),
                           np.exp(np.linspace(np.log(0.0001001), np.log(x_range[1]), num=n_x))))
        else:
            x_vals = np.exp(np.linspace(np.log(x_range[0]), np.log(x_range[1]), num=n_x))
        cdf_vals = pRW_me_interp(x_vals, xp, surv_p, tau_sqd, phi, gamma)
    else:
        if np.any(np.isnan(cdf_vals)):
            cdf_vals = pRW_me_interp(x_vals, xp, surv_p, tau_sqd, phi, gamma)
    
    # (3) Obtain the quantile level using the interpolated function
    if not large_delta_large_x:
        zeros = sum(cdf_vals==0)
        try:
            tck = interp.pchip(cdf_vals[zeros:], x_vals[zeros:]) # 1-D monotonic cubic interpolation.
        except ValueError:
            ir = sklearn.isotonic.IsotonicRegression(increasing=True)
            ir.fit(x_vals[zeros:], cdf_vals[zeros:])
            cdf_vals_1 = ir.predict(x_vals[zeros:])
            indices = np.where(np.diff(cdf_vals_1)==0)[0]+1
            tck = interp.pchip(np.delete(cdf_vals_1,indices), np.delete(x_vals[zeros:],indices))
        q_vals = tck(p)
    else:
        which = p>cdf_vals[-1]
        q_vals = np.repeat(np.nan, np.shape(p)[0])
        q_vals[which] = x_range[1]
        if np.any(~which):
            # tck = interp.interp1d(cdf_vals, x_vals, kind = 'cubic')
            tck = interp.pchip(cdf_vals, x_vals)
            q_vals[~which] = tck(p[~which])  
    return q_vals




    
# ---------------------------  2. Use Scipy: Slower -------------------------- #
## Slow especially for phi>1.5
def diff(x, cdf_target, xp, surv_p, tau_sqd, phi, gamma):
    cdf_val = pRW_me_interp(x, xp, surv_p, tau_sqd, phi, gamma)
    return (cdf_val - cdf_target)**2


def qRW_me_optim(p, xp, surv_p, tau_sqd, phi, gamma):
    if(isinstance(p, (int, np.int64, float))): p=np.array([p])
    q_vals = np.zeros(p.shape)
    for idx,p_value in enumerate(p):
        res = optim.minimize(diff, 1.0, args=(p_value,xp, surv_p, tau_sqd, phi, gamma), 
                         method='Nelder-Mead', tol=1e-6)
        q_vals[idx] = res.x[0]   
    return q_vals


# import matplotlib.pyplot as plt
# phi=0.7; gamma=1.2; tau_sqd = 10
# grid = survival_interp_grid(phi, gamma)
# xp = grid[0]; surv_p = grid[1]
# P_vals = np.linspace(0.001,0.99,num=300)

# import time
# start_time = time.time()
# Q_interp_py = qRW_me_interp_py(P_vals, xp, surv_p, tau_sqd, phi, gamma) #0.02546 secs
# time.time() - start_time

# start_time = time.time()
# Q_interp = qRW_me_interp(P_vals, xp, surv_p, tau_sqd, phi, gamma) #0.00863 secs
# time.time() - start_time

# start_time = time.time()
# Q_optim = qRW_me_optim(P_vals, xp, surv_p, tau_sqd, phi, gamma) # 0.99205 secs
# time.time() - start_time

# fig, ax = plt.subplots()
# ax.plot(Q_interp_py, P_vals, 'b', label="Ben's interp method")
# ax.plot(Q_interp, P_vals, 'b', label="Ben's interp method in C++")
# ax.plot(Q_optim, P_vals, 'r',linestyle='--', label="Scipy's optim")
# legend = ax.legend(loc = "lower right",shadow=True)
# plt.show()  
## -------------------------------------------------------------------------- ##







## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##              Calculate the exact density function for R^phi*W
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
##
def RW_density_uni(x,phi,gamma, log=True):
    tmp1 = gamma/(2*(x**(1/phi)))
    tmp2 = (gamma/2)**phi/sc.gamma(0.5)
    res = tmp2*gammaincc_unregulized(0.5-phi,tmp1)/(x**2)
    if log:
        return np.log(res)
    else:
        return res

RW_density = np.vectorize(RW_density_uni)

def RW_density_asymp(x,phi,gamma):
    if phi<0.5:
        moment = ((2*gamma)**phi)*sc.gamma(1-2*phi)/sc.gamma(1-phi)
        return moment/(x**2)
    elif phi>0.5:
        return np.sqrt(2*gamma/np.pi)*(x**(-1/(2*phi)-1))/(2*phi-1)
    else:
        return np.sqrt(2*gamma/np.pi)*(np.log(x)-1)/(x**2)

# # Compare the exact and asymptotic CDF
# gamma = 1.2; x =10; phi=0.3

# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# X_vals = np.linspace(100,1500,num=200)
# # X_vals = np.linspace(150,350,num=200) # For phi=0.3
# P_vals = RW_density(X_vals,phi,gamma,log=False)
# P_asym = RW_density_asymp(X_vals,phi,gamma)
# ax.plot(X_vals, P_vals, 'b', label="R^phi*W density")
# ax.plot(X_vals, P_asym, 'r',linestyle='--', label="R^phi*W tail approx")
# legend = ax.legend(loc = "upper right",shadow=True)
# plt.title(label="phi=0.3") 
# plt.show()
## -------------------------------------------------------------------------- ##






## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##           Calculate the density function for R^phi*W + epsilon
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
# ----------------  1. Define integrand in Python: exact form ---------------- #
def mix_den_integrand(t, xval, phi, tmp1, tmp2, tau_sqd):
    diff = xval - t
    tmp3 = tmp1/(diff**(1/phi))
    res = tmp2*gammaincc_unregulized(0.5-phi,tmp3)/(diff**2)
    result = res * np.exp(-t**2/(2*tau_sqd))
    return result

def dRW_me_uni(xval, phi, gamma, tau_sqd):
    tmp1 = gamma/2
    tmp2 = ((gamma/2)**phi)/sc.gamma(0.5)
    denom = np.sqrt(2*np.pi*tau_sqd)
    
    I_1 = si.quad(mix_den_integrand, -np.inf, xval, args=(xval, phi, tmp1, tmp2, tau_sqd)) 
    tmp = I_1[0]/denom
    if tmp<0.00001 and xval>300:
        return RW_density_uni(xval,phi,gamma,log = False)
    else:
        return tmp

dRW_me = np.vectorize(dRW_me_uni)


# -----------  2. Define integrand in Python: linear interpolation ----------- #
# Actually BETTER than numerical integration because there are no singular values.
# We use the Trapezoidal rule.
## **** (0). Generate a GRIDDED set of values for f_RW(x) ****
def density_interp_grid(phi, gamma, grid_size=800):
    xp_1 = np.linspace(0.000001, 200, grid_size, endpoint = False)
    xp_2 = np.linspace(200.5, 900, int(grid_size/4), endpoint = False)
    xp_3 = np.linspace(900.5, 100000, int(grid_size/10), endpoint = False)
    xp = np.concatenate((xp_1, xp_2, xp_3))
    
    xp = xp[::-1] # reverse order
    xp = np.ascontiguousarray(xp, np.float64) #C contiguous order: xp.flags['C_CONTIGUOUS']=True?
    n_xval = len(xp); den_p = np.empty(n_xval); surv_p = np.empty(n_xval)
    tmp_int = lib.RW_density_C(xp, phi, gamma, n_xval, den_p) # density
    if tmp_int!=1: sys.exit('C implementaion failed.')
    tmp_int = lib.RW_marginal_C(xp, phi, gamma, n_xval, surv_p) #cdf
    if tmp_int!=1: sys.exit('C implementaion failed.')
    # den_p = RW_density(xp, phi, gamma, log=False)
    # surv_p = RW_marginal(xp, phi, gamma)
    return (xp, den_p, surv_p)


## **** (1). Vectorize univariate function ****
## Time: 0.05516 secs for 300 xvals.
def dRW_me_uni_interp(xval, xp, den_p, tau_sqd):
    tp = xval-xp
    integrand_p = np.exp(-tp**2/(2*tau_sqd)) * den_p
    denom = np.sqrt(2*np.pi*tau_sqd)
    
    I_1 = sum(np.diff(tp)*(integrand_p[:-1] + integrand_p[1:])/2)  # 0.00036
    tmp = I_1/denom
    return tmp
   

def dRW_me_interp_slower(xval, xp, den_p, tau_sqd):
    return np.array([dRW_me_uni_interp(xval_i, xp, den_p, tau_sqd) for xval_i in xval])

## **** (2). Broadcast matrices and vectorize columns ****
## Time: 0.05625 secs for 300 xvals. Faster 2x.
def dRW_me_interp_py(xval, xp, den_p, tau_sqd, phi, gamma, log =False):
    if(isinstance(xval, (int, np.int64, float))): xval=np.array([xval], dtype='float64')
    tmp = np.zeros(xval.size) # Store the results
    thresh_large=820
    if(tau_sqd<1): thresh_large = 50
    # Use the smooth process CDF values if tau_sqd<0.05
    if tau_sqd>0.05: 
        which = (xval<thresh_large) 
    else: 
        which = np.repeat(False, xval.shape)
    
    # Calculate for values that are less than 820
    if(np.sum(which)>0):
        xval_less = xval[which]
        tp = xval_less-xp[:,np.newaxis]
        integrand_p = np.exp(-tp**2/(2*tau_sqd)) * den_p[:,np.newaxis]
        denom = np.sqrt(2*np.pi*tau_sqd)
        
        ncol = integrand_p.shape[1]
        I_1 = np.array([np.sum(np.diff(tp[:,index])*(integrand_p[:-1,index] + integrand_p[1:,index])/2) 
              for index in np.arange(ncol)])
        tmp_res = I_1/denom
        # Numerical negative when xval is very small
        # if(np.any(tmp_res<0)): tmp_res[tmp_res<0] = 0
        tmp[which] = tmp_res 
        
    
    # Calculate for values that are greater than 820
    if(xval.size-np.sum(which)>0):
        tmp[np.invert(which) & (xval>0)] = RW_density(xval[np.invert(which) & (xval>0)],phi,gamma,log = False)
    
    if log:
        return np.log(tmp)
    else:
        return tmp
   


## **** (3). Use the C implementation ****
def dRW_me_interp(xval, xp, den_p, tau_sqd, phi, gamma, log=False):
    if(isinstance(xval, (int, np.int64, float))): xval=np.array([xval], dtype='float64')
    n_xval = len(xval); n_grid = len(xp)
    result = np.zeros(n_xval) # Store the results
    tmp_int = lib.dRW_me_interp_C(xval, xp, den_p, tau_sqd, phi, gamma, n_xval, n_grid, result)
    if tmp_int!=1: sys.exit('C implementaion failed.')
    if log:
        return np.log(result)
    else:
        return result


# import matplotlib.pyplot as plt
# axes = plt.gca()
# # axes.set_ylim([0,0.125])
# X_vals = np.linspace(0.001,100,num=300)

# import time
# phi=0.7; gamma=1.2; tau_sqd=10
# start_time = time.time()
# D_vals = RW_density(X_vals,phi, gamma, log = False)
# time.time() - start_time

# start_time = time.time()
# D_mix = dRW_me(X_vals,phi,gamma,tau_sqd)
# time.time() - start_time

# grid = density_interp_grid(phi, gamma)
# xp = grid[0]; den_p = grid[1]
# start_time = time.time()
# D_interp_slower = dRW_me_interp_slower(X_vals, xp, den_p, tau_sqd)
# time.time() - start_time

# start_time = time.time()
# D_interp_py = dRW_me_interp_py(X_vals, xp, den_p, tau_sqd, phi, gamma)
# time.time() - start_time

# start_time = time.time()
# D_interp = dRW_me_interp(X_vals, xp, den_p, tau_sqd, phi, gamma)
# time.time() - start_time

# fig, ax = plt.subplots()
# ax.plot(X_vals, D_vals, 'b', label="Smooth R^phi*W")
# ax.plot(X_vals, D_mix, 'r',linestyle='--', label="With nugget: numerical int")
# ax.plot(X_vals, D_interp_slower, 'g',linestyle=':', label="With nugget: Linear interp 1")
# ax.plot(X_vals, D_interp_py, 'y',linestyle='-.', label="With nugget: Linear interp 2")
# ax.plot(X_vals, D_interp, 'm',linestyle='-.', label="With nugget: Linear interp 2 in C++")
# legend = ax.legend(loc = "upper right",shadow=True)
# plt.title(label="phi=0.3") 
# plt.show()
# # ----- Compared to 0.01391 secs for 300 values when using pmixture_me() --- ##
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# |  RW_density  |    dRW_me    | dRW_me_interp_slower  |   dRW_me_interp_py   |  dRW_me_interp  |
# ------------------------------------------------------------------------------------------------
# |   smooth RW  |  exact marg  |  interp w/ vectorize  | interp w/ broadcast  |  interp in C++  |
# ------------------------------------------------------------------------------------------------
# | 0.00423 secs | 0.80787 secs |      0.05088 secs     |     0.01348 secs     |   0.00444 secs  |
# ------------------------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------------------- ##









## scalemix_utils.R

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                     Compute the Matern correlation function
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
## Input from a matrix of pairwise distances and a vector of parameters
##

def corr_fn(r, theta):
    if type(r).__module__!='numpy' or isinstance(r, np.float64):
      r = np.array(r)
    if np.any(r<0):
      sys.exit('Distance argument must be nonnegative.')
    r[r == 0] = 1e-10

    # range = theta[0]
    range = np.sqrt(theta[0])  # Mark's generic Matern
    nu = theta[1]
    part1 = 2 ** (1 - nu) / sc.gamma(nu)
    part2 = (r / range) ** nu
    part3 = sc.kv(nu, r / range)
    # part2 = (np.sqrt(2 * nu) * r / range) ** nu
    # part3 = sc.kv(nu, np.sqrt(2 * nu) * r / range)
    return part1 * part2 * part3
## -------------------------------------------------------------------------- ##




## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                                    For MVN
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.

## Computes A^{-1}x
##
def eig2inv_times_vector(V, d_inv, x):
  return V@np.diag(d_inv)@V.T@x

## Computes y=A^{-1}x via solving linear system Ay=x
from scipy.linalg import lapack
def inv_times_vector(A, x):
  inv = lapack.dposv(A,x)
  return inv


## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
##
## log(|A|)
##
def eig2logdet(d):
  return sum(np.log(d))


## Multivariate normal log density of R, where each column of 
## R iid N(mean,VDV'), where VDV' is the covariance matrix
## It essentially computes the log density of each column of R, then takes
## the sum.  Faster than looping over the columns, but not as transparent.
##                                     
## Ignore the coefficient: -p/2*log(2*pi)
##
def dmvn_eig(R, V, d_inv, mean=0):
  if len(R.shape)==1: 
    n_rep = 1 
  else: 
    n_rep = R.shape[1]
  res = -0.5*n_rep*eig2logdet(1/d_inv) - 0.5 * np.sum((R-mean) * eig2inv_times_vector(V, d_inv, R-mean))
  return res

def dmvn(R, Cor, mean=0, cholesky_inv = None):## cholesky_inv is the output of inv_times_vector()
  if len(R.shape)==1: 
    n_rep = 1 
  else: 
    n_rep = R.shape[1]
    
  if cholesky_inv is None: 
      inv = lapack.dposv(Cor, R-mean)
  else:
      sol = lapack.dpotrs(cholesky_inv[0],R-mean) #Solve Ax = b using factorization
      inv = (cholesky_inv[0],sol[0])
  logdet = 2*sum(np.log(np.diag(inv[0])))
  res = -0.5*n_rep*logdet - 0.5 * np.sum((R-mean) * inv[1])
  return res

## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
##
## Computes x'A^{-1}x
##
def eig2inv_quadform_vector(V, d_inv, x):
  cp = V@np.diag(d_inv)@V.T@x
  return sum(x*cp)

def inv_quadform_vector(Cor, x, cholesky_inv = None):
  if cholesky_inv is None:
      inv = lapack.dposv(Cor, x)
      cp = inv[1]
  else:
      sol = lapack.dpotrs(cholesky_inv[0],x) #Solve Ax = b using factorization
      cp = sol[0]
  return sum(x*cp)

## -------------------------------------------------------------------------- ##





## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                For generalized extreme value (GEV) distribution
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## 
## Negative shape parametrization in scipy.genextreme
## 

def dgev(yvals, Loc, Scale, Shape, log=False):
    if log:
        return genextreme.logpdf(yvals, c=-Shape, loc=Loc, scale=Scale)  # Opposite shape
    else:
        return genextreme.pdf(yvals, c=-Shape, loc=Loc, scale=Scale)  # Opposite shape

def pgev(yvals, Loc, Scale, Shape, log=False):
    if log:
        return genextreme.logcdf(yvals, c=-Shape, loc=Loc, scale=Scale)  # Opposite shape
    else:
        return genextreme.cdf(yvals, c=-Shape, loc=Loc, scale=Scale)  # Opposite shape

def qgev(p, Loc, Scale, Shape):
    if type(p).__module__!='numpy':
      p = np.array(p)  
    return genextreme.ppf(p, c=-Shape, loc=Loc, scale=Scale)  # Opposite shape


## -------------------------------------------------------------------------- ##






## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                    Transform Normal to Standard Pareto
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##
##
## i.e., Stable variables with alpha=1/2
##
def norm_to_Pareto(z):
    if(isinstance(z, (int, np.int64, float))): z=np.array([z])
    tmp = norm.cdf(z)
    if np.any(tmp==1): tmp[tmp==1]=1-1e-9
    return 1/(1-tmp)


def pareto_to_Norm(W):
    if(isinstance(W, (int, np.int64, float))): W=np.array([W])
    if np.any(W<1): sys.exit("W must be greater than 1")
    tmp = 1-1/W
    return norm.ppf(tmp)
        
## -------------------------------------------------------------------------- ##








## scalemix_likelihoods.R 

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                           Marginal transformations
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## 
## Transforms observations from a Gaussian scale mixture to a GPD, or vice versa
## 

def RW_me_2_gev(X, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape):
    unifs = pRW_me_interp(X, xp, surv_p, tau_sqd, phi, gamma)
    gevs = qgev(unifs, Loc=Loc, Scale=Scale, Shape=Shape)
    return gevs

def gev_2_RW_me(Y, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape):
    unifs = pgev(Y, Loc, Scale, Shape)
    scalemixes = qRW_me_interp(unifs, xp, surv_p, tau_sqd, phi, gamma)
    return scalemixes 

## After GEV params are updated, the 'cen' should be re-calculated.
def which_censored(Y, Loc, Scale, Shape, prob_below):
    unifs = pgev(Y, Loc, Scale, Shape)
    return unifs<prob_below

## Only calculate the un-censored elements
def X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape):
    X = np.empty(Y.shape)
    X[:] = np.nan
    
    if np.any(~cen & ~cen_above):
        X[~cen & ~cen_above] = gev_2_RW_me(Y[~cen & ~cen_above], xp, surv_p, tau_sqd, phi, gamma, 
                                    Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
    return X

## -------------------------------------------------------------------------- ##





## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                             Censored likelihood
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## The log likelihood of the data, where the data comes from a scale mixture
## of Gaussians, transformed to GEV (matrix/vector input)
##
## NOT ACTUALLY depending on X. X and cen need to be calculated in advance.
## 
##

def marg_transform_data_mixture_me_likelihood0(Y, X, X_s, cen, prob_below, Loc, Scale, Shape, 
                tau_sqd, phi, gamma, xp=np.nan, surv_p=np.nan, den_p=np.nan, thresh_X=np.nan):
  if np.any(np.isnan(xp)):
      grid = density_interp_grid(phi, gamma)
      xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
  if np.isnan(thresh_X):
     thresh_X = qRW_me_interp(prob_below, xp, surv_p, tau_sqd, phi, gamma)
  sd = np.sqrt(tau_sqd)  
  
  # ## Generate X_s
  # X_s = (R**phi)*norm_to_Pareto(Z)
  
  ## Initialize space to store the log-likelihoods for each observation:
  ll = np.empty(Y.shape)
  ll[:] = np.nan
  if np.any(cen):
     ll[cen] = norm.logcdf(thresh_X, loc=X_s[cen], scale=sd)
  
  if np.any(~cen):
     # # Sometimes pgev easily becomes 1, which causes the gev_2_scalemix to become nan
     # if np.any(np.isnan(X[~cen])):
     #     return -np.inf
     ll[~cen] = norm.logpdf(X[~cen], loc=X_s[~cen], scale=sd
               )+dgev(Y[~cen], Loc=Loc[~cen], Scale=Scale[~cen], Shape=Shape[~cen], log=True
               )-dRW_me_interp(X[~cen], xp, den_p, tau_sqd, phi, gamma, log =True)
     
  #which = np.isnan(ll)
  #if np.any(which):
  #   ll[which] = -np.inf  # Normal density larger order than marginal density of scalemix
  return np.sum(ll)

def marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, 
                prob_below, prob_above, Loc, Scale, Shape, tau_sqd, phi, gamma, 
                xp=np.nan, surv_p=np.nan, den_p=np.nan, thresh_X=np.nan, thresh_X_above=np.nan):
  if np.any(np.isnan(xp)):
      grid = density_interp_grid(phi, gamma)
      xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
  if np.isnan(thresh_X):
     thresh_X = qRW_me_interp(prob_below, xp, surv_p, tau_sqd, phi, gamma)
     thresh_X_above = qRW_me_interp(prob_above, xp, surv_p, tau_sqd, phi, gamma)
  sd = np.sqrt(tau_sqd)  
  
  # ## Generate X_s
  # X_s = (R**phi)*norm_to_Pareto(Z)
  
  ## Initialize space to store the log-likelihoods for each observation:
  ll = np.empty(Y.shape)
  ll[:] = np.nan
  if np.any(cen):
     ll[cen] = norm.logcdf(thresh_X, loc=X_s[cen], scale=sd)
  if np.any(cen_above):
     ll[cen_above] = norm.logsf(thresh_X_above, loc=X_s[cen_above], scale=sd)
  
  if np.any(~cen & ~cen_above):
     # # Sometimes pgev easily becomes 1, which causes the gev_2_scalemix to become nan
     # if np.any(np.isnan(X[~cen])):
     #     return -np.inf
     ll[~cen & ~cen_above] = norm.logpdf(X[~cen & ~cen_above], loc=X_s[~cen & ~cen_above], scale=sd
               )+dgev(Y[~cen & ~cen_above], Loc=Loc[~cen & ~cen_above], Scale=Scale[~cen & ~cen_above], Shape=Shape[~cen & ~cen_above], log=True
               )-dRW_me_interp(X[~cen & ~cen_above], xp, den_p, tau_sqd, phi, gamma, log =True)
     
  #which = np.isnan(ll)
  #if np.any(which):
  #   ll[which] = -np.inf  # Normal density larger order than marginal density of scalemix
  return np.sum(ll)


## Univariate version
def marg_transform_data_mixture_me_likelihood_uni(Y, X, X_s, cen, cen_above, 
                prob_below, prob_above, Loc, Scale, Shape, tau_sqd, phi, gamma, 
                xp=np.nan, surv_p=np.nan, den_p=np.nan, thresh_X=np.nan, thresh_X_above=np.nan):
  if np.any(np.isnan(xp)):
      grid = density_interp_grid(phi, gamma)
      xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
  if np.isnan(thresh_X):
     thresh_X = qRW_me_interp(prob_below, xp, surv_p, tau_sqd, phi, gamma)
     thresh_X_above = qRW_me_interp(prob_above, xp, surv_p, tau_sqd, phi, gamma)
  sd = np.sqrt(tau_sqd)  
  
  # ## Generate X_s
  # X_s = (R**phi)*norm_to_Pareto(Z)
  
  ll=np.array(np.nan)
  if cen:
     ll = norm.logcdf(thresh_X, loc=X_s, scale=sd)
  elif cen_above:
     ll = norm.logsf(thresh_X_above, loc=X_s, scale=sd)
  else:
     ll = norm.logpdf(X, loc=X_s, scale=sd
         )+dgev(Y, Loc=Loc, Scale=Scale, Shape=Shape, log=True
         )-dRW_me_interp(X, xp, den_p, tau_sqd, phi, gamma, log =True)
  #if np.isnan(ll):
  #   ll = -np.inf
  return ll

##
## -------------------------------------------------------------------------- ##




## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                           Full likelihood for phi
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians, where the mixing distribution comes from 
## dlevy.
##
##

def phi_update_mixture_me_likelihood(data, params, R, Z, cen, cen_above, 
                    prob_below, prob_above, Loc, Scale, Shape, tau_sqd, gamma):
  Y = data
  phi = params
  if phi < 0:
      return -np.inf
  
  grid = density_interp_grid(phi, gamma)
  xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  X_s = (R**phi)*norm_to_Pareto(Z)
  
  # ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, prob_below, Loc, Scale, Shape, 
  #                       tau_sqd, phi, gamma, xp, surv_p, den_p) + dR_power_phi(R,phi,m=0,s=gamma,log=True)

  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p)

  return ll
                                                                             
##
## -------------------------------------------------------------------------- ##




## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                           Full likelihood for tau
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians.
##
##

def tau_update_mixture_me_likelihood(data, params, X_s, cen, cen_above, 
                    prob_below, prob_above, Loc, Scale, Shape, 
                    phi, gamma, xp, surv_p, den_p):
  Y = data
  tau_sqd = params
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p)

  return ll

def phi_tau_update_mixture_me_likelihood(data, params, R, Z, cen, cen_above, 
                    prob_below, prob_above, Loc, Scale, Shape, gamma):
  Y = data
  phi = params[0]
  tau_sqd = params[1]
  if phi < 0:
      return -np.inf
  
  grid = density_interp_grid(phi, gamma)
  xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  X_s = (R**phi)*norm_to_Pareto(Z)
  # ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, prob_below, Loc, Scale, Shape, 
  #                       tau_sqd, phi, gamma, xp, surv_p, den_p) + dR_power_phi(R,phi,m=0,s=gamma,log=True)

  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p)

  return ll

                                                                             
##
## -------------------------------------------------------------------------- ##




## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                   Full likelihood for GEV marginals (Loc)
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians.
##
##

def loc0_gev_update_mixture_me_likelihood(data, params, Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ## Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_loc0 = params
  loc0 = data@beta_loc0  # mu = Xb
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
  Loc = Loc.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)<prob_below-0.05 or np.max(tmp)>prob_above+0.05:
      return -np.inf
  
  # cen = which_censored(Y, Loc, Scale, Shape, prob_below) # 'cen' isn't altered in Global
  # cen_above = ~which_censored(Y, Loc, Scale, Shape, prob_above)
  
  ## What if GEV params are such that all Y's are censored?
  if np.all(cen):
      return -np.inf
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above) 
  return ll

def loc0_interc_gev_update_mixture_me_likelihood(data, params, beta_loc0_1, Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ## Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_loc0_0 = params
  beta_loc0 = np.r_[beta_loc0_0,beta_loc0_1]
  loc0 = data@beta_loc0  # mu = Xb
  # loc0 = loc0.astype(float)
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
  Loc = Loc.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)< prob_below or np.max(tmp)>prob_above:
      return -np.inf
  
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above) 
  return ll

## For the slope wrt T of the location parameter
def loc1_gev_update_mixture_me_likelihood(data, params, Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc0, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ##Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_loc1 = params
  loc1 = data@beta_loc1  # mu = Xb
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
  Loc = Loc.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)<prob_below-0.05 or np.max(tmp)>prob_above+0.05:
      return -np.inf
  
  # cen = which_censored(Y, Loc, Scale, Shape, prob_below) # 'cen' isn't altered in Global
  # cen_above = ~which_censored(Y, Loc, Scale, Shape, prob_above)
  
  ## What if GEV params are such that all Y's are censored?
  if(np.all(cen)):
      return -np.inf
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above)
  return ll
                                                                             
##
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                  Full likelihood for GEV marginals (Scale)
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians.
##
##

def scale_gev_update_mixture_me_likelihood(data, params, Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ## Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_scale = params
  scale = data@beta_scale  # mu = Xb
  if np.any(scale < 0):
      return -np.inf
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Scale = np.tile(scale, n_t)
  Scale = Scale.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)<prob_below-0.05 or np.max(tmp)>prob_above+0.05:
      return -np.inf
  
  # cen = which_censored(Y, Loc, Scale, Shape, prob_below) # 'cen' isn't altered in Global
  # cen_above = ~which_censored(Y, Loc, Scale, Shape, prob_above)

  ## What if GEV params are such that all Y's are censored?
  if(np.all(cen)):
      return -np.inf
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above)
  return ll


def scale_interc_gev_update_mixture_me_likelihood(data, params, beta_scale_1, Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ## Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_scale_0 = params
  beta_scale = np.r_[beta_scale_0,beta_scale_1]
  scale = data@beta_scale  # mu = Xb
  if np.any(scale < 0):
      return -np.inf
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Scale = np.tile(scale, n_t)
  Scale = Scale.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)<prob_below or np.max(tmp)>prob_above:
      return -np.inf
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above)
  return ll

##
## -------------------------------------------------------------------------- ##







## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                  Full likelihood for GEV marginals (Scale)
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians.
##
##

def shape_gev_update_mixture_me_likelihood(data, params, Y, X_s, cen, cen_above, prob_below, prob_above,
                     tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ## Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_shape = params
  shape = data@beta_shape  # mu = Xb
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Shape = np.tile(shape, n_t)
  Shape = Shape.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)<prob_below-0.05 or np.max(tmp)>prob_above+0.05:
      return -np.inf
  
  # cen = which_censored(Y, Loc, Scale, Shape, prob_below) # 'cen' isn't altered in Global
  # cen_above = ~which_censored(Y, Loc, Scale, Shape, prob_above)

  ## What if GEV params are such that all Y's are censored?
  if(np.all(cen)):
      return -np.inf
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above)
  return ll


def shape_interc_gev_update_mixture_me_likelihood(data, params, beta_shape_1, Y, X_s, cen, cen_above, prob_below, prob_above,
                     tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above):
  
  ## Design_mat = data
  ## For the time being, assume that the intercept, slope are CONSTANTS
  beta_shape_0 = params
  beta_shape = np.r_[beta_shape_0,beta_shape_1]
  shape = data@beta_shape  # mu = Xb
  
  if len(X_s.shape)==1:
      X_s = X_s.reshape((X_s.shape[0],1))
  n_t = X_s.shape[1]
  n_s = X_s.shape[0]
  Shape = np.tile(shape, n_t)
  Shape = Shape.reshape((n_s,n_t),order='F')
  
  max_support = Loc - Scale/Shape
  max_support[Shape>0] = np.inf
  
  # When cen is not updated, the best thing we can do is to make sure the unifs is not too far from [below, above].
  tmp=pgev(Y[~cen & ~cen_above], Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if np.any(Y > max_support) or np.min(tmp)<prob_below or np.max(tmp)>prob_above:
      return -np.inf
  
  
  X = X_update(Y, cen, cen_above, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
  ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above)
  return ll

##
## -------------------------------------------------------------------------- ##






## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                 Full likelihood for latent Gaussian field
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians.
##
## Both functions in this section admit column vectors only.
##

def Z_likelihood_conditional_eigen(Z, V, d):
    # R_powered = R**phi
    part1 = -0.5*eig2inv_quadform_vector(V, 1/d, Z)-0.5*np.sum(np.log(d))
    return part1

def Z_likelihood_conditional(Z, Cor, cholesky_inv):
    # R_powered = R**phi
    part1 = -0.5*inv_quadform_vector(Cor, Z, cholesky_inv)-0.5*2*sum(np.log(np.diag(cholesky_inv[0])))
    return part1

def Z_update_onetime(Y, X, R, Z, cen, cen_above, prob_below, prob_above,
                tau_sqd, phi, gamma, Loc, Scale, Shape, xp, surv_p, den_p, 
                thresh_X, thresh_X_above, Cor, cholesky_inv, Sigma_m, random_generator):
    
    n_s = X.size
    prop_Z = np.empty(X.shape)
    accept = np.zeros(n_s)
    ## Generate X_s
    X_s = (R**phi)*norm_to_Pareto(Z)
    
    log_num=0; log_denom=0 # sd= np.sqrt(tau_sqd)
    for idx, Z_idx in enumerate(Z):
        # tripped : X = Y, changing X will change Y as well.
        prop_Z[:] = Z
        # temp = X_s(iter)+v_q(iter)*R::rnorm(0,1);
        temp = Z_idx + Sigma_m[idx]*random_generator.standard_normal(1)
        prop_Z[idx] = temp
        prop_X_s_idx = (R**phi)*norm_to_Pareto(temp)
        
        log_num = marg_transform_data_mixture_me_likelihood_uni(Y[idx], X[idx], prop_X_s_idx, 
                       cen[idx], cen_above[idx], prob_below, prob_above, Loc[idx], Scale[idx], Shape[idx], tau_sqd, phi, gamma, 
                       xp, surv_p, den_p, thresh_X, thresh_X_above) + Z_likelihood_conditional(prop_Z, Cor, cholesky_inv);
        log_denom = marg_transform_data_mixture_me_likelihood_uni(Y[idx], X[idx], X_s[idx], 
                       cen[idx], cen_above[idx], prob_below, prob_above, Loc[idx], Scale[idx], Shape[idx], tau_sqd, phi, gamma, 
                       xp, surv_p, den_p, thresh_X, thresh_X_above) + Z_likelihood_conditional(Z, Cor, cholesky_inv);
        
        with np.errstate(over='raise'):
            try:
                r = np.exp(log_num - log_denom)  # this gets caught and handled as an exception
            except FloatingPointError:
                print(" -- idx="+str(idx)+", Z="+str(Z[idx])+", prop_Z="+str(temp)+", log_num="+str(log_num)+", log_denom="+str(log_denom))
                r=0
    
        if random_generator.uniform(0,1,1)<r:
            Z[idx] = temp  # changes argument 'X_s' directly
            X_s[idx] = prop_X_s_idx
            accept[idx] = accept[idx] + 1
    
    #result = (X_s,accept)
    return accept

##
## -------------------------------------------------------------------------- ##







## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                   Full likelihood for Matern parameters
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## Update covariance parameters. 
##
##

def theta_c_update_mixture_me_likelihood_eigen(data, params, S, V=np.nan, d=np.nan):
  Z = data
  range = params[0]
  nu = params[1]
  if len(Z.shape)==1:
      Z = Z.reshape((Z.shape[0],1))
  n_t = Z.shape[1]
  
  if np.any(np.isnan(V)):
    Cor = corr_fn(S, np.array([range,nu]))
    eig_Cor = np.linalg.eigh(Cor) #For symmetric matrices
    V = eig_Cor[1]
    d = eig_Cor[0]

  ll = np.empty(n_t)
  ll[:]=np.nan
  for idx in np.arange(n_t):
    ll[idx] = Z_likelihood_conditional_eigen(Z[:,idx], V, d)
  return np.sum(ll)

def range_update_mixture_me_likelihood_eigen(data, params, nu, S, V=np.nan, d=np.nan):
  Z = data
  range = params
  
  if len(Z.shape)==1:
      Z = Z.reshape((Z.shape[0],1))
  n_t = Z.shape[1]
  
  if np.any(np.isnan(V)):
    Cor = corr_fn(S, np.array([range,nu]))
    eig_Cor = np.linalg.eigh(Cor) #For symmetric matrices
    V = eig_Cor[1]
    d = eig_Cor[0]

  ll = np.empty(n_t)
  ll[:]=np.nan
  for idx in np.arange(n_t):
    ll[idx] = Z_likelihood_conditional_eigen(Z[:,idx], V, d)
  return np.sum(ll)


def theta_c_update_mixture_me_likelihood(data, params, S, Cor=None, cholesky_inv=None):
  Z = data
  range = params[0]
  nu = params[1]
  if len(Z.shape)==1:
      Z = Z.reshape((Z.shape[0],1))
  n_t = Z.shape[1]
  
  if Cor is None:
    Cor = corr_fn(S, np.array([range,nu]))
    cholesky_inv = lapack.dposv(Cor,Z[:,0])

  ll = np.empty(n_t)
  ll[:]=np.nan
  for idx in np.arange(n_t):
    ll[idx] = Z_likelihood_conditional(Z[:,idx], Cor, cholesky_inv)
  return np.sum(ll)

def range_update_mixture_me_likelihood(data, params, nu, S, Cor=None, cholesky_inv=None):
  Z = data
  range = params
  
  if len(Z.shape)==1:
      Z = Z.reshape((Z.shape[0],1))
  n_t = Z.shape[1]
  
  if Cor is None:
    Cor = corr_fn(S, np.array([range,nu]))
    cholesky_inv = lapack.dposv(Cor,Z[:,0])

  ll = np.empty(n_t)
  ll[:]=np.nan
  for idx in np.arange(n_t):
    ll[idx] = Z_likelihood_conditional(Z[:,idx], Cor, cholesky_inv)
  return np.sum(ll)

##
## -------------------------------------------------------------------------- ##







## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                           Full likelihood for R
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians.
##
##

def Rt_update_mixture_me_likelihood(data, params, X, Z, cen, cen_above, 
                prob_below, prob_above, Loc, Scale, Shape, tau_sqd, phi, gamma,
                xp, surv_p, den_p, thresh_X, thresh_X_above):
  Y = data
  R = params
    
  if R < 0:
      return -np.inf
  else:
      ## Generate X_s
      X_s = (R**phi)*norm_to_Pareto(Z)
      ll = marg_transform_data_mixture_me_likelihood(Y, X, X_s, cen, cen_above, 
                prob_below, prob_above, Loc, Scale, Shape, tau_sqd, phi, gamma, 
                xp, surv_p, den_p, thresh_X, thresh_X_above)
      return ll
                                                                             
##
## -------------------------------------------------------------------------- ##




## New python code

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
##                           Gaussian moothing kernel
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## When d > fit radius, the weight will be zero
## h is tje bandwidth parameter
##

def weights_fun(d,radius,h=1, cutoff=True):
  if(isinstance(d, (int, np.int64, float))): d=np.array([d])
  tmp = np.exp(-d**2/(2*h))
  if cutoff: tmp[d>radius] = 0
  
  return tmp/np.sum(tmp)

                                                                             
##
## -------------------------------------------------------------------------- ##