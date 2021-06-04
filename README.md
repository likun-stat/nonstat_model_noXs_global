# Space scale-aware tail dependence modeling for high-dimensional spatial extremes

## Research ongoing with Dr. Mark Risser

[![Build Status](https://github.com/likun-stat/nonstat_model/actions/workflows/python-package-conda.yml/badge.svg?branch=main)](https://github.com/likun-stat/nonstat_model/actions)


We try to model dependence in extremes of spatial processes: 

![equation](http://latex.codecogs.com/gif.latex?%5C%7BX%28s%29%3Bs%5Cin%5Cmathcal%7BS%7D%5Csubset%5Cmathbb%7BR%7D%5E2%5C%7D)

### 1.  Flexible Spatial Modeling
We wish to fit spatial models encompassing both long-range asymptotic independence and short-range weakening dependence strength that leads to either asymptotic dependence or independence.

![equation](https://latex.codecogs.com/gif.latex?X%5E*%28%5Cboldsymbol%7Bs%7D%29%3DR%28%5Cboldsymbol%7Bs%7D%29%5E%7B%5Cphi%28%5Cboldsymbol%7Bs%7D%29%7DW%28%5Cboldsymbol%7Bs%7D%29%2C)

where

![equation](https://latex.codecogs.com/gif.latex?R%28%5Cboldsymbol%7Bs%7D%29%3D%5Csum_%7Bk%3D1%7D%5EK%20w_k%28%5Cboldsymbol%7Bs%7D%29%20S_k%20%5Ctext%7B%20with%20%7DS_k%5Csim%20%5Ctext%7BStable%7D%28%5Calpha%2C1%2C%5Cgamma_k%2C%5Cdelta%29.)

The kernel functions are centered at K mixture component locations, and they are proportional to Gaussian densities:

![equation](https://latex.codecogs.com/gif.latex?w_k%28%5Cboldsymbol%7Bs%7D%29%5Cpropto%20%5Cexp%5Cleft%5C%7B-%5Cfrac%7B%7C%7C%5Cboldsymbol%7Bs%7D-%5Cboldsymbol%7Bb%7D_k%7C%7C%5E2%7D%7B2%5Clambda_w%7D%5Cright%5C%7D).

### 2. Inference for the mixture componenet model

In contrast to Wadsworth and Tawn (2019) [[1]](#1), this approach is a coherent model with a well-defined likelihood:

![equation](https://latex.codecogs.com/gif.latex?%5Cvarphi%5Cleft%28Y_t%28%5Cboldsymbol%7Bs%7D_i%29%7CX%5E*_t%28%5Cboldsymbol%7Bs%7D%29%2C%5Cboldsymbol%7B%5Ctheta%7D%28%5Cboldsymbol%7Bs%7D%29%2Cp%2C%5Cphi%28%5Cboldsymbol%7Bs%7D%29%2C%5Ctau%5E2%5Cright%29%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20%5CPhi%5Cleft%28%5Cfrac%7BF_%7BX%7C%5Cphi_i%2C%5Ctau%5E2%7D%5E%7B-1%7D%28p%29-X%5E*_t%28%5Cboldsymbol%7Bs%7D_i%29%7D%7B%5Ctau%7D%5Cright%29%26%20%5Ctext%7Bif%20%7D%20Y_t%28%5Cboldsymbol%7Bs%7D_i%29%5Cleq%20u_t%28%5Cboldsymbol%7Bs%7D_i%29%2C%5C%5C%20%5Cphi%5Cleft%28F_X%5E%7B-1%7D%5Ccirc%20F_Y%28Y_t%28%5Cboldsymbol%7Bs%7D_i%29%29%5Crvert%20X%5E*_t%28%5Cboldsymbol%7Bs%7D_i%29%2C%5Ctau%5E2%5Cright%29%5Cfrac%7Bf_Y%28Y_t%28%5Cboldsymbol%7Bs%7D_i%29%29%7D%7Bf_X%5Cleft%28F_X%5E%7B-1%7D%5Ccirc%20F_Y%28Y_t%28%5Cboldsymbol%7Bs%7D_i%29%29%5Cright%29%7D%26%20%5Ctext%7Bif%20%7D%20Y_t%28%5Cboldsymbol%7Bs%7D_i%29%3E%20u_t%28%5Cboldsymbol%7Bs%7D_i%29.%20%5Cend%7Bmatrix%7D%5Cright.)

We sample the latent process as follows:

- Draw the smooth process _X*(s)_ given the noisy process _X(s)_.

(No truncation!)

- Draw the noise given the smooth process _X*(s)_.

(Truncated, but univariate!)

![gif](www/anime.gif)



## References
<a id="1">[1]</a> 
Wadsworth, J. L., & Tawn, J. (2019).
Higher-dimensional spatial extremes via single-site conditioning. 
arXiv preprint arXiv:1912.06560.
