import os
os.chdir("/Users/LikunZhang/Desktop/PyCode/")

# import random
import nonstat_model_noXs.model_sim as utils
import nonstat_model_noXs.generic_samplers as sampler
import nonstat_model_noXs.priors as priors
import numpy as np
from scipy.stats import uniform
from scipy.stats import norm 

# ------------ 1. Simulation settings -------------
range = 1        # Matern range
nu =  3/2        # Matern smoothness

n_s = 70         # Number of sites
n_t = 64         # Number of time points
tau_sqd = 10     # Nugget SD

phi = 0.55       # For R
gamma = 0.5
prob_below=0.8
prob_above=0.99

# -------------- 2. Generate covariance matrix -----------------
# Calculate distance between rows of 'Y', and return as distance matrix
np.random.seed(seed=1234)
from scipy.spatial import distance
Stations = np.c_[uniform.rvs(0,5,n_s),uniform.rvs(0,5,n_s)]
# plt.scatter(Stations[:,0],Stations[:,1])
S = distance.squareform(distance.pdist(Stations))
Cor = utils.corr_fn(S, np.array([range,nu]))
eig_Cor = np.linalg.eigh(Cor) #For symmetric matrices
V = eig_Cor[1]
d = eig_Cor[0]

R = utils.rlevy(n_t,s=gamma)
X = np.empty((n_s,n_t))
X[:] = np.nan
X_s = np.empty((n_s,n_t))
X_s[:] = np.nan
Z = np.empty((n_s,n_t))
Z[:] = np.nan
for idx, r in enumerate(R):
  Z_t = utils.eig2inv_times_vector(V, np.sqrt(d), norm.rvs(size=n_s))
  Z_to_W_s = 1/(1-norm.cdf(Z_t))
  tmp = (r**phi)*Z_to_W_s
  X_s[: ,idx] = tmp
  X[:,idx] = tmp + np.sqrt(tau_sqd)*norm.rvs(size=n_s)
  Z[:,idx] = Z_t


# ------------ 3. Marginal transformation -----------------
grid = utils.survival_interp_grid(phi, gamma, grid_size=800)
xp = grid[:,0]; surv_p = grid[:,1]
grid = utils.density_interp_grid(phi, gamma, grid_size=800)
den_p = grid[:,1]

Lon_lat = Stations
# Design_mat = np.c_[np.repeat(1,n_s), Lon_lat[:,1]]
Design_mat = np.c_[np.repeat(1,n_s), norm.rvs(size=n_s)]
n_covariates = Design_mat.shape[1]

beta_loc0 = np.array([0.2,-1])
loc0 = Design_mat @beta_loc0

beta_loc1 = np.array([0.1, -0.1])
loc1 = Design_mat @beta_loc1

Time = np.arange(n_t)
Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
Loc = Loc.reshape((n_s,n_t),order='F')


beta_scale = np.array([1.3,0.5])
scale = Design_mat @beta_scale
Scale = np.tile(scale, n_t)
Scale = Scale.reshape((n_s,n_t),order='F')

beta_shape = np.array([0.05,0.2])
shape = Design_mat @beta_shape
Shape = np.tile(shape, n_t)
Shape = Shape.reshape((n_s,n_t),order='F')

Y = utils.RW_me_2_gev(X, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
unifs = utils.pgev(Y, Loc, Scale, Shape)

cen = unifs < prob_below
cen_above = unifs > prob_above
thresh_X =  utils.qRW_me_interp(prob_below, xp, surv_p, tau_sqd, phi, gamma)
thresh_X_above =  utils.qRW_me_interp(prob_above, xp, surv_p, tau_sqd, phi, gamma)



# ------------ 4. Save initial values -----------------
initial_values   = {'phi':phi,
                    'gamma':gamma,
                    'tau_sqd':tau_sqd,
                    'prob_below':prob_below,
                    'prob_above':prob_above,
                    'Dist':S,
                    'theta_c':np.array([range,nu]),
                    'X':X,
                    'X_s':X_s,
                    'R':R,
                    'Z':Z,
                    'Design_mat':Design_mat,
                    'beta_loc0':beta_loc0,
                    'beta_loc1':beta_loc1,
                    'Time':Time,
                    'beta_scale':beta_scale,
                    'beta_shape':beta_shape,
                    }
n_updates = 1001    
sigma_m   = {'phi':2.4**2,
             'tau_sqd':2.4**2,
             'theta_c':2.4**2/2,
             'Z_onetime':np.repeat(np.sqrt(tau_sqd),n_s),
             'R_1t':2.4**2,
             'beta_loc0':2.4**2/n_covariates,
             'beta_loc1':2.4**2/n_covariates,
             'beta_scale':2.4**2/n_covariates,
             'beta_shape':2.4**2/n_covariates,
             }
prop_sigma   = {'theta_c':np.eye(2),
                'beta_loc0':np.eye(n_covariates),
                'beta_loc1':np.eye(n_covariates),
                'beta_scale':np.eye(n_covariates),
                'beta_shape':np.eye(n_covariates)
                }

from pickle import dump
with open('./test.pkl', 'wb') as f:
     dump(Y, f)
     dump(cen, f)
     dump(cen_above,f)
     dump(initial_values, f)
     dump(sigma_m, f)
     dump(prop_sigma, f)



## ---------------------------------------------------------
## ------------------------ For phi ------------------------
## ---------------------------------------------------------

import matplotlib.pyplot as plt

def test(phi):
    return utils.phi_update_mixture_me_likelihood(Y, phi, R, Z, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
                        tau_sqd, gamma)

Phi = np.arange(0.545,0.555,step=0.0001)
Lik = np.zeros(len(Phi))
for idx, phi_tmp in enumerate(Phi):
    Lik[idx] = test(phi_tmp)
plt.plot(Phi, Lik, color='black', linestyle='solid')
plt.axvline(phi, color='r', linestyle='--');


random_generator = np.random.RandomState()
Res = sampler.static_metr(Y, 0.58, utils.phi_update_mixture_me_likelihood, priors.interval_unif, 
                  np.array([0.1,0.7]),1000, 
                  random_generator,
                  np.nan, 0.005, True, 
                  R, Z, cen, cen_above, 
                  prob_below, prob_above, Loc, Scale, Shape, tau_sqd, gamma) 
plt.plot(np.arange(1000),Res['trace'][0,:],linestyle='solid')


Res = sampler.adaptive_metr(Y, 0.552, utils.phi_update_mixture_me_likelihood, priors.interval_unif, 
                          np.array([0.1,0.7]),5000,
                          random_generator,
                          np.nan, False, False,
                          .234, 10, .8, 10, 
                          R, Z, cen, cen_above, 
                          prob_below, prob_above, Loc, Scale, Shape, tau_sqd, gamma)
plt.plot(np.arange(5000)[3:],Res['trace'][0,3:],linestyle='solid')
plt.hlines(phi, 0, 5000, colors='r', linestyles='--');


## -------------------------------------------------------
## ----------------------- For tau -----------------------
## -------------------------------------------------------
# t_chosen = 24
# def test(tau_sqd):
#     return utils.tau_update_mixture_me_likelihood(Y[:,t_chosen], tau_sqd, X_s[:,t_chosen], cen[:,t_chosen], cen_above[:,t_chosen], 
#                     prob_below, prob_above, Loc[:,t_chosen], Scale[:,t_chosen], Shape[:,t_chosen], 
#                     phi, gamma, xp, surv_p, den_p)

# Tau = np.arange(8,12,step=0.1)
# Lik = np.zeros(len(Tau))
# for idx, t in enumerate(Tau):
#     Lik[idx] = test(t) 
# plt.plot(Tau, Lik, color='black', linestyle='solid')
# plt.axvline(tau_sqd, color='r', linestyle='--');

def test(tau_sqd):
    return utils.tau_update_mixture_me_likelihood(Y, tau_sqd, X_s, cen, cen_above, 
                    prob_below, prob_above, Loc, Scale, Shape, 
                    phi, gamma, xp, surv_p, den_p)

Tau = np.arange(8,12,step=0.1)
Lik = np.zeros(len(Tau))
for idx, t in enumerate(Tau):
    Lik[idx] = test(t) 
plt.plot(Tau, Lik, color='black', linestyle='solid')
plt.axvline(tau_sqd, color='r', linestyle='--');



Res = sampler.static_metr(Y, 2.0, utils.tau_update_mixture_me_likelihood, priors.invGamma_prior, 
                          np.array([0.1,0.1]),1000, 
                          random_generator,
                          np.nan, 2.1, True, 
                          X_s, cen, cen_above, 
                          prob_below, prob_above, Loc, Scale, Shape, 
                          phi, gamma, xp, surv_p, den_p)
plt.plot(np.arange(1000),Res['trace'][0,:],linestyle='solid')


Res = sampler.adaptive_metr(Y, 2, utils.tau_update_mixture_me_likelihood, priors.invGamma_prior, 
                          np.array([0.1,0.1]),5000,
                          random_generator,
                          np.nan, False, False,
                          .234, 10, .8, 10, 
                          X_s, cen, cen_above, 
                          prob_below, prob_above, Loc, Scale, Shape, 
                          phi, gamma, xp, surv_p, den_p)
plt.plot(np.arange(5000)[3:],Res['trace'][0,3:],linestyle='solid')
plt.hlines(tau_sqd, 0, 5000, colors='r', linestyles='--');





## --------------------------------------------------------------
## ----------------------- For GEV params -----------------------
## --------------------------------------------------------------
# (1) loc0: 0.2,-1
# ## Update cen and cen_above or not?
# def X_update_test(coef):
#     beta_loc0 = np.array([coef,-1])
#     loc0 = Design_mat@beta_loc0 
#     n_t = X_s.shape[1]
#     n_s = X_s.shape[0]
#     Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
#     Loc = Loc.reshape((n_s,n_t),order='F')
#     # cen = utils.which_censored(Y, Loc, Scale, Shape, prob_below) # 'cen' isn't altered in Global
#     # cen_above = ~utils.which_censored(Y, Loc, Scale, Shape, prob_above)
                              
#     X_tmp = utils.gev_2_RW_me(Y[~cen & ~cen_above], xp, surv_p, tau_sqd, phi, gamma, 
#                             Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above])
    
#     ll=utils.marg_transform_data_mixture_me_likelihood(Y[~cen & ~cen_above], X_tmp, X_s[~cen & ~cen_above], cen[~cen & ~cen_above], cen_above[~cen & ~cen_above], prob_below, prob_above, Loc[~cen & ~cen_above], Scale[~cen & ~cen_above], Shape[~cen & ~cen_above], 
#                         tau_sqd, phi, gamma, xp, surv_p, den_p, thresh_X, thresh_X_above) 
#     return ll


# Coef = np.arange(0.18,0.22,step=0.001)
# X_ups = np.zeros(len(Coef))
# for idx, coef in enumerate(Coef):
#     X_ups[idx] = X_update_test(coef)
# plt.plot(Coef, X_ups, linestyle='solid')    

def test(x):
    return utils.loc0_gev_update_mixture_me_likelihood(Design_mat, np.array([x,-1]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)


Coef = np.arange(0.18,0.22,step=0.001)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(0.2, color='r', linestyle='--');

def test(x):
    return utils.loc0_gev_update_mixture_me_likelihood(Design_mat, np.array([0.2,x]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)
Coef = np.arange(-1.2,-0.8,step=0.01)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(-1, color='r', linestyle='--');


Res = sampler.adaptive_metr(Design_mat, np.array([0.2,-1]), utils.loc0_gev_update_mixture_me_likelihood, 
                            priors.unif_prior, 20, 5000, random_generator, np.nan, True,
                            False, .234, 10, .8,  10,
                            Y, X_s, cen, cen_above, prob_below, prob_above, 
                            tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)

prop_Sigma=np.cov(Res['trace'])



plt.plot(np.arange(5000),Res['trace'][0,:], linestyle='solid')
plt.hlines(0.2, 0, 5000, colors='r', linestyles='--');

plt.plot(np.arange(5000),Res['trace'][1,:], linestyle='solid')
plt.hlines(-1, 0, 5000, colors='r', linestyles='--');
plt.plot(*Res['trace'])



def tmpf(x,y):
    return utils.loc0_gev_update_mixture_me_likelihood(Design_mat, np.array([x,y]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)
try_size = 50
x = np.linspace(0.18, 0.22, try_size)
y = np.linspace(-1.1, -0.92, try_size)

func = np.empty((try_size,try_size))
for idy,yi in enumerate(y):
    for idx,xi in enumerate(x):
         func[idy,idx] = tmpf(xi,yi)

plt.contourf(x, y, func, 20, cmap='RdGy')
plt.colorbar();

## Use prop_Sigma
Res = sampler.adaptive_metr_ratio(Design_mat, np.array([0.2,-1]), utils.loc0_gev_update_mixture_me_likelihood, 
                            priors.unif_prior, 20, 5000, random_generator, prop_Sigma, -0.2262189, 0.2827393557113686, True,
                            False, .234, 10, .8,  10,
                            Y, X_s, cen, cen_above, prob_below, prob_above, 
                            tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)


# (2) loc1: 0.1, -0.1
def test(x):
    return utils.loc1_gev_update_mixture_me_likelihood(Design_mat, np.array([x,-0.1]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc0, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)
Coef = np.arange(0.098,0.14,step=0.0005)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(0.1, color='r', linestyle='--');


def test(x):
    return utils.loc1_gev_update_mixture_me_likelihood(Design_mat, np.array([0.1,x]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, loc0, Scale, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)

Coef = np.arange(-0.11,-0.08,step=0.001)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(-0.1, color='r', linestyle='--');



Res = sampler.adaptive_metr(Design_mat, np.array([0.1,-0.1]), utils.loc1_gev_update_mixture_me_likelihood, 
                            priors.unif_prior, 20, 5000, random_generator, 
                            np.nan, True, False, .234, 10, .8,  10,
                            Y, X_s, cen, cen_above, prob_below, prob_above, 
                            tau_sqd, phi, gamma, loc0, Scale, Shape, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)



plt.plot(np.arange(5000),Res['trace'][0,:],linestyle='solid')
plt.hlines(0.1, 0, 5000, colors='r', linestyles='--');

plt.plot(np.arange(5000),Res['trace'][1,:],linestyle='solid')
plt.hlines(-0.1, 0, 5000, colors='r', linestyles='--');
plt.plot(*Res['trace'])



def tmpf(x,y):
    return utils.loc1_gev_update_mixture_me_likelihood(Design_mat, np.array([x,y]), 
                            Y, X_s, cen, cen_above, prob_below, prob_above, 
                            tau_sqd, phi, gamma, loc0, Scale, Shape, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)
try_size = 50
x = np.linspace(0.098, 0.103, try_size)
y = np.linspace(-0.103, -0.096, try_size)

func  = np.empty((try_size,try_size))
for idy,yi in enumerate(y):
    for idx,xi in enumerate(x):
         func[idy,idx] = tmpf(xi,yi)

plt.contourf(x, y, func, 20, cmap='RdGy')
plt.colorbar();



# (3) scale: 0.1,1
def test(x):
    return utils.scale_gev_update_mixture_me_likelihood(Design_mat, np.array([x,beta_scale[1]]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)
Coef = np.arange(beta_scale[0]-0.01,beta_scale[0]+0.01,step=0.0005)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(beta_scale[0], color='r', linestyle='--');


def test(x):
    return utils.scale_gev_update_mixture_me_likelihood(Design_mat, np.array([beta_scale[0],x]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                     tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)

Coef = np.arange(beta_scale[1]-0.01,beta_scale[1]+0.01,step=0.0005)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(beta_scale[1], color='r', linestyle='--');



Res = sampler.adaptive_metr(Design_mat, np.array([0.1,1]), utils.scale_gev_update_mixture_me_likelihood, 
                            priors.unif_prior, 20, 5000, random_generator, 
                            np.nan, True, False, .234, 10, .8,  10, 
                            Y, X_s, cen, cen_above, prob_below, prob_above, 
                            tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)


plt.plot(np.arange(5000),Res['trace'][0,:], linestyle='solid')
plt.hlines(0.1, 0, 5000, colors='r', linestyles='--');

plt.plot(np.arange(5000),Res['trace'][1,:], linestyle='solid')
plt.hlines(1, 0, 5000, colors='r', linestyles='--');
plt.plot(*Res['trace'])


def tmpf(x,y):
    return utils.scale_gev_update_mixture_me_likelihood(Design_mat, np.array([x,y]), Y, X_s, cen, cen_above, prob_below, prob_above, 
                            tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)
try_size = 50
x = np.linspace(0.097, 0.104, try_size)
y = np.linspace(0.996, 1.005, try_size)

func = np.empty((try_size,try_size))
for idy,yi in enumerate(y):
    for idx,xi in enumerate(x):
         func[idy,idx] = tmpf(xi,yi)

plt.contourf(x, y, func, 20, cmap='RdGy')
plt.colorbar();



# (4) shape: -0.02,0.2
def test(x):
    return utils.shape_gev_update_mixture_me_likelihood(Design_mat, np.array([x,0.2]), Y, X_s, cen, cen_above, prob_below, prob_above,
                     tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)
Coef = np.arange(-0.03,0.,step=0.0005)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(-0.02, color='r', linestyle='--');


def test(x):
    return utils.shape_gev_update_mixture_me_likelihood(Design_mat, np.array([-0.02,x]), Y, X_s, cen, cen_above, prob_below, prob_above,
                     tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                     thresh_X, thresh_X_above)

Coef = np.arange(0.18,0.22,step=0.001)
Lik = np.zeros(len(Coef))
for idx, coef in enumerate(Coef):
    Lik[idx] = test(coef)
plt.plot(Coef, Lik, color='black', linestyle='solid')
plt.axvline(0.2, color='r', linestyle='--');



Res = sampler.adaptive_metr(Design_mat, np.array([-0.02,0.2]), utils.shape_gev_update_mixture_me_likelihood, 
                            priors.unif_prior, 20, 5000, random_generator, 
                            np.nan, True, False, .234, 10, .8,  10,
                            Y, X_s, cen, cen_above, prob_below, prob_above,
                            tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)


plt.plot(np.arange(5000),Res['trace'][0,:], linestyle='solid')
plt.hlines(-0.02, 0, 5000, colors='r', linestyles='--');

plt.plot(np.arange(5000),Res['trace'][1,:], linestyle='solid')
plt.hlines(0.2, 0, 5000, colors='r', linestyles='--');
plt.plot(*Res['trace'])



def tmpf(x,y):
    return utils.shape_gev_update_mixture_me_likelihood(Design_mat, np.array([x,y]), Y, X_s, cen, cen_above, prob_below, prob_above,
                            tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                            thresh_X, thresh_X_above)
try_size = 50
x = np.linspace(-0.0215, -0.0185, try_size)
y = np.linspace(0.1985, 0.2020, try_size)

func = np.empty((try_size,try_size))
for idy,yi in enumerate(y):
    for idx,xi in enumerate(x):
         func[idy,idx] = tmpf(xi,yi)

plt.contourf(x, y, func, 20, cmap='RdGy')
plt.colorbar();









## -------------------------------------------------------
## --------------------- For theta_c ---------------------
## -------------------------------------------------------

# 1, 1.5
def test(x):
    return utils.theta_c_update_mixture_me_likelihood(Z, np.array([x,1.5]), S)

Range = np.arange(0.5,1.3,step=0.01)
Lik = np.zeros(len(Range))
for idx, r in enumerate(Range):
    Lik[idx] = test(r) 
plt.plot(Range, Lik, color='black', linestyle='solid')
plt.axvline(range, color='r', linestyle='--');


def test(x):
    return utils.theta_c_update_mixture_me_likelihood(Z, np.array([1,x]), S)

Nu = np.arange(0.9,1.8,step=0.01)
Lik = np.zeros(len(Nu))
for idx, r in enumerate(Nu):
    Lik[idx] = test(r) 
plt.plot(Nu, Lik, color='black', linestyle='solid')
plt.axvline(nu, color='r', linestyle='--');


Res = sampler.adaptive_metr(Z, np.array([1,1.5]), utils.theta_c_update_mixture_me_likelihood, 
                            priors.unif_prior, 20, 5000, random_generator,
                            np.nan, True, False, .234, 10, .8,  10,
                            S)

plt.plot(np.arange(5000),Res['trace'][0,:], linestyle='solid')
plt.hlines(1, 0, 5000, colors='r', linestyles='--');

plt.plot(np.arange(5000),Res['trace'][1,:], linestyle='solid')
plt.hlines(1.5, 0, 5000, colors='r', linestyles='--');
plt.plot(*Res['trace'])



def tmpf(x,y):
    return utils.theta_c_update_mixture_me_likelihood(Z, np.array([x,y]), S)

try_size = 50
x = np.linspace(0.85, 1.15, try_size)
y = np.linspace(1.37, 1.7, try_size)

func  = np.empty((try_size,try_size))
for idy,yi in enumerate(y):
    for idx,xi in enumerate(x):
         func[idy,idx] = tmpf(xi,yi)

plt.contourf(x, y, func, 20, cmap='RdGy')
plt.colorbar();




## -------------------------------------------------------
## ------------------------ For Rt -----------------------
## -------------------------------------------------------

# R[0] = 0.074
t_chosen = 0
def test(x):
    return utils.Rt_update_mixture_me_likelihood(Y[:,t_chosen], x, X[:,t_chosen], Z[:,t_chosen], cen[:,t_chosen], cen_above[:,t_chosen], 
                prob_below, prob_above, Loc[:,t_chosen], Scale[:,t_chosen], Shape[:,t_chosen], tau_sqd, phi, gamma,
                xp, surv_p, den_p, thresh_X, thresh_X_above) + priors.R_prior(x, gamma)

Rt = np.arange(0.01,0.1,step=0.001)
Lik = np.zeros(len(Rt))
for idx, r in enumerate(Rt):
    Lik[idx] = test(r) 
plt.plot(Rt, Lik, linestyle='solid')
plt.axvline(R[t_chosen], color='r', linestyle='--');



Res = sampler.adaptive_metr(Y[:,t_chosen], R[t_chosen], utils.Rt_update_mixture_me_likelihood, 
                            priors.R_prior, gamma, 5000, random_generator,
                            np.nan, True, False, .234, 10, .8,  10,
                            X[:,t_chosen], Z[:,t_chosen], cen[:,t_chosen], cen_above[:,t_chosen], 
                            prob_below, prob_above, Loc[:,t_chosen], Scale[:,t_chosen], Shape[:,t_chosen], tau_sqd, phi, gamma,
                            xp, surv_p, den_p, thresh_X, thresh_X_above)
plt.plot(np.arange(5000), Res['trace'][0,:], linestyle='solid')
plt.hlines(R[t_chosen], 0, 5000, colors='r', linestyles='--');


## -------------------------------------------------------
## ------------------------ For Z -----------------------
## -------------------------------------------------------

t_chosen = 1; idx = 20
prop_Z = np.empty(n_s)
prop_Z[:] = Z[:,t_chosen]
def test(x):
    prop_Z[idx] =x
    prop_X_s_idx = (R[t_chosen]**phi)*utils.norm_to_Pareto(x)
    return utils.marg_transform_data_mixture_me_likelihood_uni(Y[idx, t_chosen], X[idx, t_chosen], prop_X_s_idx, 
                       cen[idx, t_chosen], cen_above[idx, t_chosen], prob_below, prob_above, 
                       Loc[idx, t_chosen], Scale[idx, t_chosen], Shape[idx, t_chosen], tau_sqd, phi, gamma, 
                       xp, surv_p, den_p, thresh_X, thresh_X_above) + utils.Z_likelihood_conditional(prop_Z, V, d)
Zt = np.arange(Z[idx, t_chosen]-3,Z[idx, t_chosen]+1,step=0.01)
Lik = np.zeros(len(Zt))
for idy, x in enumerate(Zt):
    Lik[idy] = test(x) 
plt.plot(Zt, Lik, color='black', linestyle='solid')
plt.axvline(Z[idx, t_chosen], color='r', linestyle='--');


n_updates = 5000
Z_trace = np.empty((3,n_updates))

Z_new = np.empty(n_s)
Z_new[:] = Z[:, t_chosen]
K=10; k=3
r_opt = .234; c_0 = 10; c_1 = .8
accept = np.zeros(n_s)
Sigma_m = np.repeat(np.sqrt(tau_sqd),n_s)
for idx in np.arange(n_updates):
    tmp = utils.Z_update_onetime(Y[:,t_chosen], X[:,t_chosen], R[t_chosen],  Z_new, cen[:,t_chosen], cen_above[:,t_chosen], prob_below, prob_above,
                                   tau_sqd, phi, gamma, Loc[:,t_chosen], Scale[:,t_chosen], Shape[:,t_chosen], xp, surv_p, den_p,
                                   thresh_X, thresh_X_above, V, d, Sigma_m, random_generator)
    Z_trace[:,idx] = np.array([Z_new[4],Z_new[6],Z_new[7]])
    accept = accept + tmp
    
    if (idx % K) == 0:
        print('Finished ' + str(idx) + ' out of ' + str(n_updates) + ' iterations ')
        gamma2 = 1 / ((idx/K) + k)**(c_1)
        gamma1 = c_0*gamma2
        R_hat = accept/K
        Sigma_m = np.exp(np.log(Sigma_m) + gamma1*(R_hat - r_opt))
        accept[:] = 0


plt.plot(np.arange(n_updates),Z_trace[0,:],linestyle='solid')
plt.hlines(Z[4, t_chosen], 0, n_updates, colors='r', linestyles='--');

plt.plot(np.arange(n_updates),Z_trace[1,:],linestyle='solid')
plt.hlines(Z[6, t_chosen], 0, n_updates, colors='r', linestyles='--');

plt.plot(np.arange(n_updates),Z_trace[2,:],linestyle='solid')
plt.hlines(Z[7, t_chosen], 0, n_updates, colors='r', linestyles='--');



## -------------------------------------------------------
## --------------------- For Sampler ---------------------
## -------------------------------------------------------
from pickle import load
with open('nonstat_progress_0.pkl', 'rb') as f:
     Y_tmp=load(f)
     cen_tmp=load(f)
     cen_above_tmp=load(f)
     initial_values_tmp=load(f)
     sigma_m=load(f)
     prop_sigma=load(f)
     iter_tmp=load(f)
     phi_trace_tmp=load(f)
     tau_sqd_trace=load(f)
     theta_c_trace_tmp=load(f)
     beta_loc0_trace_tmp=load(f)
     beta_loc1_trace_tmp=load(f)
     beta_scale_trace_tmp=load(f)
     beta_shape_trace_tmp=load(f)
                   
     X_s_1t_trace=load(f)
     R_1t_trace=load(f)
     Y_onetime=load(f)
     X_onetime=load(f)
     X_s_onetime=load(f)
     R_onetime=load(f)





loc0=initial_values_tmp['Design_mat']@initial_values_tmp['beta_loc0']
loc1=initial_values_tmp['Design_mat']@initial_values_tmp['beta_loc1']
scale=initial_values_tmp['Design_mat']@initial_values_tmp['beta_scale']
shape=initial_values_tmp['Design_mat']@initial_values_tmp['beta_shape']

n_t = Y_tmp.shape[1]
n_s = Y_tmp.shape[0]
Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(initial_values_tmp['Time'],n_s)
Loc = Loc.reshape((n_s,n_t),order='F')
Scale = np.tile(scale, n_t)
Scale = Scale.reshape((n_s,n_t),order='F')
Shape = np.tile(shape, n_t)
Shape = Shape.reshape((n_s,n_t),order='F')