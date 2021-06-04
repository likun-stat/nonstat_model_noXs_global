


###################################################################################
## Main sampler

## Must provide data input 'data_input.pkl' to initiate the sampler.
## In 'data_input.pkl', one must include 
##      Y ........................................... censored observations on GEV scale
##      cen ........................................................... indicator matrix
##      initial.values ........ a dictionary: phi, tau_sqd, prob_below, prob_above, Dist, 
##                                             theta_c, X, X_s, R, Design_mat, beta_loc0, 
##                                             beta_loc1, Time, beta_scale, beta_shape
##      n_updates .................................................... number of updates
##      thinning ......................................... number of runs in each update
##      experiment_name
##      echo_interval ......................... echo process every echo_interval updates
##      sigma_m
##      prop_Sigma
##      true_params ....................... a dictionary: phi, rho, tau_sqd, theta_gpd, 
##                                              prob_below, X_s, R
##

 
      

if __name__ == "__main__":
    import os
    from pickle import load
    from pickle import dump
    import nonstat_model_noXs.model_sim as utils
    import nonstat_model_noXs.generic_samplers as sampler
    import nonstat_model_noXs.priors as priors
    import numpy as np
    from scipy.spatial import distance
    from scipy.stats import norm 
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.linalg import lapack
    
    # Check whether the 'mpi4py' is installed
    test_mpi = os.system("python -c 'from mpi4py import *' &> /dev/null")
    if test_mpi != 0:
        import sys
        sys.exit("mpi4py import is failing, aborting...")
       
    # get rank and size
    from mpi4py import MPI
      
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    ## ---------------------------------------------------------------------------
    ##            Automatically make a directory for the new simulation
    ## ---------------------------------------------------------------------------
    run = 1
    save_directory = "Simulation_"+str(run)
    if rank==0:
        dirs = os.listdir()
        while save_directory in dirs:
            run+=1
            save_directory = "Simulation_"+str(run)
        os.mkdir(save_directory)
    run = comm.bcast(run,root=0)
    save_directory = "Simulation_"+str(run)   
    
    ## -------------------------------------------------------
    ##                     General setup
    ## -------------------------------------------------------
    # size=64;rank=0
    thinning = 10; echo_interval = 20; n_updates = 50001
      
    # Generate multiple independent random streams
    random_generator = np.random.RandomState()
    
    # Constants to control adaptation of the Metropolis sampler
    c_0 = 10
    c_1 = 0.8
    offset = 3  # the iteration offset
    r_opt_1d = .41
    r_opt_2d = .35
    eps = 1e-6 # a small number
      
    # Hyper parameters for the prior of the mixing distribution parameters and 
    hyper_params_phi = np.array([0.1,0.7])
    hyper_params_tau_sqd = np.array([0.1,0.1])
    hyper_params_theta_c = np.array([0, 20])
    hyper_params_theta_gev = 25
    # hyper_params_range = np.array([0.5,1.5]) # in case where roughness is not updated
      
    # Load simulated data
    data_filename = 'data_sim'+str(run)+'.pkl'
    with open(data_filename, 'rb') as f:
        Y_all = load(f)
        cen_all = load(f)
        cen_above_all = load(f)
        data_all = load(f)
        sigma_m = load(f)
        prop_sigma = load(f)
        f.close()
    
    ## -------------------------------------------------------
    ##            Get the data for the local fit
    ## -------------------------------------------------------
    local_fit_no = 0
    radius = data_all['radius']
    
    # Filename for storing the intermediate results
    filename = './'+save_directory + '/local_' + str(local_fit_no) + '_progress_' + str(rank) + '.pkl'
    
    # Subset the stations within radius of the knot_of_interest
    knot_of_interest = data_all['Knots'][local_fit_no,:]
    Dist_from_knot = distance.cdist(knot_of_interest.reshape((-1,2)),data_all['Stations'])
    subset_indices = (Dist_from_knot[0,:] <= radius)
    Stations_local = data_all['Stations'][subset_indices,:]
    
    # import matplotlib.pyplot as plt
    # circle = plt.Circle((knot_of_interest[0],knot_of_interest[1]), radius, color='r', fill=False)
    # ax = plt.gca()
    # ax.cla() # clear things for fresh plot
    # ax.set_xlim((0, 10))
    # ax.set_ylim((0, 10))
    # ax.scatter(data_all['Stations'][:,0], data_all['Stations'][:,1], c='gray')
    # ax.scatter(Stations_local[:,0], Stations_local[:,1], c='r')
    # ax.add_patch(circle)
    
    # Subset the observations
    Y = Y_all[subset_indices,:]
    # cen = cen_all[subset_indices,:]
    # cen_above = cen_above_all[subset_indices,:]
    
    # Bookkeeping
    n_s = Y.shape[0]
    n_t = Y.shape[1]
    if n_t != size:
       import sys
       sys.exit("Make sure the number of cpus (N) = number of time replicates (n_t), i.e.\n     srun -N python nonstat_sampler.py")
    
    n_updates_thinned = np.int(np.ceil(n_updates/thinning))
    wh_to_plot_Xs = n_s*np.array([0.25,0.5,0.75])
    wh_to_plot_Xs = wh_to_plot_Xs.astype(int)
    sigma_m['Z_onetime'] = sigma_m['Z_onetime'][:n_s]
    
    
    # plt.scatter(data_all['Stations'][:,0], data_all['Stations'][:,1],c=data_all['phi_vec'], marker='o', alpha=0.5, cmap='jet')
    # plt.colorbar()
    
    # plt.scatter(data_all['Stations'][:,0], data_all['Stations'][:,1],c=data_all['range_vec'], marker='o', alpha=0.5, cmap='jet')
    # plt.colorbar()
    
    ## -------------------------------------------------------
    ##                  Set initial values
    ## -------------------------------------------------------
    phi = data_all['phi_at_knots'][local_fit_no]
    gamma = data_all['gamma']
    tau_sqd = data_all['tau_sqd']
    prob_below = data_all['prob_below']
    prob_above = data_all['prob_above']
    range = data_all['range_at_knots'][local_fit_no]
    nu = data_all['nu']
    
    # 1. For current values of phi and gamma, obtain grids of survival probs and densities
    grid = utils.density_interp_grid(phi, gamma, grid_size=800)
    xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
    thresh_X =  utils.qRW_me_interp(prob_below, xp, surv_p, tau_sqd, phi, gamma)
    thresh_X_above =  utils.qRW_me_interp(prob_above, xp, surv_p, tau_sqd, phi, gamma)
    
    # 2. Marginal GEV parameters: per location x time
    Design_mat = data_all['Design_mat'][subset_indices,:]
    beta_loc0 = data_all['beta_loc0']
    beta_loc1 = data_all['beta_loc1']
    Time = data_all['Time']
    beta_scale = data_all['beta_scale']
    beta_shape = data_all['beta_shape']
    
    loc0 = Design_mat @beta_loc0
    loc1 = Design_mat @beta_loc1
    Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
    Loc = Loc.reshape((n_s,n_t),order='F')
    
    scale = Design_mat @beta_scale
    Scale = np.tile(scale, n_t)
    Scale = Scale.reshape((n_s,n_t),order='F')
       
    # Design_mat1 = np.c_[np.repeat(1,n_s), np.log(Design_mat[:,1])]
    shape = Design_mat @beta_shape
    Shape = np.tile(shape, n_t)
    Shape = Shape.reshape((n_s,n_t),order='F')
    
    unifs = utils.pgev(Y, Loc, Scale, Shape)
    
    cen = unifs < prob_below
    cen_above = unifs > prob_above
    
    
    # 3. Eigendecomposition of the correlation matrix
    n_covariates = len(beta_loc0)
    theta_c = np.array([range,nu])
    Dist = distance.squareform(distance.pdist(Stations_local))
    tmp_vec = np.ones(n_s)
    Cor = utils.corr_fn(Dist, theta_c)
    # eig_Cor = np.linalg.eigh(Cor) #For symmetric matrices
    # V = eig_Cor[1]
    # d = eig_Cor[0]
    cholesky_inv = lapack.dposv(Cor,tmp_vec)
    
    # 4. Process data given initial values
    # X = data_all['X'][subset_indices,:]
    # X_s = data_all['X_s'][subset_indices,:]
    # Z = data_all['Z'][subset_indices,:]
    X = utils.gev_2_RW_me(Y, xp, surv_p, tau_sqd, phi, gamma, Loc, Scale, Shape)
    R = data_all['R_at_knots'][local_fit_no,:]
    
    Z = np.empty((n_s,n_t))
    Z[:] = np.nan
    for idx in np.arange(n_t):
        X_s_tmp = X[:,idx]-np.sqrt(tau_sqd)*norm.rvs(size=n_s)
        lower_limit = R[idx]**phi
        X_s_tmp[X_s_tmp<lower_limit] = lower_limit + 0.01
        Z[:,idx] = norm.ppf(1-1/(X_s_tmp/(R[idx]**phi)))
    
    # import matplotlib.pyplot as plt
    # plt.scatter(Stations_local[:,0], Stations_local[:,1],c=Z[:,rank], marker='o', alpha=0.5, cmap='jet')
    # plt.colorbar()
    # plt.title("Z onetime");
    
    v_q=np.repeat(2.4**2,n_s)
    for idx in np.arange(n_t):
        tmp = utils.Z_update_onetime(Y[:,idx], X[:,idx], R[idx], Z[:,idx], cen[:,idx], cen_above[:,idx], prob_below, prob_above,
                                       tau_sqd, phi, gamma, Loc[:,idx], Scale[:,idx], Shape[:,idx], xp, surv_p, den_p,
                                       thresh_X, thresh_X_above, Cor, cholesky_inv, v_q, random_generator)
    
    
    Y_onetime = Y[:,rank]
    X_onetime = X[:,rank]
    R_onetime = R[rank]
    X_s_onetime = (R_onetime**phi)*utils.norm_to_Pareto(Z[:,rank])
    Z_onetime = Z[:,rank]
    
    
    
    # Initial trace objects
    Z_1t_accept = np.zeros(n_s)
    R_accept = 0
    Z_1t_trace = np.empty((n_s,n_updates_thinned)); Z_1t_trace[:] = np.nan
    Z_1t_trace[:,0] = Z_onetime  
    R_1t_trace = np.empty(n_updates_thinned); R_1t_trace[:] = np.nan
    R_1t_trace[0] = R_onetime
    if rank == 0:
       print("Number of time replicates = %d"%size)
       X_s = np.empty((n_s,n_t))
       phi_trace = np.empty(n_updates_thinned); phi_trace[:] = np.nan
       phi_trace[0] = phi
       tau_sqd_trace = np.empty(n_updates_thinned); tau_sqd_trace[:] = np.nan
       tau_sqd_trace[0] = tau_sqd
       theta_c_trace_within_thinning = np.empty((2,thinning)); theta_c_trace_within_thinning[:] = np.nan
       theta_c_trace = np.empty((2,n_updates_thinned)); theta_c_trace[:] = np.nan
       theta_c_trace[:,0] = theta_c
       beta_loc0_trace_within_thinning = np.empty((n_covariates,thinning)); beta_loc0_trace_within_thinning[:] = np.nan
       beta_loc0_trace = np.empty((n_covariates,n_updates_thinned)); beta_loc0_trace[:] = np.nan
       beta_loc0_trace[:,0] = beta_loc0
       beta_loc1_trace_within_thinning = np.empty((n_covariates,thinning)); beta_loc1_trace_within_thinning[:] = np.nan
       beta_loc1_trace = np.empty((n_covariates,n_updates_thinned)); beta_loc1_trace[:] = np.nan
       beta_loc1_trace[:,0] = beta_loc1
       beta_scale_trace_within_thinning = np.empty((n_covariates,thinning)); beta_scale_trace_within_thinning[:] = np.nan
       beta_scale_trace = np.empty((n_covariates,n_updates_thinned)); beta_scale_trace[:] = np.nan
       beta_scale_trace[:,0] = beta_scale
       beta_shape_trace_within_thinning = np.empty((n_covariates,thinning)); beta_shape_trace_within_thinning[:] = np.nan
       beta_shape_trace = np.empty((n_covariates,n_updates_thinned)); beta_shape_trace[:] = np.nan
       beta_shape_trace[:,0] = beta_shape
        
       phi_accept = 0
       tau_sqd_accept = 0
       theta_c_accept = 0
       beta_loc0_accept = 0
       beta_loc1_accept = 0
       beta_scale_accept = 0
       beta_shape_accept = 0
    
    
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    # --------------------------- Start Metropolis Updates ------------------------------
    # -----------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------
    for iter in np.arange(1,n_updates):
        # Update X
        # print(str(rank)+" "+str(iter)+" Gathered? "+str(np.where(~cen)))
        X_onetime = utils.X_update(Y_onetime, cen[:,rank], cen_above[:,rank], xp, surv_p, tau_sqd, phi, gamma, Loc[:,rank], Scale[:,rank], Shape[:,rank])
          
        # Update Z
        tmp = utils.Z_update_onetime(Y_onetime, X_onetime, R_onetime, Z_onetime, cen[:,rank], cen_above[:,rank], prob_below, prob_above,
                                    tau_sqd, phi, gamma, Loc[:,rank], Scale[:,rank], Shape[:,rank], xp, surv_p, den_p,
                                    thresh_X, thresh_X_above, Cor, cholesky_inv, sigma_m['Z_onetime'], random_generator)
        Z_1t_accept = Z_1t_accept + tmp
          
        # Update R
        Metr_R = sampler.static_metr(Y_onetime, R_onetime, utils.Rt_update_mixture_me_likelihood, 
                        priors.R_prior, gamma, 2, 
                        random_generator,
                        np.nan, sigma_m['R_1t'], False, 
                        X_onetime, Z_onetime, cen[:,rank], cen_above[:,rank], 
                        prob_below, prob_above, Loc[:,rank], Scale[:,rank], Shape[:,rank], tau_sqd, phi, gamma,
                        xp, surv_p, den_p, thresh_X, thresh_X_above)
        R_accept = R_accept + Metr_R['acc_prob']
        R_onetime = Metr_R['trace'][0,1]
           
        X_s_onetime = (R_onetime**phi)*utils.norm_to_Pareto(Z_onetime)
    
        # *** Gather items ***
        X_s_recv = comm.gather(X_s_onetime,root=0)
        X_recv = comm.gather(X_onetime, root=0)
        Z_recv = comm.gather(Z_onetime, root=0)
        R_recv = comm.gather(R_onetime, root=0)
    
        if rank==0:
            X_s[:] = np.vstack(X_s_recv).T
            X[:] = np.vstack(X_recv).T
            # Check whether X is negative
            if np.any(X[~cen & ~cen_above]<0):
                sys.exit("X value abnormalty "+str(phi)+" "+str(tau_sqd))
                   
            Z[:] = np.vstack(Z_recv).T
            R[:] = R_recv
            index_within = (iter-1)%thinning
            # print('beta_shape_accept=',beta_shape_accept, ', iter=', iter)
    
            # Update phi
            Metr_phi = sampler.static_metr(Y, phi, utils.phi_update_mixture_me_likelihood, priors.interval_unif, 
                    hyper_params_phi, 2, 
                    random_generator,
                    np.nan, sigma_m['phi'], False, 
                    R, Z, cen, cen_above, 
                    prob_below, prob_above, Loc, Scale, Shape, tau_sqd, gamma)
            phi_accept = phi_accept + Metr_phi['acc_prob']
            phi = Metr_phi['trace'][0,1]
               
            # Update gamma (TBD)
            #
               
            grid = utils.density_interp_grid(phi, gamma, grid_size=800)
            xp = grid[0]; den_p = grid[1]; surv_p = grid[2]
            X_s = (R**phi)*utils.norm_to_Pareto(Z)
               
            # Update tau_sqd
            Metr_tau_sqd = sampler.static_metr(Y, tau_sqd, utils.tau_update_mixture_me_likelihood, priors.invGamma_prior, 
                        hyper_params_tau_sqd, 2, 
                        random_generator,
                        np.nan, sigma_m['tau_sqd'], False,
                        X_s, cen, cen_above, 
                        prob_below, prob_above, Loc, Scale, Shape, 
                        phi, gamma, xp, surv_p, den_p)
            tau_sqd_accept = tau_sqd_accept + Metr_tau_sqd['acc_prob']
            tau_sqd = Metr_tau_sqd['trace'][0,1]
              
            thresh_X =  utils.qRW_me_interp(prob_below, xp, surv_p, tau_sqd, phi, gamma)
            thresh_X_above =  utils.qRW_me_interp(prob_above, xp, surv_p, tau_sqd, phi, gamma)
    
               
            # Update theta_c
            # Metr_theta_c = sampler.static_metr(Z, theta_c, utils.theta_c_update_mixture_me_likelihood, 
            #                   priors.interval_unif_multi, hyper_params_theta_c, 2,
            #                   random_generator,
            #                   prop_sigma['theta_c'], sigma_m['theta_c'], False,
            #                   Dist)
            # theta_c_accept = theta_c_accept + Metr_theta_c['acc_prob']
            # theta_c = Metr_theta_c['trace'][:,1]
            # theta_c_trace_within_thinning[:,index_within] = theta_c
            Metr_theta_c = sampler.static_metr(Z, theta_c[0], utils.range_update_mixture_me_likelihood, 
                              priors.interval_unif, hyper_params_theta_c, 2,
                              random_generator,
                              np.nan, sigma_m['range'], False,
                              theta_c[1],Dist)
            theta_c_accept = theta_c_accept + Metr_theta_c['acc_prob']
            theta_c = np.array([Metr_theta_c['trace'][0,1],theta_c[1]])
            theta_c_trace_within_thinning[:,index_within] = theta_c
              
            if Metr_theta_c['acc_prob']>0:
                Cor = utils.corr_fn(Dist, theta_c)
                # eig_Cor = np.linalg.eigh(Cor) #For symmetric matrices
                # V = eig_Cor[1]
                # d = eig_Cor[0]
                cholesky_inv = lapack.dposv(Cor,tmp_vec)
               
            # Update beta_loc0
            # Metr_beta_loc0 = sampler.static_metr(Design_mat, beta_loc0, utils.loc0_gev_update_mixture_me_likelihood, 
            #                    priors.unif_prior, hyper_params_theta_gev, 2,
            #                    random_generator,
            #                    prop_sigma['beta_loc0'], sigma_m['beta_loc0'], False, 
            #                    Y, X_s, cen, cen_above, prob_below, prob_above, 
            #                    tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
            #                    thresh_X, thresh_X_above)
            # beta_loc0_accept = beta_loc0_accept + Metr_beta_loc0['acc_prob']
            # beta_loc0 = Metr_beta_loc0['trace'][:,1]
            # beta_loc0_trace_within_thinning[:,index_within] = beta_loc0
            # loc0 = Design_mat @beta_loc0
            Metr_beta_loc0 = sampler.static_metr(Design_mat, beta_loc0[0], utils.loc0_interc_gev_update_mixture_me_likelihood, 
                               priors.unif_prior_1dim, hyper_params_theta_gev, 2,
                               random_generator,
                               np.nan, sigma_m['beta_loc0'], False, 
                               beta_loc0[1], Y, X_s, cen, cen_above, prob_below, prob_above, 
                               tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
                               thresh_X, thresh_X_above)
            beta_loc0_accept = beta_loc0_accept + Metr_beta_loc0['acc_prob']
            beta_loc0 = np.array([Metr_beta_loc0['trace'][0,1],beta_loc0[1]])
            beta_loc0_trace_within_thinning[:,index_within] = beta_loc0
            loc0 = Design_mat @beta_loc0
              
            # Update beta_loc1
            # Metr_beta_loc1 = sampler.static_metr(Design_mat, beta_loc1, utils.loc1_gev_update_mixture_me_likelihood, 
            #                    priors.unif_prior, hyper_params_theta_gev, 2,
            #                    random_generator,
            #                    prop_sigma['beta_loc1'], sigma_m['beta_loc1'], False, 
            #                    Y, X_s, cen, cen_above, prob_below, prob_above, 
            #                    tau_sqd, phi, gamma, loc0, Scale, Shape, Time, xp, surv_p, den_p, 
            #                    thresh_X, thresh_X_above)
            # beta_loc1_accept = beta_loc1_accept + Metr_beta_loc1['acc_prob']
            # beta_loc1 = Metr_beta_loc1['trace'][:,1]
            # beta_loc1_trace_within_thinning[:,index_within] = beta_loc1
            # loc1 = Design_mat @beta_loc1
            Loc = np.tile(loc0, n_t) + np.tile(loc1, n_t)*np.repeat(Time,n_s)
            Loc = Loc.reshape((n_s,n_t),order='F')
               
            # Update beta_scale
            # Metr_beta_scale = sampler.static_metr(Design_mat, beta_scale, utils.scale_gev_update_mixture_me_likelihood, 
            #                     priors.unif_prior, hyper_params_theta_gev, 2,
            #                     random_generator,
            #                     prop_sigma['beta_scale'], sigma_m['beta_scale'], False,
            #                     Y, X_s, cen, cen_above, prob_below, prob_above, 
            #                     tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
            #                     thresh_X, thresh_X_above)
            # beta_scale_accept = beta_scale_accept + Metr_beta_scale['acc_prob']
            # beta_scale = Metr_beta_scale['trace'][:,1]
            # beta_scale_trace_within_thinning[:,index_within] = beta_scale
            # scale = Design_mat @beta_scale
            # Scale = np.tile(scale, n_t)
            # Scale = Scale.reshape((n_s,n_t),order='F')
            Metr_beta_scale = sampler.static_metr(Design_mat, beta_scale[0], utils.scale_interc_gev_update_mixture_me_likelihood, 
                                priors.unif_prior_1dim, hyper_params_theta_gev, 2,
                                random_generator,
                                np.nan, sigma_m['beta_scale'], False,
                                beta_scale[1], Y, X_s, cen, cen_above, prob_below, prob_above, 
                                tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
                                thresh_X, thresh_X_above)
            beta_scale_accept = beta_scale_accept + Metr_beta_scale['acc_prob']
            beta_scale = np.array([Metr_beta_scale['trace'][0,1],beta_scale[1]])
            beta_scale_trace_within_thinning[:,index_within] = beta_scale
            scale = Design_mat @beta_scale
            Scale = np.tile(scale, n_t)
            Scale = Scale.reshape((n_s,n_t),order='F')
              
            # Update beta_shape
            # Metr_beta_shape = sampler.static_metr(Design_mat, beta_shape, utils.shape_gev_update_mixture_me_likelihood, 
            #                     priors.unif_prior, hyper_params_theta_gev, 2, 
            #                     random_generator,
            #                     prop_sigma['beta_shape'], sigma_m['beta_shape'], False,
            #                     Y, X_s, cen, cen_above, prob_below, prob_above,
            #                     tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
            #                     thresh_X, thresh_X_above)
            # beta_shape_accept = beta_shape_accept + Metr_beta_shape['acc_prob']
            # beta_shape = Metr_beta_shape['trace'][:,1]
            # beta_shape_trace_within_thinning[:,index_within] = beta_shape
            # # shape = Design_mat1 @beta_shape
            # shape = Design_mat @beta_shape
            # Shape = np.tile(shape, n_t)
            # Shape = Shape.reshape((n_s,n_t),order='F')
            Metr_beta_shape = sampler.static_metr(Design_mat, beta_shape[0], utils.shape_interc_gev_update_mixture_me_likelihood, 
                                priors.unif_prior_1dim, hyper_params_theta_gev, 2, 
                                random_generator,
                                np.nan, sigma_m['beta_shape'], False,
                                beta_shape[1], Y, X_s, cen, cen_above, prob_below, prob_above,
                                tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
                                thresh_X, thresh_X_above)
            beta_shape_accept = beta_shape_accept + Metr_beta_shape['acc_prob']
            beta_shape = np.array([Metr_beta_shape['trace'][0,1],beta_shape[1]])
            beta_shape_trace_within_thinning[:,index_within] = beta_shape
            # shape = Design_mat1 @beta_shape
            shape = Design_mat @beta_shape
            Shape = np.tile(shape, n_t)
            Shape = Shape.reshape((n_s,n_t),order='F')
              
            # cen[:] = utils.which_censored(Y, Loc, Scale, Shape, prob_below)
            # cen_above[:] = ~utils.which_censored(Y, Loc, Scale, Shape, prob_above)
               
               
        # *** Broadcast items ***
        phi = comm.bcast(phi,root=0)
        xp = comm.bcast(xp,root=0)
        den_p = comm.bcast(den_p,root=0)
        surv_p = comm.bcast(surv_p,root=0)
        tau_sqd = comm.bcast(tau_sqd,root=0)
        thresh_X = comm.bcast(thresh_X,root=0)
        thresh_X_above = comm.bcast(thresh_X_above,root=0)
        theta_c = comm.bcast(theta_c,root=0)
        # V = comm.bcast(V,root=0)
        # d = comm.bcast(d,root=0)
        Cor = comm.bcast(Cor,root=0)
        cholesky_inv = comm.bcast(cholesky_inv,root=0)
        Loc = comm.bcast(Loc,root=0)
        Scale = comm.bcast(Scale,root=0)
        Shape = comm.bcast(Shape,root=0)
        # cen = comm.bcast(cen,root=0)
        # cen_above = comm.bcast(cen_above,root=0)
          
        # ----------------------------------------------------------------------------------------
        # --------------------------- Summarize every 'thinning' steps ---------------------------
        # ----------------------------------------------------------------------------------------
        if (iter % thinning) == 0:
            index = np.int(iter/thinning)
               
            # Fill in trace objects
            Z_1t_trace[:,index] = Z_onetime  
            R_1t_trace[index] = R_onetime
            if rank == 0:
                phi_trace[index] = phi
                tau_sqd_trace[index] = tau_sqd
                theta_c_trace[:,index] = theta_c
                beta_loc0_trace[:,index] = beta_loc0
                beta_loc1_trace[:,index] = beta_loc1
                beta_scale_trace[:,index] = beta_scale
                beta_shape_trace[:,index] = beta_shape
              
                
            # Adapt via Shaby and Wells (2010)
            gamma2 = 1 / (index + offset)**(c_1)
            gamma1 = c_0*gamma2
            sigma_m['Z_onetime'] = np.exp(np.log(sigma_m['Z_onetime']) + gamma1*(Z_1t_accept/thinning - r_opt_1d))
            Z_1t_accept[:] = 0
            sigma_m['R_1t'] = np.exp(np.log(sigma_m['R_1t']) + gamma1*(R_accept/thinning - r_opt_1d))
            R_accept = 0
              
            if rank == 0:
                sigma_m['phi'] = np.exp(np.log(sigma_m['phi']) + gamma1*(phi_accept/thinning - r_opt_1d))
                phi_accept = 0
                sigma_m['tau_sqd'] = np.exp(np.log(sigma_m['tau_sqd']) + gamma1*(tau_sqd_accept/thinning - r_opt_1d))
                tau_sqd_accept = 0
              
                sigma_m['range'] = np.exp(np.log(sigma_m['range']) + gamma1*(theta_c_accept/thinning - r_opt_1d))
                theta_c_accept = 0
                # sigma_m['theta_c'] = np.exp(np.log(sigma_m['theta_c']) + gamma1*(theta_c_accept/thinning - r_opt_2d))
                # theta_c_accept = 0
                # prop_sigma['theta_c'] = prop_sigma['theta_c'] + gamma2*(np.cov(theta_c_trace_within_thinning) - prop_sigma['theta_c'])
                # check_chol_cont = True
                # while check_chol_cont:
                #     try:
                #         # Initialize prop_C
                #         np.linalg.cholesky(prop_sigma['theta_c'])
                #         check_chol_cont = False
                #     except  np.linalg.LinAlgError:
                #         prop_sigma['theta_c'] = prop_sigma['theta_c'] + eps*np.eye(2)
                #         print("Oops. Proposal covariance matrix is now:\n")
                #         print(prop_sigma['theta_c'])
                       
                
                sigma_m['beta_loc0'] = np.exp(np.log(sigma_m['beta_loc0']) + gamma1*(beta_loc0_accept/thinning - r_opt_1d))
                beta_loc0_accept = 0      
                # sigma_m['beta_loc0'] = np.exp(np.log(sigma_m['beta_loc0']) + gamma1*(beta_loc0_accept/thinning - r_opt_2d))
                # beta_loc0_accept = 0
                # prop_sigma['beta_loc0'] = prop_sigma['beta_loc0'] + gamma2*(np.cov(beta_loc0_trace_within_thinning) - prop_sigma['beta_loc0'])
                # check_chol_cont = True
                # while check_chol_cont:
                #     try:
                #         # Initialize prop_C
                #         np.linalg.cholesky(prop_sigma['beta_loc0'])
                #         check_chol_cont = False
                #     except  np.linalg.LinAlgError:
                #         prop_sigma['beta_loc0'] = prop_sigma['beta_loc0'] + eps*np.eye(n_covariates)
                #         print("Oops. Proposal covariance matrix is now:\n")
                #         print(prop_sigma['beta_loc0'])
                
                # sigma_m['beta_loc1'] = np.exp(np.log(sigma_m['beta_loc1']) + gamma1*(beta_loc1_accept/thinning - r_opt_2d))
                # beta_loc1_accept = 0
                # prop_sigma['beta_loc1'] = prop_sigma['beta_loc1'] + gamma2*(np.cov(beta_loc1_trace_within_thinning) - prop_sigma['beta_loc1'])
                # check_chol_cont = True
                # while check_chol_cont:
                #     try:
                #         # Initialize prop_C
                #         np.linalg.cholesky(prop_sigma['beta_loc1'])
                #         check_chol_cont = False
                #     except  np.linalg.LinAlgError:
                #         prop_sigma['beta_loc1'] = prop_sigma['beta_loc1'] + eps*np.eye(n_covariates)
                #         print("Oops. Proposal covariance matrix is now:\n")
                #         print(prop_sigma['beta_loc1'])
                      
                sigma_m['beta_scale'] = np.exp(np.log(sigma_m['beta_scale']) + gamma1*(beta_scale_accept/thinning - r_opt_1d))
                beta_scale_accept = 0
                # sigma_m['beta_scale'] = np.exp(np.log(sigma_m['beta_scale']) + gamma1*(beta_scale_accept/thinning - r_opt_2d))
                # beta_scale_accept = 0
                # prop_sigma['beta_scale'] = prop_sigma['beta_scale'] + gamma2*(np.cov(beta_scale_trace_within_thinning) - prop_sigma['beta_scale'])
                # check_chol_cont = True
                # while check_chol_cont:
                #     try:
                #         # Initialize prop_C
                #         np.linalg.cholesky(prop_sigma['beta_scale'])
                #         check_chol_cont = False
                #     except  np.linalg.LinAlgError:
                #         prop_sigma['beta_scale'] = prop_sigma['beta_scale'] + eps*np.eye(n_covariates)
                #         print("Oops. Proposal covariance matrix is now:\n")
                #         print(prop_sigma['beta_scale'])
                     
                # sigma_m['beta_shape'] = np.exp(np.log(sigma_m['beta_shape']) + gamma1*(beta_shape_accept/thinning - r_opt_2d))
                # beta_shape_accept = 0
                sigma_m['beta_shape'] = np.exp(np.log(sigma_m['beta_shape']) + gamma1*(beta_shape_accept/thinning - r_opt_1d))
                beta_shape_accept = 0
                # prop_sigma['beta_shape'] = prop_sigma['beta_shape'] + gamma2*(np.cov(beta_shape_trace_within_thinning) - prop_sigma['beta_shape'])
                # check_chol_cont = True
                # while check_chol_cont:
                #     try:
                #         # Initialize prop_C
                #         np.linalg.cholesky(prop_sigma['beta_shape'])
                #         check_chol_cont = False
                #     except  np.linalg.LinAlgError:
                #         prop_sigma['beta_shape'] = prop_sigma['beta_shape'] + eps*np.eye(n_covariates)
                #         print("Oops. Proposal covariance matrix is now:\n")
                #         print(prop_sigma['beta_shape'])
              
        # ----------------------------------------------------------------------------------------                
        # -------------------------- Echo & save every 'thinning' steps --------------------------
        # ----------------------------------------------------------------------------------------
        if (iter / thinning) % echo_interval == 0:
            print(rank, iter, phi, tau_sqd)
            if rank == 0:
                print('Done with '+str(index)+" updates while thinned by "+str(thinning)+" steps,\n")
                   
                # Save the intermediate results to filename
                initial_values = {'phi':phi,
                     'gamma':gamma,
                     'tau_sqd':tau_sqd,
                     'prob_below':prob_below,
                     'prob_above':prob_above,
                     'Dist':Dist,
                     'theta_c':theta_c,
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
                with open(filename, 'wb') as f:
                    dump(Y, f)
                    dump(cen, f)
                    dump(cen_above,f)
                    dump(initial_values, f)
                    dump(sigma_m, f)
                    dump(prop_sigma, f)
                    dump(iter, f)
                    dump(phi_trace, f)
                    dump(tau_sqd_trace, f)
                    dump(theta_c_trace, f)
                    dump(beta_loc0_trace, f)
                    dump(beta_loc1_trace, f)
                    dump(beta_scale_trace, f)
                    dump(beta_shape_trace, f)
                       
                    dump(Z_1t_trace, f)
                    dump(R_1t_trace, f)
                    dump(Y_onetime, f)
                    dump(X_onetime, f)
                    dump(X_s_onetime, f)
                    dump(R_onetime, f)
                    dump(Z_onetime, f)
                    f.close()
                       
                # Echo trace plots
                pdf_pages = PdfPages('./'+save_directory+'/progress.pdf')
                grid_size = (4,2)
                #-page-1
                fig = plt.figure(figsize = (8.75, 11.75))
                plt.subplot2grid(grid_size, (0,0)) # phi
                plt.plot(phi_trace, color='gray', linestyle='solid')
                plt.ylabel(r'$\phi$')
                plt.subplot2grid(grid_size, (0,1)) # tau_sqd
                plt.plot(tau_sqd_trace, color='gray', linestyle='solid')
                plt.ylabel(r'$\tau^2$')
                plt.subplot2grid(grid_size, (1,0)) # rho
                plt.plot(theta_c_trace[0,:], color='gray', linestyle='solid')
                plt.ylabel(r'Matern $\rho$')
                plt.subplot2grid(grid_size, (1,1)) # nu
                plt.plot(theta_c_trace[1,:], color='gray', linestyle='solid')
                plt.ylabel(r'Matern $\nu$')
                plt.subplot2grid(grid_size, (2,0)) # mu0: beta_0
                plt.plot(beta_loc0_trace[0,:], color='gray', linestyle='solid')
                plt.ylabel(r'Location $\mu_0$: $\beta_0$')
                plt.subplot2grid(grid_size, (2,1)) # mu0: beta_1
                plt.plot(beta_loc0_trace[1,:], color='gray', linestyle='solid')
                plt.ylabel(r'Location $\mu_0$: $\beta_1$') 
                plt.subplot2grid(grid_size, (3,0)) # mu1: beta_0
                plt.plot(beta_loc1_trace[0,:], color='gray', linestyle='solid')
                plt.ylabel(r'Location $\mu_1$: $\beta_0$')
                plt.subplot2grid(grid_size, (3,1)) # mu1: beta_1
                plt.plot(beta_loc1_trace[1,:], color='gray', linestyle='solid')
                plt.ylabel(r'Location $\mu_1$: $\beta_1$')
                plt.tight_layout()
                pdf_pages.savefig(fig)
                plt.close()
                       
                #-page-2
                fig = plt.figure(figsize = (8.75, 11.75))
                plt.subplot2grid(grid_size, (0,0)) # scale: beta_0
                plt.plot(beta_scale_trace[0,:], color='gray', linestyle='solid')
                plt.ylabel(r'Scale $\sigma$: $\beta_0$')
                plt.subplot2grid(grid_size, (0,1)) # scale: beta_1
                plt.plot(beta_scale_trace[1,:], color='gray', linestyle='solid')
                plt.ylabel(r'Scale $\sigma$: $\beta_1$')
                plt.subplot2grid(grid_size, (1,0))  # shape: beta_0
                plt.plot(beta_shape_trace[0,:], color='gray', linestyle='solid')
                plt.ylabel(r'Shape $\xi$: $\beta_0$')
                plt.subplot2grid(grid_size, (1,1))  # shape: beta_1
                plt.plot(beta_shape_trace[1,:], color='gray', linestyle='solid')
                plt.ylabel(r'Shape $\xi$: $\beta_1$')
                plt.subplot2grid(grid_size, (2,0))   # X^*
                plt.plot(Z_1t_trace[1,:], color='gray', linestyle='solid')
                plt.ylabel(r'$Z$'+'['+str(1)+","+str(rank)+']')
                where = [(2,1),(3,0),(3,1)]
                for wh_sub,i in enumerate(wh_to_plot_Xs):
                    plt.subplot2grid(grid_size, where[wh_sub]) # X^*
                    plt.plot(Z_1t_trace[i,:], color='gray', linestyle='solid')
                    plt.ylabel(r'$Z$'+'['+str(i)+","+str(rank)+']')        
                plt.tight_layout()
                pdf_pages.savefig(fig)
                plt.close()
                pdf_pages.close()                 
            else:
                with open(filename, 'wb') as f:
                    dump(Y, f)
                    dump(cen, f)
                    dump(cen_above,f)
                    # dump(initial_values, f)
                    dump(sigma_m, f)
                    dump(iter, f)
                    dump(Z_1t_trace, f)
                    dump(R_1t_trace, f)
                    dump(Y_onetime, f)
                    dump(X_onetime, f)
                    dump(X_s_onetime, f)
                    dump(R_onetime, f)
                    dump(Z_onetime, f)
                    f.close()
                   
    
    
# import matplotlib.pyplot as plt
# def test(phi):
#     return utils.phi_update_mixture_me_likelihood(Y, phi, R, Z, cen, cen_above, prob_below, prob_above, Loc, Scale, Shape, 
#                         tau_sqd, gamma)

# Phi = np.arange(phi-0.01,phi+0.005,step=0.001)
# Lik = np.zeros(len(Phi))
# for idx, phi_tmp in enumerate(Phi):
#     Lik[idx] = test(phi_tmp)
# plt.plot(Phi, Lik, color='black', linestyle='solid')
# plt.axvline(phi, color='r', linestyle='--');

# # X_s = (R**phi)*utils.norm_to_Pareto(Z)
# def test(x):
#     return utils.loc0_gev_update_mixture_me_likelihood(Design_mat, np.array([x,beta_loc0[1]]), Y, X_s, cen, cen_above, prob_below, prob_above, 
#                       tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
#                       thresh_X, thresh_X_above)


# Coef = np.arange(beta_loc0[0]-0.5,beta_loc0[0]+1,step=0.01)
# Lik = np.zeros(len(Coef))
# for idx, coef in enumerate(Coef):
#     Lik[idx] = test(coef)
# plt.plot(Coef, Lik, color='black', linestyle='solid')
# plt.axvline(beta_loc0[0], color='r', linestyle='--');


# def test(x):
#     return utils.loc0_gev_update_mixture_me_likelihood(Design_mat, np.array([beta_loc0[0],x]), Y, X_s, cen, cen_above, prob_below, prob_above, 
#                       tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
#                       thresh_X, thresh_X_above)


# Coef = np.arange(beta_loc0[1]-0.5,beta_loc0[1]+1,step=0.01)
# Lik = np.zeros(len(Coef))
# for idx, coef in enumerate(Coef):
#     Lik[idx] = test(coef)
# plt.plot(Coef, Lik, color='black', linestyle='solid')
# plt.axvline(beta_loc0[1], color='r', linestyle='--');

# def test(x):
#     return utils.loc0_interc_gev_update_mixture_me_likelihood(Design_mat, x, beta_loc0[1], Y, X_s, cen, cen_above, prob_below, prob_above, 
#                       tau_sqd, phi, gamma, loc1, Scale, Shape, Time, xp, surv_p, den_p, 
#                       thresh_X, thresh_X_above)


# Coef = np.arange(beta_loc0[0]-0.01,beta_loc0[0]+0.01,step=0.001)
# Lik = np.zeros(len(Coef))
# for idx, coef in enumerate(Coef):
#     Lik[idx] = test(coef)
# plt.plot(Coef, Lik, color='black', linestyle='solid')
# plt.axvline(beta_loc0[0], color='r', linestyle='--');
   
# def test(x):
#     return utils.scale_interc_gev_update_mixture_me_likelihood(Design_mat, x, beta_scale[1], Y, X_s, cen, cen_above, prob_below, prob_above, 
#                       tau_sqd, phi, gamma, Loc, Shape, Time, xp, surv_p, den_p, 
#                       thresh_X, thresh_X_above)


# Coef = np.arange(beta_scale[0]-0.01,beta_scale[0]+0.01,step=0.001)
# Lik = np.zeros(len(Coef))
# for idx, coef in enumerate(Coef):
#     Lik[idx] = test(coef)
# plt.plot(Coef, Lik, color='black', linestyle='solid')
# plt.axvline(beta_scale[0], color='r', linestyle='--');

# def test(x):
#     return utils.shape_interc_gev_update_mixture_me_likelihood(Design_mat, x, beta_shape[1], Y, X_s, cen, cen_above, prob_below, prob_above, 
#                       tau_sqd, phi, gamma, Loc, Scale, Time, xp, surv_p, den_p, 
#                       thresh_X, thresh_X_above)

# # def test(x):
# #     return np.sum(utils.dgev(Y[~cen & ~cen_above],20,1,x,log=True))

# Coef = np.arange(beta_shape[0]-0.01,beta_shape[0]+0.01,step=0.001)
# Lik = np.zeros(len(Coef))
# for idx, coef in enumerate(Coef):
#     Lik[idx] = test(coef)
# plt.plot(Coef, Lik, color='black', linestyle='solid')
# plt.axvline(beta_shape[0], color='r', linestyle='--');

     
# X_s = (R**phi)*utils.norm_to_Pareto(Z) 
# def test(tau_sqd):
#     return utils.tau_update_mixture_me_likelihood(Y, tau_sqd, X_s, cen, cen_above, 
#                     prob_below, prob_above, Loc, Scale, Shape, 
#                     phi, gamma, xp, surv_p, den_p)

# Tau = np.arange(1000,10000,step=90)
# Lik = np.zeros(len(Tau))
# for idx, t in enumerate(Tau):
#     Lik[idx] = test(t) 
# plt.plot(Tau, Lik, color='black', linestyle='solid')
# plt.axvline(tau_sqd, color='r', linestyle='--');   
    
# def test(x):
#     return utils.theta_c_update_mixture_me_likelihood(Z, np.array([x,1.5]), Dist)

# Range = np.arange(range-0.02,range+0.04,step=0.004)
# Lik = np.zeros(len(Range))
# for idx, r in enumerate(Range):
#     Lik[idx] = test(r) 
# plt.plot(Range, Lik, color='black', linestyle='solid')
# plt.axvline(range, color='r', linestyle='--');

# def test(x):
#     return utils.theta_c_update_mixture_me_likelihood(Z, np.array([range,x]), Dist)

# Nu = np.arange(0.1,1.8,step=0.01)
# Lik = np.zeros(len(Nu))
# for idx, r in enumerate(Nu):
#     Lik[idx] = test(r) 
# plt.plot(Nu, Lik, color='black', linestyle='solid')
# plt.axvline(nu, color='r', linestyle='--');

# t_chosen = 3
# def test(x):
#     return utils.Rt_update_mixture_me_likelihood(Y[:,t_chosen], x, X[:,t_chosen], Z[:,t_chosen], cen[:,t_chosen], cen_above[:,t_chosen], 
#                 prob_below, prob_above, Loc[:,t_chosen], Scale[:,t_chosen], Shape[:,t_chosen], tau_sqd, phi, gamma,
#                 xp, surv_p, den_p, thresh_X, thresh_X_above) + priors.R_prior(x, gamma)

# Rt = np.arange(R[t_chosen]-0.1,R[t_chosen]+0.1,step=0.001)
# Lik = np.zeros(len(Rt))
# for idx, r in enumerate(Rt):
#     Lik[idx] = test(r) 
# plt.plot(Rt, Lik, linestyle='solid')
# plt.axvline(R[t_chosen], color='r', linestyle='--');

from pickle import load
with open('local_0_progress_0.pkl', 'rb') as f:
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
                   
      Z_1t_trace=load(f)
      R_1t_trace=load(f)
      Y_onetime=load(f)
      X_onetime=load(f)
      X_s_onetime=load(f)
      R_onetime=load(f)
      Z_onetime=load(f)
    
plt.plot(R_1t_trace[:],linestyle='solid')
plt.hlines(data_all['R_at_knots'][0,0], 0, n_updates_thinned, colors='r', linestyles='--');

plt.plot(beta_loc0_trace_tmp[0,1000:],linestyle='solid')

# from pickle import load
# with open('local_0_progress_3.pkl', 'rb') as f:
#      Y_tmp=load(f)
#      cen_tmp=load(f)
#      cen_above_tmp=load(f)
#      # initial_values_tmp=load(f)
#      sigma_m=load(f)
#      # prop_sigma=load(f)
#      iter_tmp=load(f)
#      Z_1t_trace=load(f)
#      R_1t_trace=load(f)
#      Y_onetime=load(f)
#      X_onetime=load(f)
#      X_s_onetime=load(f)
#      R_onetime=load(f)
#      Z_onetime=load(f)
    
# plt.plot(R_1t_trace[:],linestyle='solid')
# plt.hlines(data_all['R_at_knots'][0,3], 0, n_updates_thinned, colors='r', linestyles='--');