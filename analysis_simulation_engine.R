library(GPfit)
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
#library(LaplacesDemon)
library(mvnfast)
library(RColorBrewer)
library(pheatmap)
library(philentropy)

peaks_significance_threshold=0.3

mean_func <- function(x) {
  rep(0, length(x))
}

sample_prior <- function(
    x, mean_func, cov_func, ..., random_seed = -1, n_samples = 5
) {
  x_mean <- mean_func(x)
  x_cov <- cov_func(get_r(x, x), ...)
  random_seed <- as.integer(random_seed)
  if (random_seed > 0) set.seed(random_seed)
  return(MASS::mvrnorm(n_sample, x_mean, x_cov))
}



sample_posterior <- function(X, 
                             y, 
                             X_star, 
                             mean_func, 
                             sigma_n = 5, 
                             random_seed = -1, 
                             n_samples = 5,
                             K_xsxs = NA,
                             plot_kernel=T) 
  {

  # k_xx = cov_func(get_r(x, x), l,v)
  # k_xxs = cov_func(get_r(x, x_star), l,v)
  # k_xsx = cov_func(get_r(x_star, x), l,v)
  # k_xsxs = cov_func(get_r(x_star, x_star), l,v)

  # 
  # 
  K_xx <- rbf_kernel(X, sigma=sigma_n)
  K_xxs <- rbf_kernel(X, Xtag=X_star, sigma=sigma_n)
  K_xsx <-  rbf_kernel(X_star, Xtag=X, sigma=sigma_n)
  
  I = diag(1, dim(K_xx)[1])
  K_xx_noise = solve(K_xx + sigma_n ^ 2 * I)
  Kxsx_kxxNoise = K_xsx %*% K_xx_noise
  

  mu_conditioned = Kxsx_kxxNoise %*% y
  cov_mat_conditioned = K_xsxs - Kxsx_kxxNoise %*% K_xxs
  
  if (plot_kernel) {
    pheatmap(cov_mat_conditioned, cluster_cols=F, cluster_rows=F) #col=mako(100))
  }
  
  # I = diag(1, dim(k_xx)[1])
  # k_xx_noise = solve(k_xx + sigma_n ^ 2 * I)
  # kxsx_kxxNoise = k_xsx %*% k_xx_noise
  # # Eq.2.23, 24
  # fsb = kxsx_kxxNoise %*% y
  # cov_fs = k_xsxs - kxsx_kxxNoise %*% k_xxs
  
  A <- matrix(0, n_samples, nrow(X_star))
  rmvn(n_samples, mu = mu_conditioned, cov_mat_conditioned, A=A)
  return(A)
  # random_seed <- as.integer(random_seed)
  
  
  # if (random_seed > 0) set.seed(random_seed)
  # return(MASS::mvrnorm(n_samples, mu_conditioned, cov_mat_conditioned))
}



# par(mfrow = c(3, 3))
#   
# x1_train <- c(-10, -8, -5, 3, 7, 15)
# x2_trains <- c(1, 3, 4, 6, 9, 0)
# y_train <- c(-5, -2, 3, 3, 2, 5)
# x1_star <- seq(-20, 20, length = 200)
# x2_star <- seq(-8,8, length=200)
# 
# X_train <- cbind(sort(runif(50, -5,5 )), 
#                  runif(50, -5,5))
# y_train <- apply(X_train, 1 , function(r) {test_objective(r[1], r[2])})
# 
# 
#   dt <- sample_posterior(X_train, 
#                          y_train, 
#                          X_star, 
#                          mean_func, 
#                          n_samples = 5,
#                          sigma_n = 0.5,
#                          random_seed = -1)
#   
  

euclidean <- function(a, b) sqrt(sum((a - b)^2))

rbf_kernel <- function(X, Xtag = NA, sigma=1) {
  
  if (!all(is.na(Xtag))) {
    rbf_kernel_matrix <- matrix(rep(0, times=nrow(X) * nrow(Xtag)), nrow = nrow(X), ncol=nrow(Xtag))
    
    for (i in 1:nrow(X)) {
      for(j in 1:nrow(Xtag)) {
        rbf_kernel_matrix[i,j] <- exp(-((euclidean(X[i,], Xtag[j,]))/(2 * sigma)))
      }
    }    
    
  } else {
    rbf_kernel_matrix <- matrix(rep(0, times=nrow(X) ** 2), nrow = nrow(X))
    
    for (i in 1:nrow(X)) {
      for(j in 1:nrow(X)) {
        rbf_kernel_matrix[i,j] <- exp(-((euclidean(X[i,], X[j,]))/(2 * sigma)))
      }
    }    
  }
  
  return(rbf_kernel_matrix)
}

test_objective <- function(x1,x2) {
  return (x1 ** 2 + x2 ** 2)
}

plot_gaussian_process_prediction <- function(res, X_star, x_train, y_train, improval_step=6) { 
  
  plot(res[1,], col="gray60", lty=2, type="l")
  
  for (i in 2:nrow(res)) {
    lines(res[i,], col="gray60", lty=2, type="l")
  }
  
  x_train_idx <-     
    apply(X_train, 1 , function(p1) {
      return(which.min(apply(X_star, 1, function(p2) {euclidean(as.numeric(p1), as.numeric(p2))})))  
    })
  
  points(x_train_idx, y_train, col="black", type="p", pch=19)
}

# Dummy function to estimate
f <- function(x) {
  return((6 * x - 2)^2 * sin(12 * x - 4))
}

# objective
bayesian_optimizer <- function(objective, 
                               parameter_set, 
                               cov_mat_params,
                               n0, 
                               max_iterations=10, 
                               improval_step=5,
                               sigma_n=0.35,
                               n_samples=20,
                               aquisition="prob",
                               true_spike_train,
                               stim_trace,
                               percentage,
                               use_SI,
                               remove_edges=F,
                               dt=0.05,
                               dist_fun="JSD") {
  
  
  print("##### FITTING SIMULATION PARAMETERS ON SPIKE TRAIN #####")
  print(sprintf("Spike train dimensions: neurons (%d) X frames (%d)", nrow(true_spike_train), ncol(true_spike_train)))
  print(sprintf("Stimulus trace frames: (%d)", len(stim_trace)))
  print(sprintf("Iterations: %d, Initial points (n0): %d, RBF Kernel sigma: %f, True percentage: %f",
                max_iterations,
                n0,
                sigma_n,
                percentage))
  print(sprintf("Samples drawn per predicted point (n_samples): %d, Point aquisition: '%s', Using SI %d", 
                n_samples, 
                aquisition, 
                as.numeric(use_SI)))
  print(sprintf("Parameter space dimensions: points (%d) X dimensions (%d)", nrow(parameter_set), ncol(parameter_set)))
  print(sprintf("Using %d Hz", 1/dt))
  
  
  cov_mat_pred <- cov_mat_params
  X_star <- parameter_set
  X_train <- parameter_set[sample(1:nrow(parameter_set), n0),]
  
  print(X_train)
  print ("###################")
  
  y <- objective(X_train, 
                 true_spike_train = true_spike_train,
                 stim_trace = stim_trace,
                 place_cell_pct = percentage,
                 use_SI = use_SI,
                 remove_edges=remove_edges,
                 dt=dt,
                 dist_fun = dist_fun)
  
  
  for (iter in 1:max_iterations){
    res <- sample_posterior(as.matrix(X_train),
                            y,
                            as.matrix(X_star),
                            mean_func,
                            sigma_n=sigma_n,
                            n_sample=n_samples,
                            K_xsxs = cov_mat_pred,
                            plot_kernel = F)
    
    
    
    #fl <- sprintf("C:\\Users\\itayta\\Desktop\\allcells_figures\\supp_figure_simulations\\for_movie\\iter_%d.csv", iter)
    #write.csv(file=fl, res)
    
    min_per_x <- apply(res,2,min)
    max_per_x <- apply(res,2,max)
    
    mean_per_pred <- apply(res,2,mean)
    sd_per_pred <- apply(res,2,sd)
    # plot(x=X_train, y, col="black", pch=19, type="p", ylim=c(-5,5),
    #      xlab="X to predict", ylab="Predicted Y")
    # polygon(c(X_star, rev(X_star)), c(mean_per_pred - sd_per_pred, rev(mean_per_pred + sd_per_pred)), col="gray50",border = "gray50")  
    # lines(x=X_train, y, col="black", pch=19, type="p", ylim=c(-5,5))
    # lines(x=X_star, f(X_star), col="red")
    # lines(x=X_star, mean_per_pred, col="green")
    
    gamma <- (min(y) - mean_per_pred)/sd_per_pred
    prob_to_improve <- pnorm(gamma)
    top_to_improve <- order(prob_to_improve, decreasing=T)
    
    
    # plot(y=prob_to_improve, x=X_star, type="l", xlab="X to predict", ylab="Probability to improve minimum", lwd=2)
    # points(X_star[top_to_improve], prob_to_improve[top_to_improve], col="red", pch=19)
    
    expected_improve <- sd_per_pred * (gamma * pnorm(gamma) + dnorm(gamma))
    top_to_improve_exp <- order(expected_improve, decreasing=T)
    # plot(y=expected_improve, x=X_star, type="l", xlab="X to predict", ylab="Expected improvement", lwd=2)
    # points(X_star[top_to_improve_exp], expected_improve[top_to_improve_exp], col="purple", pch=19)
    
    
    addition_idx = 1
    left_to_add = improval_step
    new_X <- c()
    
    while (left_to_add > 0) {
      
      addition_idx <- addition_idx + 1
      
      if (sum(apply(X_train, 1, function(r) {all(r == X_star[top_to_improve[addition_idx],])})) > 0) {
        next
      }
      
      if (aquisition == "prob") {
        new_X <- rbind(new_X, X_star[top_to_improve[addition_idx],])
      } else {
        new_X <- rbind(new_X, X_star[top_to_improve_exp[addition_idx],])
      }
      
      left_to_add <- left_to_add - 1
    }
    
    
    # Add new points
    X_train <- rbind(X_train, new_X)
    y <- c(y, objective(new_X, 
                        true_spike_train = true_spike_train,
                        stim_trace = stim_trace,
                        place_cell_pct = percentage,
                        use_SI = use_SI,
                        remove_edges=remove_edges,
                        dt=dt,
                        dist_fun=dist_fun))
    
    
    print(sprintf("------- Added %d new points (%d Total)", nrow(new_X), nrow(X_train)))
    
  }
  
  res <- sample_posterior(as.matrix(X_train),
                          y,
                          as.matrix(X_star),
                          mean_func,
                          sigma_n=sigma_n,
                          n_sample=n_samples,
                          K_xsxs = cov_mat_pred,
                          plot_kernel = F)
  
  
  #fl <- sprintf("C:\\Users\\itayta\\Desktop\\allcells_figures\\supp_figure_simulations\\for_movie\\iter_final.csv")
  #write.csv(fl, res)
  
  return(res)
}


get_peaks <- function(tuning_curve, threshold_from_max=NA) {
  
  tc <- c(0,tuning_curve,0)
  d <- diff(tc); d[d<0] <- -1; d[d>0] <- 1
  peaks <- which(diff(d) == -2)
  
  
  if (!is.na(threshold_from_max)) {
    peak_values <- tuning_curve[peaks]
    return(peaks[which(peak_values > threshold_from_max * max(peak_values))])
  }
  
  return(peaks)
}

generate_tuning_curves_cost <- function(n, 
                                        percentage, 
                                        plot=T, 
                                        average_width=0.14,
                                        sd_width=0.05,
                                        fixed_fr=NA, 
                                        noise=0.015, 
                                        double_peak_pct=0, 
                                        triple_peak_pct=0,
                                        average_width_3=0.07,
                                        average_width_2=0.07,
                                        n_bins = 24,
                                        first_bin=1,
                                        return_non_place=F) {
  
  n_neurons = n
  mu_rates=0.08
  std_rates=0.02
  
  noise_level=noise
  average_sigma=average_width
  std_sigma_f=sd_width
  average_sigma2=average_width_2
  average_sigma3=average_width_3 
  
  pref_positions = seq(0.5, n_bins + 0.5, length.out=n_neurons)
  pref_positions_2 = seq(0.5, n_bins + 0.5, length.out=n_neurons)
  pref_positions_2 <- pref_positions_2[sample(1:n_neurons, n_neurons)]
  
  pref_positions_3 = seq(0.5, n_bins + 0.5, length.out=n_neurons)
  pref_positions_3 <- pref_positions_3[sample(1:n_neurons, n_neurons)]
  average_rates=fixed_fr[sample(1:length(fixed_fr), len(fixed_fr))]
  
  
  if (plot) {
    par(mfrow=c(3,2))
    hist(average_rates, breaks=50)
  }
  
  
  mu_sigma=log((average_sigma ** 2) / sqrt(std_sigma_f + average_sigma ** 2))
  std_sigma=sqrt(log(std_sigma_f / (average_sigma ** 2) + 1));
  normed_sigma=rlnorm(n_neurons, meanlog=mu_sigma, sdlog=std_sigma);
  sigma=normed_sigma*n_bins
  
  mu_sigma2=log((average_sigma2 ** 2) / sqrt(std_sigma_f + average_sigma2 ** 2))
  std_sigma2=sqrt(log(std_sigma_f / (average_sigma2 ** 2) + 1));
  normed_sigma2=rlnorm(n_neurons, meanlog=mu_sigma2, sdlog=std_sigma2);
  sigma2=normed_sigma*n_bins
  
  mu_sigma3=log((average_sigma3 ** 2) / sqrt(std_sigma_f + average_sigma3 ** 2))
  std_sigma3=sqrt(log(std_sigma_f / (average_sigma3 ** 2) + 1));
  normed_sigma3=rlnorm(n_neurons, meanlog=mu_sigma3, sdlog=std_sigma3);
  sigma3=normed_sigma3*n_bins
  
  get_tuning <- function(n) {return((1/sqrt(2*pi*sigma[n])) * exp((-(1:n_bins - pref_positions[n]) ** 2) / (2 * sigma[n] ** 2)))}
  get_tuning2 <- function(n) {return((1/sqrt(2*pi*sigma[n])) * (exp((-(1:n_bins - pref_positions[n]) ** 2) / (2 * sigma[n] ** 2)) +
                                                                  runif(1,0.1,1) * exp((-(1:n_bins - pref_positions_2[n]) ** 2) / (2 * sigma2[n] ** 2))))}
  
  get_tuning3 <- function(n) {return((1/sqrt(2*pi*sigma[n])) * (exp((-(1:n_bins - pref_positions[n]) ** 2) / (2 * sigma[n] ** 2)) +
                                                                  runif(1,0.4,0.8) * exp((-(1:n_bins - pref_positions_2[n]) ** 2) / (2 * sigma2[n] ** 2)) + 
                                                                  runif(1,0.4,0.8) * exp((-(1:n_bins - pref_positions_3[n]) ** 2) / (2 * sigma3[n] ** 2))))}
  
  
  double_peak_ind <- -1
  triple_peak_ind <- -1
  
  if (double_peak_pct > 0) {
    double_peaked <- floor(rnorm(1, mean=n_neurons * double_peak_pct, sd=3))
    
    while(double_peaked < 1) {
      print("resampling")
      double_peaked <- floor(rnorm(1, mean=n_neurons * double_peak_pct, sd=3))
    }
    
    double_peak_ind <- sample(1:n_neurons, double_peaked)
  }
  
  if(triple_peak_pct > 0) {
    triple_peaked <- floor(rnorm(1, mean=n_neurons * triple_peak_pct, sd=1))
    
    while(triple_peaked < 1) {
      print("resampling")
      triple_peaked <- floor(rnorm(1, mean=n_neurons * triple_peak_pct, sd=1))
    }
    
    indices_to_sample <- 1:n_neurons
    
    if (!all(double_peak_ind == -1)) {
      indices_to_sample <- which(!1:n_neurons %in% double_peak_ind)
    }
    
    triple_peak_ind <- sample(indices_to_sample, triple_peaked)
  }
  
  
  
  tuning_curves <- 
    lapply(1:n_neurons, 
           function(n) {
             
             if ((!all(double_peak_ind == -1)) && (n %in% double_peak_ind)) {
               tuning <- get_tuning2(n)
             }  else if ((!all(triple_peak_ind == -1)) && (n %in% triple_peak_ind)) {
               tuning <- get_tuning3(n)
             } else {
               tuning <- get_tuning(n)  
             }
             
             factor <- average_rates[n]/(mean(tuning) + 10^-30)
             tuning <- tuning * factor + abs(rnorm(n_bins) * noise_level)
             return(tuning)
           })
  
  
  
  if (percentage != 1) {
    n_non_place_neurons <- n_neurons - floor(rnorm(1 ,mean=n_neurons * percentage, sd=10))
    
    while(n_non_place_neurons < 0 || n_non_place_neurons > n_neurons) {
      n_non_place_neurons <- n_neurons - floor(rnorm(1 ,mean=n_neurons * percentage, sd=10))
    }
    
    # Randomly draw n_non neurons and set them with random tuning curve
    non_place = sample(1:n_neurons, n_non_place_neurons)
    for (i in non_place) {
      #tuning_curves[[i]] <-  abs(rnorm(n_bins) * runif(1,0.01, 1.5) * average_rates[i])
      #tuning_curves[[i]] <- runif(n_bins, 0.8, 1) * average_rates[i]
      tuning_curves[[i]] <- rep(1, times=n_bins) + abs(rnorm(n_bins) * noise_level)
    }
  }
  
  if (plot) {
    barplot(tuning_curves[[runif(1,1,n_neurons)]])
    barplot(tuning_curves[[runif(1,1,n_neurons)]])
    
    if ((!all(double_peak_ind == -1))) {
      barplot(tuning_curves[[sample(double_peak_ind, 1)]], col="purple")
    }
    
    if ((!all(triple_peak_ind == -1))) {
      barplot(tuning_curves[[sample(triple_peak_ind, 1)]], col="orange")
    }
  }
  
  # assert(all(unlist(lapply(sample(1:400, 50), function(i) { all(tuning_curves[[i]] == mt[i,])}))))
  tc2 <- lapply(tuning_curves, function(r) {r / sum(r)})
  mixed_rates <- fixed_fr[sample(1:length(fixed_fr), len(fixed_fr))]
  tc3 <- lapply(1:length(tuning_curves), function(i) {tc2[[i]] * mixed_rates[i]})
  mt <- t(matrix(unlist(tc3), nrow = n_bins))
  
  if (return_non_place && percentage != 1) {
    return(list(mt, non_place))
  }
  
  return(mt)
  
}
currate_spike_train_cost <- function(tc, true_tb, stim_trace, by_pval_only=F,verbose=F, jump=0.05) {
  
  if (dt == 0.05) {
    fact <- seq(.05, 3.6, by=jump)
  } else {
    fact <- seq(0.4, 3.6, by=jump)
    
  }
  
  
  
  pvals <- c()
  lengths <- c()
  for (f in fact) {
    gst <- generate_spike_trains_cost(tc, stim_trace, fs=1, factor=f)
    sim <- preprocess_spike_train(gst, stim_trace, verbose=F)
    pvals <- c(pvals, ks.test(sim$working_time_bins_per_cells, true_tb)$p.value)
    lengths <- c(lengths, (len(sim$working_time_bins_per_cells) / len(true_tb)))
    if (verbose) {
      print(f)
      print(ks.test(sim$working_time_bins_per_cells, true_tb)$p.value)
      print(len(sim$working_time_bins_per_cells) / len(true_tb))
    }
  }
  
  o1 <- order(pvals, decreasing = T)
  o2 <- order((lengths - 1) ** 2)
  
  rank <- rep(0, times=len(fact))
  names(rank) <- fact
  
  rank[o1] <- rank[o1] + 1:len(fact)
  rank[o2] <- rank[o2] + 1:len(fact)
  
  if (by_pval_only) {
    return(fact[which.max(pvals)])
  }
  
  return(fact[which.min(rank)])
}


get_peaks <- function(tuning_curve, threshold_from_max=NA) {
  
  tc <- c(0,tuning_curve,0)
  d <- diff(tc); d[d<0] <- -1; d[d>0] <- 1
  peaks <- which(diff(d) == -2)
  
  
  if (!is.na(threshold_from_max)) {
    peak_values <- tuning_curve[peaks]
    return(peaks[which(peak_values > threshold_from_max * max(peak_values))])
  }
  
  return(peaks)
}


generate_spike_trains_cost <- function(tuning_curves, stim_trace, fs, plot=F, factor=1) {
  
  st <- t(apply(tuning_curves, 1 , function(tc) {rpois(len(stim_trace), tc[stim_trace] / (fs * factor))}))
  
  if (plot) {
    a <- preprocess_spike_train(st, stim_trace, verbose=F)
    hist(a$working_time_bins_per_cells, breaks=50)
  }
  
  return(st)
}

spike_train_similiarity_cost_new <- 
  function(true_spike_train, 
           stim_trace,
           distance_function="DKL",
           average_width,
           sd_width,
           noise_level,
           double_peak_percent,
           triple_peak_percent,
           place_cell_percentage,
           dt=0.05,
           n_reps=20,
           pois_factor_reps=0,
           remove_edges=F,
           use_SI=F) {
    
    print(sprintf("Calculating cost for mu(%f), sd(%f), noise(%f), double_p(%f), triple_p(%f), place_pct(%f)",
                  average_width,
                  sd_width,
                  noise_level,
                  double_peak_percent,
                  triple_peak_percent,
                  place_cell_percentage))
    
    additive_cost <- c()
    
    firing_rate <-rowMeans(true_spike_train) / dt
    processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
    true_cells_spike_train <- processed_real$working_cells_spike_train
    true_firing_rate <- processed_real$working_firing_rate
    true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
    
    tmp <- compute_tuning(true_cells_spike_train, stim_trace)
    stim_prob <- tmp[[1]]
    true_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    smoothed_tuning <- 
      t(apply(true_tuning_curve, 1, function(trc) {
        barp <- barplot(trc, plot=F)
        smt <- smooth.spline(barp, trc, all.knots = T, lambda=1e-5)
        smt$y[smt$y < 0] <- 0
        return(smt$y);
      }))
    
    true_peaks <- unlist(apply(smoothed_tuning, 1, 
                               function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
    true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
    
    final_pois_factor = 0.6
    pois_calc_reps_left = pois_factor_reps
    
    for (nrep in 1:n_reps) {
      tuning_curve <-
        generate_tuning_curves_cost(n = nrow(true_spike_train),
                                    percentage = place_cell_percentage,
                                    average_width = average_width, 
                                    sd_width = sd_width,
                                    fixed_fr=firing_rate,
                                    noise=noise_level,
                                    n_bins=ifelse(remove_edges, 20, 24),
                                    double_peak_pct = double_peak_percent,
                                    triple_peak_pct = triple_peak_percent,
                                    plot=F)
      
      if (pois_calc_reps_left > 0){
        pois_factor <- currate_spike_train_cost(tuning_curve, 
                                                true_time_bins_per_cells,
                                                stim_trace, verbose=T)
        final_pois_factor = final_pois_factor + pois_factor
        pois_calc_reps_left = pois_calc_reps_left - 1
        
        if (pois_calc_reps_left == 0){
          final_pois_factor = final_pois_factor / pois_factor_reps
          #print(sprintf("Using final pois factor %f", final_pois_factor))
        }
      } else {
        pois_factor = final_pois_factor
      }
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)#fs=(1/dt))
      
      
      ### Generate a spike train
      processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
      generated_cells_active_st <- processed_generated$working_cells_spike_train
      generated_firing_rate <- processed_generated$working_firing_rate
      generated_time_bins <- processed_generated$working_time_bins_per_cells
      
      
      tmp <- compute_tuning(generated_cells_active_st, stim_trace)
      gen_stim_prob <- tmp[[1]]
      gen_tuning_curve <- tmp[[2]]
      
      gen_smoothed_tuning <- 
        t(apply(gen_tuning_curve, 1, function(trc) {
          barp <- barplot(trc, plot=F)
          smt <- smooth.spline(barp, trc, all.knots = T, lambda=1e-5)
          smt$y[smt$y < 0] <- 0
          return(smt$y);
        }))
      
      rm(tmp)
      
      generated_peaks <- unlist(apply(gen_smoothed_tuning, 1, 
                                      function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
      generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
      generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
      
      
      
      if (distance_function == "KS") {
        cost <- KS_distance(generated_SI = generated_SI,
                            generated_time_bins = generated_time_bins,
                            generated_peaks = generated_peaks,
                            generated_active_spatial_bins_per_cell = generated_active_spatial_bins_per_cell,
                            true_SI = true_SI,
                            true_time_bins_per_cells = true_time_bins_per_cells,
                            true_peaks = true_peaks,
                            true_active_spatial_bins_per_cell = true_active_spatial_bins_per_cell,
                            use_SI = use_SI)
        
      }  else if (distance_function == "JSD") {
        
        cost <- JSD_distance(generated_SI = generated_SI,
                             generated_time_bins = generated_time_bins,
                             generated_peaks = generated_peaks,
                             generated_active_spatial_bins_per_cell = generated_active_spatial_bins_per_cell,
                             true_SI = true_SI,
                             true_time_bins_per_cells = true_time_bins_per_cells,
                             true_peaks = true_peaks,
                             true_active_spatial_bins_per_cell = true_active_spatial_bins_per_cell,
                             use_SI = use_SI)
      } else {
        cost <- DKL_distance(generated_SI = generated_SI,
                             generated_time_bins = generated_time_bins,
                             generated_peaks = generated_peaks,
                             generated_active_spatial_bins_per_cell = generated_active_spatial_bins_per_cell,
                             true_SI = true_SI,
                             true_time_bins_per_cells = true_time_bins_per_cells,
                             true_peaks = true_peaks,
                             true_active_spatial_bins_per_cell = true_active_spatial_bins_per_cell,
                             use_SI = use_SI)  
      }
      
      #print(sprintf("Cost is: %f", cost))
      
      additive_cost <- c(additive_cost , cost)
    }
    
    
    sum(additive_cost)
    return(additive_cost)
  }

spike_train_similiarity_cost <- function(true_spike_train, 
                                         stim_trace,
                                         distance_function="DKL",
                                         average_width,
                                         sd_width,
                                         noise_level,
                                         double_peak_percent,
                                         place_cell_percentage,
                                         remove_edges=F,
                                         dt=0.05,
                                         n_reps=20,
                                         pois_factor_reps=1,
                                         use_SI=F) {
  
  print(sprintf("Calculating cost for mu(%f), sd(%f), noise(%f), double_p(%f), place_pct(%f)",
                average_width,
                sd_width,
                noise_level,
                double_peak_percent,
                place_cell_percentage))
  
  additive_cost <- c()
  
  firing_rate <-rowMeans(true_spike_train) / dt
  processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(true_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  true_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  true_peaks <- unlist(apply(true_tuning_curve, 1, 
                             function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
  true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
  true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
  
  final_pois_factor = 0
  pois_calc_reps_left = pois_factor_reps
  
  for (nrep in 1:n_reps) {
    tuning_curve <-
      generate_tuning_curves_cost(n = nrow(true_spike_train),
                                  percentage = place_cell_percentage,
                                  average_width = average_width, 
                                  sd_width = sd_width,
                                  fixed_fr=firing_rate,
                                  noise=noise_level,
                                  n_bins=ifelse(remove_edges, 20, 24),
                                  double_peak_pct = double_peak_percent,
                                  plot=F)
    
    if (pois_calc_reps_left > 0){
      pois_factor <- currate_spike_train_cost(tuning_curve, 
                                              true_time_bins_per_cells,
                                              stim_trace, verbose=F)
      final_pois_factor = final_pois_factor + pois_factor
      pois_calc_reps_left = pois_calc_reps_left - 1
      
      if (pois_calc_reps_left == 0){
        final_pois_factor = final_pois_factor / pois_factor_reps
        #print(sprintf("Using final pois factor %f", final_pois_factor))
      }
    } else {
      pois_factor = final_pois_factor
    }
    generated_spike_train <- 
      generate_spike_trains_cost(tuning_curves = tuning_curve,
                                 stim_trace = stim_trace,
                                 factor=pois_factor,
                                 fs=1)#fs=(1/dt))
    
    
    ### Generate a spike train
    processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
    generated_cells_active_st <- processed_generated$working_cells_spike_train
    generated_firing_rate <- processed_generated$working_firing_rate
    generated_time_bins <- processed_generated$working_time_bins_per_cells
    
    
    tmp <- compute_tuning(generated_cells_active_st, stim_trace)
    gen_stim_prob <- tmp[[1]]
    gen_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    generated_peaks <- unlist(apply(gen_tuning_curve, 1, 
                                    function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
    generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
    
    
    
    if (distance_function == "KS") {
      cost <- KS_distance(generated_SI = generated_SI,
                          generated_time_bins = generated_time_bins,
                          generated_peaks = generated_peaks,
                          generated_active_spatial_bins_per_cell = generated_active_spatial_bins_per_cell,
                          true_SI = true_SI,
                          true_time_bins_per_cells = true_time_bins_per_cells,
                          true_peaks = true_peaks,
                          true_active_spatial_bins_per_cell = true_active_spatial_bins_per_cell,
                          use_SI = use_SI)
      
    }  else if (distance_function == "JSD") {
      
      cost <- JSD_distance(generated_SI = generated_SI,
                           generated_time_bins = generated_time_bins,
                           generated_peaks = generated_peaks,
                           generated_active_spatial_bins_per_cell = generated_active_spatial_bins_per_cell,
                           true_SI = true_SI,
                           true_time_bins_per_cells = true_time_bins_per_cells,
                           true_peaks = true_peaks,
                           true_active_spatial_bins_per_cell = true_active_spatial_bins_per_cell,
                           use_SI = use_SI)
    } else {
      cost <- DKL_distance(generated_SI = generated_SI,
                           generated_time_bins = generated_time_bins,
                           generated_peaks = generated_peaks,
                           generated_active_spatial_bins_per_cell = generated_active_spatial_bins_per_cell,
                           true_SI = true_SI,
                           true_time_bins_per_cells = true_time_bins_per_cells,
                           true_peaks = true_peaks,
                           true_active_spatial_bins_per_cell = true_active_spatial_bins_per_cell,
                           use_SI = use_SI)  
    }
    
    #print(sprintf("Cost is: %f", cost))
    
    additive_cost <- c(additive_cost , cost)
  }
  
  
  sum(additive_cost)
  return(additive_cost)
}


DKL_distance <- function(generated_SI,
                         generated_time_bins,
                         generated_peaks,
                         generated_active_spatial_bins_per_cell,
                         true_SI,
                         true_time_bins_per_cells,
                         true_peaks,
                         true_active_spatial_bins_per_cell,
                         use_SI) {
  # Compare SI distribution distances
  max_SI <- max(c(generated_SI[[1]], true_SI[[1]]))
  true_si_hist <- hist(true_SI[[1]], breaks=seq(0, max_SI, length.out=100), plot=F)
  true_si_prob_dist <- true_si_hist$density / sum(true_si_hist$density)
  
  gen_si_hist <- hist(generated_SI[[1]], breaks=seq(0, max_SI, length.out=100), plot=F)
  gen_si_prob_dist <- gen_si_hist$density / sum(gen_si_hist$density)
  
  si_kld <- KLD(true_si_prob_dist, gen_si_prob_dist)
  
  
  # Compare time bins distribution distances
  max_time_bins <- max(generated_time_bins, true_time_bins_per_cells)
  true_time_bins_hist <- hist(true_time_bins_per_cells, breaks=seq(0, max_time_bins, length.out=100), plot=F)
  true_time_bins_prob_dist <- true_time_bins_hist$density / sum(true_time_bins_hist$density)
  
  gen_time_bins_hist <- hist(generated_time_bins, breaks=seq(0, max_time_bins, length.out=100), plot=F)
  gen_time_bins_prob_dist <- gen_time_bins_hist$density / sum(gen_time_bins_hist$density)
  
  time_bins_kld <- KLD(true_time_bins_prob_dist, gen_time_bins_prob_dist)
  
  # Compare active spatial bins distribution distances
  max_spatial_bins <- max(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)
  true_spatial_bins_hist <- hist(true_active_spatial_bins_per_cell, 
                                 breaks=seq(0, max_spatial_bins, length.out=20), plot=F)
  true_spatial_bins_prob_dist <- true_spatial_bins_hist$density / sum(true_spatial_bins_hist$density)
  
  gen_spatial_bins_hist <- hist(generated_active_spatial_bins_per_cell, 
                                breaks=seq(0, max_spatial_bins, length.out=20), plot=F)
  gen_spatial_bins_prob_dist <- gen_spatial_bins_hist$density / sum(gen_spatial_bins_hist$density)
  
  spatial_bins_kld <- KLD(true_spatial_bins_prob_dist, gen_spatial_bins_prob_dist)
  
  
  # Compare peaks bins distribution distances
  max_peaks <- max(generated_peaks, true_peaks)
  true_peaks_hist <- hist(true_peaks, breaks=seq(0, max_peaks, length.out=10), plot=F)
  true_peaks_prob_dist <- true_peaks_hist$density / sum(true_peaks_hist$density)
  
  gen_peaks_hist <- hist(generated_peaks, breaks=seq(0, max_peaks, length.out=10), plot=F)
  gen_peaks_prob_dist <- gen_peaks_hist$density / sum(gen_peaks_hist$density)
  
  peaks_kld <- KLD(true_peaks_prob_dist, gen_peaks_prob_dist)
  
  
  cost <- sum(c(peaks_kld$sum.KLD.px.py,
                spatial_bins_kld$sum.KLD.px.py,
                si_kld$sum.KLD.px.py))
  
  return(cost)
}


JSD_distance <- function(generated_SI,
                         generated_time_bins,
                         generated_peaks,
                         generated_active_spatial_bins_per_cell,
                         true_SI,
                         true_time_bins_per_cells,
                         true_peaks,
                         true_active_spatial_bins_per_cell,
                         use_SI,
                         weight_vector=c(SI=1.4, 
                                         peaks=1, 
                                         sb=1)) {
  # Compare SI distribution distances
  true_SI <- true_SI[[1]]
  generated_SI <- generated_SI[[1]]
  max_SI <- max(c(generated_SI, true_SI))
  true_si_hist <- hist(true_SI, breaks=seq(0, max_SI, length.out=31), plot=F)
  true_si_prob_dist <- true_si_hist$density / sum(true_si_hist$density)
  
  gen_si_hist <- hist(generated_SI, breaks=seq(0, max_SI, length.out=31), plot=F)
  gen_si_prob_dist <- gen_si_hist$density / sum(gen_si_hist$density)
  
  si_JSD <- JSD(rbind(gen_si_prob_dist, true_si_prob_dist),)
  
  
  # # Compare time bins distribution distances
  # max_time_bins <- max(generated_time_bins, true_time_bins_per_cells)
  # true_time_bins_hist <- hist(true_time_bins_per_cells, breaks=seq(0, max_time_bins, length.out=100), plot=F)
  # true_time_bins_prob_dist <- true_time_bins_hist$density / sum(true_time_bins_hist$density)
  # 
  # gen_time_bins_hist <- hist(generated_time_bins, breaks=seq(0, max_time_bins, length.out=100), plot=F)
  # gen_time_bins_prob_dist <- gen_time_bins_hist$density / sum(gen_time_bins_hist$density)
  # 
  # time_bins_kld <- JSD(rbind(true_time_bins_prob_dist, gen_time_bins_prob_dist))
  
  # Compare active spatial bins distribution distances
  max_spatial_bins <- max(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)
  true_spatial_bins_hist <- hist(true_active_spatial_bins_per_cell, 
                                 breaks=seq(0, max_spatial_bins, length.out=15), plot=F)
  true_spatial_bins_prob_dist <- true_spatial_bins_hist$density / sum(true_spatial_bins_hist$density)
  
  gen_spatial_bins_hist <- hist(generated_active_spatial_bins_per_cell, 
                                breaks=seq(0, max_spatial_bins, length.out=15), plot=F)
  gen_spatial_bins_prob_dist <- gen_spatial_bins_hist$density / sum(gen_spatial_bins_hist$density)
  
  spatial_bins_JSD <- JSD(rbind(true_spatial_bins_prob_dist, gen_spatial_bins_prob_dist))
  
  
  # Compare peaks bins distribution distances
  max_peaks <- max(generated_peaks, true_peaks)
  true_peaks_hist <- hist(true_peaks, breaks=seq(0, max_peaks, length.out=8), plot=F)
  true_peaks_prob_dist <- true_peaks_hist$density / sum(true_peaks_hist$density)
  
  gen_peaks_hist <- hist(generated_peaks, breaks=seq(0, max_peaks, length.out=8), plot=F)
  gen_peaks_prob_dist <- gen_peaks_hist$density / sum(gen_peaks_hist$density)
  
  peaks_JSD <- JSD(rbind(true_peaks_prob_dist, gen_peaks_prob_dist))
  
  cost <- sum(c(weight_vector["peaks"] * peaks_JSD, 
                weight_vector["sb"] * spatial_bins_JSD, 
                weight_vector["SI"] * si_JSD))
  
  return(cost)
}



KS_distance  <- function(generated_SI,
                         generated_time_bins,
                         generated_peaks,
                         generated_active_spatial_bins_per_cell,
                         true_SI,
                         true_time_bins_per_cells,
                         true_peaks,
                         true_active_spatial_bins_per_cell,
                         use_SI,
                         weight_vector=c(SI=1, 
                                         peaks=1, 
                                         sb=1)) {
  
  #print(weight_vector[["SI"]])
  
  cost <- sum(c(weight_vector["SI"] * ks.test(generated_SI[[1]], true_SI[[1]])$statistic,
                weight_vector["peaks"] * ks.test(generated_peaks, true_peaks)$statistic,
                weight_vector["sb"] * ks.test(generated_active_spatial_bins_per_cell, true_active_spatial_bins_per_cell)$statistic))
  
  
  return(cost)
}

place_cell_params_sim_objective <- 
  function(X_train,
           true_spike_train,
           stim_trace,
           place_cell_pct,
           use_SI=F,
           remove_edges=F,
           dt=0.05,
           dist_fun="KS") {
    
    y <- c()
    for (row in 1:nrow(X_train)){
      
      cost_vec <- 
        spike_train_similiarity_cost(true_spike_train,
                                     stim_trace,
                                     dist_fun,
                                     X_train[row,"average_width"],
                                     X_train[row,"sd_width"],
                                     X_train[row,"noise"],
                                     X_train[row,"double_peak_pct"],
                                     place_cell_percentage = place_cell_pct,
                                     remove_edges=remove_edges,
                                     dt=dt,
                                     use_SI=use_SI)
      
      
      print(sum(cost_vec))
      y <- c(y, sum(cost_vec))
    }
    
    return(y)
  }

place_cell_params_sim_objective_new <- 
  function(X_train,
           true_spike_train,
           stim_trace,
           place_cell_pct,
           use_SI=F,
           remove_edges=F,
           dt=0.05,
           dist_fun="KS") {
    
    y <- c()
    for (row in 1:nrow(X_train)){
      
      cost_vec <- 
        spike_train_similiarity_cost_new(true_spike_train,
                                         stim_trace,
                                         dist_fun,
                                         X_train[row,"average_width"],
                                         X_train[row,"sd_width"],
                                         X_train[row,"noise"],
                                         X_train[row,"double_peak_pct"],
                                         remove_edges=remove_edges,
                                         0,
                                         place_cell_percentage = place_cell_pct,
                                         dt=dt,
                                         use_SI=use_SI)
      
      
      print(sum(cost_vec))
      y <- c(y, sum(cost_vec))
    }
    
    return(y)
  }


contour_maps <- function(predicted, X_star, folder="heatmaps_new/") {
  dir.create(folder)
  values <- apply(X_star, 2, unique)
  
  ind_mat <-  combn(1:len(values), 2)  
  
  for (idx in 1:ncol(ind_mat)) {
    ind <- ind_mat[,idx]
    
    v1 <- unlist(values[ind[1]])
    v2 <- unlist(values[ind[2]])
    
    mean_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                            nrow=len(v1))
    
    sd_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                          nrow=len(v1))
    
    min_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                           nrow=len(v1))
    for (i in 1:len(v1))  {
      for (j in 1:len(v2)) {
        
        mean_cont_mat[i,j] <-  
          mean(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
        
        sd_cont_mat[i,j] <-  
          sd(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
        
        min_cont_mat[i,j] <- 
          min(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
      }
    }
    
    color = colorRampPalette(rev(brewer.pal(n = 7,  name = "Spectral")))
    
    rownames(mean_cont_mat) <- round(v1, digits=3)
    colnames(mean_cont_mat) <- round(v2, digits=3)
    rownames(sd_cont_mat) <- round(v1, digits=3)
    colnames(sd_cont_mat) <- round(v2, digits=3)
    rownames(min_cont_mat) <- round(v1, digits=3)
    colnames(min_cont_mat) <- round(v2, digits=3)
    
    print("PLOTTING INNNNNNNNNNNN")
    print(folder)
    png(filename=sprintf("%s/MEAN_%s_%s_heatmap.png",
                         folder,
                         colnames(X_star)[ind[1]],
                         colnames(X_star)[ind[2]]),
        units="in",
        width=3.5,
        height=3.5,
        res=250)
    pheatmap(mean_cont_mat, cluster_rows=F, cluster_cols=F, col=color(200))
    dev.off()
    
    png(filename=sprintf("%s/SD_%s_%s_heatmap.png",
                         folder,
                         colnames(X_star)[ind[1]],
                         colnames(X_star)[ind[2]]),
        units="in",
        width=3.5,
        height=3.5,
        res=250)
    pheatmap(sd_cont_mat, cluster_rows=F, cluster_cols=F, col=color(200))
    dev.off()
    
    png(filename=sprintf("%s/MIN_%s_%s_heatmap.png",
                         folder,
                         colnames(X_star)[ind[1]],
                         colnames(X_star)[ind[2]]),
        units="in",
        width=3.5,
        height=3.5,
        res=250)
    pheatmap(min_cont_mat, cluster_rows=F, cluster_cols=F, col=color(200))
    dev.off()
    
  }
}

generated_simulated_activity_non_stationary <- 
  function(percentage,
           jump_lambda=0.1,
           jump_size=1,
           num_of_jumps=1,
           plot=F, 
           average_width=0.14,
           sd_width=0.05,
           noise=0.015, 
           double_peak_pct=0, 
           triple_peak_pct=0,
           average_width_2=0.07,
           n_bins = 24,
           first_bin=1,
           stim_trace,
           spike_train,
           return_non_place=F,
           fixed_pois=-1) {
    
    n = nrow(spike_train)
    fixed_fr <- rowMeans(spike_train) / dt
    n_neurons = n
    mu_rates=0.08
    std_rates=0.02
    
    noise_level=noise
    average_sigma=average_width
    std_sigma_f=sd_width
    average_sigma2=average_width_2
    
    
    pref_positions = seq(0.5, n_bins + 0.5, length.out=n_neurons)
    pref_positions_2 = seq(0.5, n_bins + 0.5, length.out=n_neurons)
    pref_positions_2 <- pref_positions_2[sample(1:n_neurons, n_neurons)]
    
    
    average_rates=fixed_fr[sample(1:length(fixed_fr), len(fixed_fr))]
    
    mu_sigma=log((average_sigma ** 2) / sqrt(std_sigma_f + average_sigma ** 2))
    std_sigma=sqrt(log(std_sigma_f / (average_sigma ** 2) + 1));
    normed_sigma=rlnorm(n_neurons, meanlog=mu_sigma, sdlog=std_sigma);
    sigma=normed_sigma*n_bins
    
    mu_sigma2=log((average_sigma2 ** 2) / sqrt(std_sigma_f + average_sigma2 ** 2))
    std_sigma2=sqrt(log(std_sigma_f / (average_sigma2 ** 2) + 1));
    normed_sigma2=rlnorm(n_neurons, meanlog=mu_sigma2, sdlog=std_sigma2);
    sigma2=normed_sigma*n_bins
    
    
    switch_positions <- matrix(c(pref_positions), nrow=1)
    current_poisition <- pref_positions
    
    if (num_of_jumps > 0 ) {
      for (i in 1:num_of_jumps) {
        
        
        switch_pos <- unlist(lapply(current_poisition, 
                                    function(pos) {
                                      delta <- rnorm(1, mean=jump_size, sd=3); new_pos <- pos + delta; 
                                      
                                      if (new_pos < 0.5) {
                                        new_pos <- 0.5 + abs(0.5 - new_pos)
                                      } else if (new_pos > (n_bins + 0.5)) {
                                        new_pos <- (n_bins + 0.5) - (new_pos - (n_bins + 0.5))
                                      } 
                                      return(new_pos)
                                    }))
        
        switch_positions <- rbind(switch_positions,switch_pos)
        current_position <- switch_pos
        
      }
      
    }
    
    get_tuning <- function(n,i) {return((1/sqrt(2*pi*sigma[n])) * exp((-(1:n_bins - switch_positions[i,n]) ** 2) / (2 * sigma[n] ** 2)))}
    get_tuning2 <- function(n,i) {return((1/sqrt(2*pi*sigma[n])) * (exp((-(1:n_bins - switch_positions[i,n]) ** 2) / (2 * sigma[n] ** 2)) +
                                                                      runif(1,0.1,1) * exp((-(1:n_bins - pref_positions_2[n]) ** 2) / (2 * sigma2[n] ** 2))))}
    
    double_peak_ind <- -1
    
    if (double_peak_pct > 0) {
      double_peaked <- floor(rnorm(1, mean=n_neurons * double_peak_pct, sd=3))
      
      while(double_peaked < 1) {
        print("resampling")
        double_peaked <- floor(rnorm(1, mean=n_neurons * double_peak_pct, sd=3))
      }
      
      double_peak_ind <- sample(1:n_neurons, double_peaked)
    }
    
    
    non_place = -1
    
    if (percentage != 1) {
      n_non_place_neurons <- n_neurons - floor(rnorm(1 ,mean=n_neurons * percentage, sd=10))
      
      while(n_non_place_neurons < 0 || n_non_place_neurons > n_neurons) {
        n_non_place_neurons <- n_neurons - floor(rnorm(1 ,mean=n_neurons * percentage, sd=10))
      }
      
      non_place = sample(1:n_neurons, n_non_place_neurons)
    }
    
    mixed_rates <- fixed_fr[sample(1:length(fixed_fr), len(fixed_fr))]
    
    switch_tuning_list <- list()
    
    for (idx in 1:nrow(switch_positions)) { 
      
      tuning_curves <- 
        lapply(1:n_neurons, 
               function(n) {
                 
                 if ((!all(double_peak_ind == -1)) && (n %in% double_peak_ind)) {
                   
                   tuning <- get_tuning2(n, idx)
                 } else {
                   
                   tuning <- get_tuning(n, idx)  
                 }
                 
                 factor <- average_rates[n]/(mean(tuning) + 10^-30)
                 tuning <- tuning * factor + abs(rnorm(n_bins) * noise_level)
                 return(tuning)
               })
      
      
      if (percentage != 1) {
        
        # Randomly draw n_non neurons and set them with random tuning curve
        for (i in non_place) {
          tuning_curves[[i]] <- rep(1, times=n_bins) + abs(rnorm(n_bins) * noise_level)
        }
      }
      
      
      # assert(all(unlist(lapply(sample(1:400, 50), function(i) { all(tuning_curves[[i]] == mt[i,])}))))
      
      tc2 <- lapply(tuning_curves, function(r) {r / sum(r)})
      tc3 <- lapply(1:length(tuning_curves), function(i) {tc2[[i]] * mixed_rates[i]})
      mt <- t(matrix(unlist(tc3), nrow = n_bins))
      
      switch_tuning_list <- append(switch_tuning_list, list(mt))
    }
    
    n_id <- sample(1:n, 1)
    
    barplot(switch_tuning_list[[1]][n_id,], col=adjustcolor("red", alpha=0.2))
    
    if (num_of_jumps > 0) {
      for (t_id in 2:len(switch_tuning_list)) {
        
        barplot(switch_tuning_list[[t_id]][n_id,], col=adjustcolor("red", alpha=0.2), add=T)  
        
      }
    }
    
    if (percentage != 1)  {
      place_neurons <- which(!1:n %in% non_place)
    } else {
      place_neurons <- 1:n
    }
    
    switch_vector <- rep(0, times=n)
    switch_vector[place_neurons] <- rpois(len(place_neurons),
                                          lambda=jump_lambda)
    
    if (percentage != 1) {
      assert(all(switch_vector[non_place] == 0))
    }
    
    processed <- preprocess_spike_train(spike_train, stim_trace, verbose=F)
    
    if(fixed_pois != -1) {
      pois_factor <- fixed_pois
    } else {
      pois_factor <- 
        currate_spike_train_non_stationary(switch_tuning_list,
                                           switch_vector,
                                           processed$working_time_bins_per_cells,
                                           stim_trace,
                                           verbose=T)
    }
    
    gst <- generate_spike_trains_non_stationary(switch_tuning_list = switch_tuning_list,
                                                switch_vector = switch_vector,
                                                stim_trace,
                                                factor=pois_factor)
    
    
    plot_neur_raster(stim_trace, gst,  which(switch_vector != 0)[floor(runif(1,1,sum(switch_vector != 0)))])
    
    return(list(tuning=switch_tuning_list,
                place_neurons=place_neurons,
                generated_spike_train=gst,
                switch_vector=switch_vector,
                pois_factor=pois_factor))
  }

currate_spike_train_non_stationary <- function(switch_tuning_list, 
                                               switch_vector,
                                               true_tb, 
                                               stim_trace, 
                                               by_pval_only=F,
                                               verbose=F) {
  
  if (dt == 0.05) {
    fact <- seq(1.2, 3.6, by=0.05)
  } else {
    fact <- seq(0.4, 3.6, by=0.05)
    
  }
  
  
  
  pvals <- c()
  lengths <- c()
  for (f in fact) {
    gst <- generate_spike_trains_non_stationary(switch_tuning_list=switch_tuning_list,
                                                switch_vector=switch_vector,
                                                stim_trace=stim_trace, 
                                                factor=f,
                                                verbose=F)
    
    sim <- preprocess_spike_train(gst, stim_trace, verbose=F)
    pvals <- c(pvals, ks.test(sim$working_time_bins_per_cells, true_tb)$p.value)
    lengths <- c(lengths, (len(sim$working_time_bins_per_cells) / len(true_tb)))
    if (verbose) {
      print(f)
      print(ks.test(sim$working_time_bins_per_cells, true_tb)$p.value)
      print(len(sim$working_time_bins_per_cells) / len(true_tb))
    }
  }
  
  o1 <- order(pvals, decreasing = T)
  o2 <- order((lengths - 1) ** 2)
  
  rank <- rep(0, times=len(fact))
  names(rank) <- fact
  
  rank[o1] <- rank[o1] + 1:len(fact)
  rank[o2] <- rank[o2] + 1:len(fact)
  
  if (by_pval_only) {
    return(fact[which.max(pvals)])
  }
  
  return(fact[which.min(rank)])
}

generate_spike_trains_non_stationary <- function(switch_tuning_list, 
                                                 switch_vector,
                                                 stim_trace, 
                                                 factor=1,
                                                 min_session_dur_jump=0.15,
                                                 max_session_dur_jump=0.85,
                                                 jump_step=0.05,
                                                 verbose=F) {
  
  
  st <- t(apply(switch_tuning_list[[1]], 
                1, 
                function(tc) {rpois(len(stim_trace), tc[stim_trace] / (1 * factor))}))
  
  
  start_frame <- floor(min_session_dur_jump * len(stim_trace))
  end_frame <- floor(max_session_dur_jump * len(stim_trace))
  
  switch_neurons <- which(switch_vector != 0)
  
  
  switch_frames <- 
    lapply(switch_vector[switch_neurons], 
           function(num_switches) {
             jump_slots <- seq(min_session_dur_jump, max_session_dur_jump, by=jump_step)
             if (num_switches > len(jump_slots)) {
               num_switches <- len(jump_slots) 
             }
             sort(floor(sample(jump_slots, num_switches) * len(stim_trace)))
           })
  
  names(switch_frames) <- switch_neurons
  
  
  for (switch_neur in switch_neurons) {
    
    switch_frames_list <- switch_frames[[as.character(switch_neur)]]
    
    
    curr_frame <- 1
    switch_index = 1
    curr_tuning <- switch_tuning_list[[switch_index]]
    
    
    #plot_neur_raster(stim_trace, st, switch_neur)
    #barplot(switch_tuning_curves[switch_neur,], col=adjustcolor("red", alpha=0.3))
    #barplot(tuning_curves[switch_neur,], col=adjustcolor("blue", alpha=0.3), add=T)
    switch_st <- c()
    
    #plot_neur_raster(stim_trace, st, switch_neur)
    for (switch_frame in switch_frames_list) {
      #print(switch_frames_list)
      
      ind <- curr_frame:(switch_frame - 1)
      switch_st <- c(switch_st, rpois(len(ind), curr_tuning[switch_neur,][stim_trace[ind]] / (1 * factor)))
      
      # Switch tuning_curve
      if (verbose) { 
        print(sprintf("Switching from TC %d", switch_index))
      }
      
      if (switch_index == 1) {
        # First tuning curve, jump to next
        switch_index = switch_index + 1
      } else if(switch_index == len(switch_tuning_list)) {
        # Last tuning curve, jump to previos
        switch_index = switch_index - 1
      } else {
        # Randomly jump to next / previous tuning curve
        switch_index = switch_index + sample(c(-1,1), 1)
      }
      
      if (verbose) { 
        print(sprintf("Switching to TC %d", switch_index))
      }
      
      
      curr_tuning <- switch_tuning_list[[switch_index]]
      curr_frame <- switch_frame
    }
    
    ind <- curr_frame:ncol(st)
    switch_st <- c(switch_st, rpois(len(ind),  curr_tuning[switch_neur,][stim_trace[ind]] / (1 * factor)))
    # 
    # print(len(switch_st))
    # print(ncol(st))
    # 
    #  plot(stim_trace, type="l")
    #  points(y=stim_trace[which(switch_st != 0)], x=which(switch_st != 0), pch=19, col="red")
    # 
    st[switch_neur,] <- switch_st
    
  }
  
  return(st)
}