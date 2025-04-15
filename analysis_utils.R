dt=0.05
shuffle_type='cyclic'

active_bins_threshold=5
firing_rate_threshold=0
num_shuffles=1000;
num_of_subsample_repetitions=500
tuning_significance_threshold=0.05
bin_size_active_time_bins = 10
plot_res = T
verbose = T
save_res = T
spatial_bins_t = 1:24
spatial_bins_t_no_edges <- 1:20

ph <- function(mt, ...) {pheatmap(mt, cluster_cols = F, cluster_rows=F, ...)}

spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
rdylbu_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "RdYlBu")))

base_plot_theme <- theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black"),
                         legend.position="NA",
                         panel.border = element_blank(),
                         panel.background = element_blank())


##### ONE LINER FUNCTIONS 
bin  <- function(v, bin_size) {return(cut(v, seq(0,max(v) + bin_size, by = bin_size)))}
gen_num_frames <- function(stim_trace) {return(round(rnorm(1, mean = mean(table(stim_trace)[3:22]),  sd = sd(table(stim_trace)[3:22]) / len(3:22))))}
len <- length

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,1), 
         symbols = c("****", "***", "**", "*", "n.s."))
}

sem <- function(r) {
  ln <- sum(!is.na(r))
  if (ln <= 1) {
    return(0)
  }
  return(sqrt(var(as.numeric(r), na.rm=T))/sqrt(ln))
}

compute_tuning <- function(spike_train, stim_trace, remove_edges=F) {
  #spatial_bins <- unique(stim_trace)
  spatial_bins <- unlist(ifelse(remove_edges,
                                list(spatial_bins_t_no_edges), 
                                list(spatial_bins_t)))
  num_of_frames <- length(stim_trace)
  num_of_cells <- nrow(spike_train)
  bin_prob <- table(stim_trace) / length(stim_trace)
  
  
  if (len(bin_prob) < len(spatial_bins)) {
    missing <- spatial_bins[which(!spatial_bins %in% as.numeric(names(bin_prob)))]
    names_all <- c(as.numeric(names(bin_prob)), missing)
    bin_prob <- c(bin_prob, rep(0, times=length(missing)))
    names(bin_prob) <- names_all
  }
  
  bin_prob <- bin_prob[spatial_bins]
  
  #print(spike_train)
  tuning_curve <- 
    lapply(spatial_bins, 
           function(sb) {
             
             bin_ind <- which(stim_trace == sb)
             
             if (len(bin_ind) == 0) {
               return(0)
             }
             
             if (nrow(spike_train) == 1){
               return(mean(spike_train[,bin_ind])/dt)
             }
             
             if (len(bin_ind) == 1) {
               spike_train <- as.matrix(spike_train[,bin_ind])
               return(rowMeans(spike_train) / dt)
             } 
             
             if (len(bin_ind) == 0) {
               return(rep(0, times=nrow(spike_train)))
             }
             
             return(rowMeans(spike_train[,bin_ind])/dt)
           })
  
  tuning_curve <- do.call(cbind,tuning_curve)
  colnames(tuning_curve) <- spatial_bins
  
  return(list(bin_prob, tuning_curve))
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


compute_SI <- function(stim_prob, tuning_curve, firing_rate) {
  
  eps <- 10^-30
  
  # For each cell, this is r_i / r_avg
  # The result matrix, has N cells x S spatial bins
  normalized_tuning <- apply(tuning_curve, 
                             2, 
                             function(c) {
                               c/(firing_rate + eps)
                             })
  
  if (nrow(tuning_curve) == 1) {
    normalized_tuning <- t(as.matrix(normalized_tuning))
  }
  
  # for each cell, compute log2(r_i / r_avg)
  log_normalized_tuning <- log(normalized_tuning + eps, 2)
  
  
  # for reach cell, compute r_i / r_avg * log2(r_i / r_avg)
  tmp <- lapply(1:ncol(normalized_tuning),
                function(i) {
                  normalized_tuning[,i] * 
                    log_normalized_tuning[,i]
                })
  
  
  tmp <- do.call(cbind,tmp)
  
  # for each cell, compute p_i * (r_i / r_avg) * log2(r_i / r_avg)
  final <- t(apply(tmp, 1, function(r) {r * stim_prob}))
  information <- rowSums(final)
  
  if (min(information) < 0) {
    print("ERROR! Sanity check violated!!")
    print(information)
    print(stim_prob)
    print(tuning_curve)
    print(firing_rate)
    assert(0==1)
  }
  
  return(list(information, information * firing_rate))
}

compute_pois_product <- function(spike_train, stim_trace) {
  
  rm <- rowMeans(spike_train)
  rounded <- round(spike_train)
  for (sb in sort(unique(stim_trace))) {
    
    ind <- which(stim_trace == sb)
    rs <- rounded[,ind]
    res <- cbind(res, 
                 unlist(lapply(1:nrow(spike_train), function(idx) {sum(ppois(rs[idx,], lambda=mean(rs[idx,])))})))
  }
  
  return(apply(res, 1, sum))
}

compute_MI <- function(stim_trace, spike_train) {
  eps <- 10^-30
  spatial_bins <- spatial_bins_t
  rounded_spike_train <- ceiling(spike_train)
  
  MI <- lapply(1:nrow(rounded_spike_train),
               
               function (n_indx) {
                 neur <- rounded_spike_train[n_indx,]
                 spikes <- c(0,1)
                 #spikes <- sort(unique(neur))
                 
                 joint_prob_mat <- matrix(rep(0, times=len(spikes) * len(spatial_bins)), 
                                          nrow=len(spatial_bins))
                 rownames(joint_prob_mat) <- spatial_bins
                 colnames(joint_prob_mat) <- spikes
                 
                 for (sbin in spatial_bins) {
                   ind <- which(stim_trace == sbin)
                   # for (s_count in spikes) {
                   #   s_count_idx <- which(s_count == spikes)
                   #   joint_prob_mat[sbin, s_count_idx] <- (sum(neur[ind] == s_count))
                   # }
                   
                   joint_prob_mat[sbin, 1] <- (sum(neur[ind] == 0))
                   joint_prob_mat[sbin, 2] <-  len(neur[ind]) - joint_prob_mat[sbin, 1]
                 }
                 
                 sm <- sum(joint_prob_mat)
                 spike_prob <- colSums(joint_prob_mat) / sm
                 bin_prob <- rowSums(joint_prob_mat) / sm 
                 
                 MI <- c()
                 
                 
                 for (sbin in spatial_bins) {
                   for (s_count in spikes) {
                     s_count_idx <- which(s_count == spikes)
                     pso <- joint_prob_mat[sbin, s_count_idx] / sm
                     MI <- c(MI,
                             pso * log((pso/(spike_prob[s_count_idx] * bin_prob[sbin])) + eps))
                     
                   }
                 }
                 
                 return(sum(MI))
               })
  
  return(unlist(MI))
  
}

compute_place_signif_pois_product <- function(i,
                                              working_cells_spike_train,
                                              stim_trace,
                                              working_PP,
                                              shuffle_type="r",
                                              num_shuffles=1000,
                                              return_pval=T,
                                              verbose=T) {0
  
  shuffle <- generate_shuffled_spike_train(spike_train = working_cells_spike_train,
                                           idx = i,
                                           num_of_shuffles = num_shuffles,
                                           shuffle_type = shuffle_type,
                                           verbose=verbose)
  
  
  # spike train in the size of Nshuffles x T 
  shuffled_spike_train <- t(array_reshape(shuffle, c(ncol(working_cells_spike_train),num_shuffles)))
  rm(shuffle)
  

  
  
  shuffled_tuning <- compute_tuning(shuffled_spike_train, stim_trace)
  
  shuffled_PP <- compute_pois_product(shuffled_spike_train, stim_trace)
  
  if (!return_pval) {
    return(shuffled_PP)
  }
  
  f <- ecdf(shuffled_PP)
  pval <- 1 - f(working_PP[i])
  
  return(pval)
}


compute_place_signif_MI <- function(i,
                                    working_cells_spike_train,
                                    stim_trace,
                                    working_MI,
                                    shuffle_type="r",
                                    num_shuffles=1000,
                                    return_pval=T,
                                    verbose=T) {
  
  #print(working_cells_spike_train[i,1:20])
  print(working_MI[i])
  print(rowSums(working_cells_spike_train)[i])
  
  shuffle <- generate_shuffled_spike_train(spike_train = working_cells_spike_train,
                                           idx = i,
                                           num_of_shuffles = num_shuffles,
                                           shuffle_type = shuffle_type,
                                           verbose=F)
  
  
  # spike train in the size of Nshuffles x T 
  shuffled_spike_train <- t(array_reshape(shuffle, c(ncol(working_cells_spike_train),num_shuffles)))
  
  print(shuffled_spike_train[5,1:20])
  print(shuffled_spike_train[8,1:20])
  rm(shuffle)
  
  
  shuffled_MI <- compute_MI(stim_trace, shuffled_spike_train)
  
  
  if (!return_pval) {
    return(shuffled_MI)
  }
  
  #print(shuffled_MI)
  f <- ecdf(shuffled_MI)
  pval <- 1 - f(working_MI[i])
  
  return(pval)
}



compute_place_signif <- function(i,
                                 working_cells_spike_train,
                                 stim_trace,
                                 working_firing_rate,
                                 working_SI,
                                 shuffle_type="r",
                                 num_shuffles=1000,
                                 return_pval=T,
                                 rate_to_use=1,
                                 verbose=F) {
  shuffle <- generate_shuffled_spike_train(spike_train = working_cells_spike_train,
                                           idx = i,
                                           num_of_shuffles = num_shuffles,
                                           shuffle_type = shuffle_type,
                                           verbose=verbose)
  
  
  # spike train in the size of Nshuffles x T 
  shuffled_spike_train <- t(array_reshape(shuffle, c(ncol(working_cells_spike_train),num_shuffles)))
  rm(shuffle)
  
  
  # Sanity check
  shuffled_rate <- rowMeans(shuffled_spike_train) / dt
  assert(all(shuffled_rate == shuffled_rate[1]))
  assert(shuffled_rate[1] == working_firing_rate[i])
  
  tmp <- compute_tuning(shuffled_spike_train, stim_trace)
  shuffled_stim_prob <- tmp[[1]]
  shuffled_tuning <- tmp[[2]]
  
  # Sanity check
  # assert(all(shuffled_stim_prob == stim_prob))
  # assert(rowSums())
  
  # Compute SI for all shuffles
  shuffled_SI <- compute_SI(shuffled_stim_prob,
                            shuffled_tuning,
                            shuffled_rate)
  
  if (!return_pval) {
    return(shuffled_SI)
  }
  
  
  # Compute significance, by estimating cdf
  f <- ecdf(shuffled_SI[[rate_to_use]])
  pval <- 1 - f(working_SI[[rate_to_use]][i])
  
  
  rm(shuffled_SI)
  rm(shuffled_tuning)
  rm(shuffled_stim_prob)
  rm(shuffled_spike_train)
  return(pval)
}

generate_shuffled_spike_train <- 
  function(spike_train, 
           idx,
           num_of_shuffles=1000,
           shuffle_type="r",
           granularity=1,
           verbose=F) {
    n_frames <- ncol(spike_train)
    #n_cells <- nrow(spike_train)
    0
    working_spike_train <- t(as.matrix(spike_train[idx:(idx + granularity - 1),]))
    if (granularity > 1) working_spike_train <- t(working_spike_train)
    
    shuffle <- array(rep(0, times=n_frames * granularity * num_of_shuffles),
                     dim=c(granularity, n_frames, num_of_shuffles))
    
    if (verbose) {
      if(shuffle_type == "c") {
        print(sprintf("%d. ##### Running %d cyclic shuffles", idx, num_of_shuffles))
      } else {
        print(sprintf("%d. ##### Running %d random shuffles", idx, num_of_shuffles))
      }
    }
    if (shuffle_type == "c") {
      # I equals to shuffle idx
      for (i in 1:num_of_shuffles) {
        for(cell_idx in 1:granularity) {
          # random shift sized n whereas n << m.
          shift <- round(runif(n=1, min=1, max=(n_frames - 1)))
          
          # copy first 1 ... m - n to n + 1 ..... m
          shuffle[cell_idx,(shift + 1):n_frames,i] <-  
            working_spike_train[cell_idx,1:(n_frames-shift)]
          
          # then copy m - n + 1 ..... m to 1 ... n
          shuffle[cell_idx,1:shift,i] <-  
            working_spike_train[cell_idx,(n_frames-shift + 1):n_frames]
        }
      }
    } else {
      for (i in 1:num_of_shuffles) {
        for(cell_idx in 1:granularity) {
          # Randomly sample indices 
          shuffle[cell_idx,1:n_frames,i] <-  
            working_spike_train[cell_idx,sample(1:n_frames, n_frames)]
        }
      }
    }
      return(shuffle)
  }

preprocess_spike_train <- function(spike_train, stim_trace, verbose=T, active_bins_threshold=5) {
  time_bins_activity_per_cell <- 
    sapply(1:nrow(spike_train), function(i) {sum(as.numeric(spike_train[i,] > 0))})
  
  activity_per_cell <- rowMeans(spike_train) / dt
  selected_cells_ind <- activity_per_cell > firing_rate_threshold & time_bins_activity_per_cell > active_bins_threshold
  
  if (verbose) {
    print(sprintf("#### There are %d sufficiently active cells", sum(selected_cells_ind)))
  }
  
  selected_cells_ind <- which(selected_cells_ind)
  
  # Filter out sufficiently active cells 
  working_cells_spike_train <- spike_train[selected_cells_ind,]
  working_firing_rate <- activity_per_cell[selected_cells_ind]
  working_time_bins_per_cells <- time_bins_activity_per_cell[selected_cells_ind]
  
  return(list(working_cells_spike_train=working_cells_spike_train,
              working_firing_rate=working_firing_rate,
              working_time_bins_per_cells=working_time_bins_per_cells,
              ind=selected_cells_ind))
}

get_spike_train_and_stim_trace_from_path <- function(path, session, old=F, equalize_prior=T, simulate=F, remove_edges=F) {
  
  if (old) {
    m <- readMat(paste(path, "/stimulus_trace.mat", sep=""))
    stim_trace <- m$stimulus.trace
    
    
    m <- readMat(paste(path, spike_train_filename, sep=""))
    spike_train <- m$spike.train
  } else {  
    m <- readMat(paste(path, stimulus_trace_filename, sep=""))
    
    if(len(grep("Left", path)) > 0) {
      stim_trace <- m$position.left.trials
    } else {
      stim_trace <- m$position.right.trials
    }
    
    
    m <- readMat(paste(path, spike_train_filename, sep=""))
    
    if(len(grep("Left", path)) > 0) {
      spike_train <- m$spikes.left.trials
    } else {
      spike_train <- m$spikes.right.trials
    }
  }
  
  tmp_path <- path
  tmp_spike_train <- spike_train
  tmp_stimulus_trace <- stim_trace
  
  
  spike_train <- tmp_spike_train
  stim_trace <- tmp_stimulus_trace
  spike_train <- spike_train[[session]][[1]]
  spike_train <- t(as.matrix(spike_train))
  stim_trace <- stim_trace[[session]][[1]]
  stim_trace <-  as.vector(stim_trace)
  
  
  
  if (!old && equalize_prior) {
    equalize_bins <- c(1:3,22:24)
    indices <- 1:len(stim_trace)
    
    
    final_ind <- rep(T, times=len(stim_trace))
    
    for (ebin in equalize_bins) {
      print(ebin)
      
      num_to_draw <- gen_num_frames(stim_trace)
      if (num_to_draw > len(which(stim_trace == ebin))) {
        
        next
      }
      
      sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
      tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin       
      final_ind <- final_ind & tmp
    }
    
    spike_train <- spike_train[,final_ind]
    stim_trace <- stim_trace[final_ind]
  }
  
  if (!old && remove_edges) {
    non_edges_ind <- which(stim_trace %in% c(3:22))
    stim_trace <- stim_trace[non_edges_ind]
    spike_train <- spike_train[,non_edges_ind]
    stim_trace <- stim_trace - 2 # Refer to bin 3 as bin 1
  }
  
  if (simulate) {
    fr <- rowMeans(spike_train) / dt
    processed_real <- preprocess_spike_train(spike_train, stim_trace)
    
    true_cells_spike_train <- processed_real$working_cells_spike_train
    true_firing_rate <- processed_real$working_firing_rate
    true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
    
    print(sprintf("%s\\equalized\\session_%d", path, 2))
    params <- get_fit_params_figures_from_path(sprintf("%s\\equalized\\session_%d", path, session) , T)
    simulated_tuning_curve <-
      generate_tuning_curves_cost(n = nrow(spike_train),
                                  percentage = params$pct,
                                  average_width = params$params["average_width"], 
                                  sd_width = params$params["sd_width"],
                                  fixed_fr=fr,
                                  noise=params$params["noise"],
                                  double_peak_pct = params$params["double_peak_pct"],
                                  n_bins=ifelse(remove_edges, 20, 24),
                                  plot=F)
    
    pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                            true_time_bins_per_cells,
                                            stim_trace,
                                            verbose = T)
    
    generated_spike_train <- 
      generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                 stim_trace = stim_trace,
                                 factor=pois_factor,
                                 fs=1)
    
    orig_st <- spike_train
    orig_path <- path
    spike_train <- generated_spike_train
    return (list(spike_train, stim_trace, simulated_tuning_curve, orig_st))
  }
  
  
  return (list(spike_train, stim_trace))
}


get_tuning_properties <- function(spike_train, stim_trace, peaks_significance_threshold=0.3) {
  processed_real <- preprocess_spike_train_no_threshold(spike_train,
                                                        stim_trace=stim_trace)
  #processed_real <- preprocess_spike_train(spike_train, stim_trace=stim_trace)
  
  working_cells_spike_train <- processed_real$working_cells_spike_train
  firing_rate <- processed_real$working_firing_rate
  time_bins <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(working_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  tuning_curves <- tmp[[2]]
  rm(tmp)
  
  
  peaks <- unlist(apply(tuning_curves, 1, 
                             function(n) {
                               return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))
                               }))
  spatial_bins <- apply(tuning_curves, 1, function(n){sum(n>0)})
  #print(true_active_spatial_bins_per_cell)
  SI <- compute_SI(stim_prob, tuning_curves, firing_rate)
  
  tuning_properties <- data.frame(peaks=peaks,
                                  time_bins=time_bins,
                                  firing_rate=firing_rate,
                                  SI=SI[[1]],
                                  SI_sec=SI[[2]],
                                  spatial_bins=spatial_bins)
  
  return(list(tuning_properties=tuning_properties,
              tuning_curves=tuning_curves))
}

get_cost_compare_vectors_by_path <- function(path, true_spike_train, stim_trace, pval_func=wilcox.test, ret=F, smooth_lambda=-1, best=T, new_simulation=F, use_JSD=F, edge_free=F) {
  
  firing_rate <-rowMeans(true_spike_train) / dt
  processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(true_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  true_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  
  if (smooth_lambda == - 1){
    tuning_for_peaks <- true_tuning_curve
  } else {
    tuning_for_peaks <- 
      t(apply(true_tuning_curve, 1, function(trc) {
        barp <- barplot(trc, plot=F)
        smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
        smt$y[smt$y < 0] <- 0
        return(smt$y);
      }))
  }
  
  true_peaks <- unlist(apply(tuning_for_peaks, 1, 
                             function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
  true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
  #print(true_active_spatial_bins_per_cell)
  true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
  
  params <- get_fit_params_figures_from_path(path, best = best, use_JSD = use_JSD, edge_free = edge_free)
  SI_f <- c()
  tb_f <- c()
  peaks_f <- c()
  sb_f <- c()
  simulated_tuning_curve <-
    generate_tuning_curves_cost(n = nrow(true_spike_train),
                                percentage = params$pct,
                                average_width = params$params["average_width"], 
                                sd_width = params$params["sd_width"],
                                fixed_fr=firing_rate,
                                noise=params$params["noise"],
                                double_peak_pct = params$params["double_peak_pct"],
                                n_bins = ifelse(edge_free, 24, 20),
                                plot=F)
  pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                          true_time_bins_per_cells,
                                          stim_trace,
                                          verbose = F)
  
  
  for (i in 1:10) { 
    simulated_tuning_curve <-
      generate_tuning_curves_cost(n = nrow(true_spike_train),
                                  percentage = params$pct,
                                  average_width = params$params["average_width"], 
                                  sd_width = params$params["sd_width"],
                                  fixed_fr=firing_rate,
                                  noise=params$params["noise"],
                                  double_peak_pct = params$params["double_peak_pct"],
                                  plot=F)
    
    
    generated_spike_train <- 
      generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                 stim_trace = stim_trace,
                                 factor=pois_factor,
                                 fs=1)
    
    
    processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
    generated_cells_active_st <- processed_generated$working_cells_spike_train
    generated_firing_rate <- processed_generated$working_firing_rate
    generated_time_bins <- processed_generated$working_time_bins_per_cells
    
    tmp <- compute_tuning(generated_cells_active_st, stim_trace)
    gen_stim_prob <- tmp[[1]]
    gen_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    if (smooth_lambda == - 1){
      tuning_for_peaks <- gen_tuning_curve
    } else {
      tuning_for_peaks <- 
        t(apply(gen_tuning_curve, 1, function(trc) {
          barp <- barplot(trc, plot=F)
          smt <- smooth.spline(barp, trc, all.knots = T, lambda=smooth_lambda)
          smt$y[smt$y < 0] <- 0
          return(smt$y);
        }))
    }
    
    generated_peaks <- unlist(apply(tuning_for_peaks, 1, 
                                    function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
    generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
    generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
    
    
    if (ret) {
      return(list(true_SI=true_SI[[1]],
                  true_tb=true_time_bins_per_cells,
                  true_peaks=true_peaks,
                  true_sb=true_active_spatial_bins_per_cell,
                  gen_SI=generated_SI[[1]],
                  gen_tb=generated_time_bins,
                  gen_peaks=generated_peaks,
                  gen_sb=generated_active_spatial_bins_per_cell))
    }
    
    SI=pval_func(true_SI[[1]], generated_SI[[1]])$p.value
    tb=pval_func(true_time_bins_per_cells, generated_time_bins)$p.value
    peaks=pval_func(true_peaks, generated_peaks)$p.value
    sb=pval_func(true_active_spatial_bins_per_cell, generated_active_spatial_bins_per_cell)$p.value
    
    print(sprintf("SI %.3f, tb %.3f, peaks %.3f, sb %.3f",
                  SI,
                  tb,
                  peaks,
                  sb))
    SI_f <- c(SI_f, SI)
    tb_f <- c(tb_f, tb)
    peaks_f <- c(peaks_f, peaks)
    sb_f <- c(sb_f, sb)
    
  }
  
  return(c(mean(SI_f), mean(tb_f), mean(peaks_f), mean(sb_f)))
  
}

get_fit_params <- function(session_path, 
                           likelihood_path="JSD_simulations_likelihood_c", 
                           estimation_path="JSD_simulations_estimate",
                           just_likelihood=F,
                           pct_range=seq(0.5,1,by=0.1),
                           percent_dataframe_filename="pct_estimate_df.R") {
  
  load(sprintf("%s\\%s\\%s", session_path, likelihood_path, percent_dataframe_filename))
  likelihood_vec <- final_df[,"likelihood"]
  
  if (just_likelihood) {
    return(final_df)
  }
  names(likelihood_vec) <- rev(pct_range)
  estimated_pct <- as.numeric(names(which.max(likelihood_vec)))  
  
  extr_df <- extract_simulation_params(session_path, estimation_path = estimation_path, pct_range=pct_range)
  
  return(list(params=extr_df[as.character(estimated_pct),],pct=estimated_pct))
}

extract_simulation_params <- function(path, 
                                      plot=F, 
                                      absolute_min=T, 
                                      estimation_path="KS_simulations_estimate", 
                                      pct_range = seq(0.5, 1, by=0.1),
                                      file_name="pred_df.R") {
  
  print(sprintf("Estimating from %s", estimation_path))
  
  df_final_mean <- c()
  df_final_min <- c()
  
  for (true_p in pct_range) { 
    pred_df_path <- sprintf("%s/%s/param_fit_%.3f\\\\%s", path, estimation_path, true_p, file_name)
    
    load(pred_df_path)
    predicted <- t(predicted_df[,5:24])
    X_star <- predicted_df[,1:4]
    
    if (plot) {
      contour_maps(predicted, X_star, folder=sprintf("%s/%s/param_fit_%.3f/", path, estimation_path, true_p))
    }
    
    if (absolute_min) {
      df_final_min <- rbind(df_final_min,
                            X_star[which.min(colMeans(predicted)),])
      df_final_mean <- rbind(df_final_mean,
                             X_star[which.min(apply(predicted, 2, min)),])
      next
    }
    
    
    values <- apply(X_star, 2, unique)
    ind_mat <-  combn(1:len(values), 2)
    
    min_matrices_list <- 
      lapply(1:ncol(ind_mat), 
             function(col_idx) {
               ind <- ind_mat[,col_idx]
               
               v1 <- unlist(values[ind[1]])
               v2 <- unlist(values[ind[2]])
               
               
               min_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                                      nrow=len(v1))
               for (i in 1:len(v1))  {
                 for (j in 1:len(v2)) {
                   min_cont_mat[i,j] <- 
                     min(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
                 }
               }
               
               rownames(min_cont_mat) <- round(v1, digits=3)
               colnames(min_cont_mat) <- round(v2, digits=3)
               return(min_cont_mat)
             })
    
    mean_matrices_list <- 
      lapply(1:ncol(ind_mat), 
             function(col_idx) {
               ind <- ind_mat[,col_idx]
               
               v1 <- unlist(values[ind[1]])
               v2 <- unlist(values[ind[2]])
               
               
               mean_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                                       nrow=len(v1))
               for (i in 1:len(v1))  {
                 for (j in 1:len(v2)) {
                   mean_cont_mat[i,j] <- 
                     mean(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
                 }
               }
               rownames(mean_cont_mat) <- round(v1, digits=3)
               colnames(mean_cont_mat) <- round(v2, digits=3)
               return(mean_cont_mat)
             })
    
    
    
    final_params_min <- rep(0, times=len(values))
    final_params_mean <- rep(0, times=len(values))
    names(final_params_mean) <- names(values)
    names(final_params_min) <- names(values)
    for (mt_idx in 1:len(min_matrices_list)) {
      mt <- min_matrices_list[[mt_idx]]
      mean_mt <- mean_matrices_list[[mt_idx]]
      
      ind <- ind_mat[,mt_idx]
      v1 <- unlist(values[ind[1]])
      v2 <- unlist(values[ind[2]])
      
      row_min <- which.min(apply(mt, 1, min))
      col_min <- apply(mt, 1, which.min)[row_min]
      
      mean_row_min <- which.min(apply(mean_mt, 1, min))
      mean_col_min <- apply(mean_mt, 1, which.min)[mean_row_min]
      
      
      final_params_min[ind[1]] <- final_params_min[ind[1]] +  v1[row_min]
      final_params_min[ind[2]] <- final_params_min[ind[2]] +  v2[col_min]
      
      final_params_mean[ind[1]] <- final_params_mean[ind[1]] +  v1[mean_row_min]
      final_params_mean[ind[2]] <- final_params_mean[ind[2]] +  v2[mean_col_min]
      #return(c(v1[row_min], v2[col_min]))
    }
    
    final_params_min <- final_params_min / (len(values) - 1)
    final_params_mean <- final_params_mean / (len(values) - 1)
    
    df_final_min <- rbind(df_final_min, final_params_min)
    df_final_mean <- rbind(df_final_mean, final_params_mean)
  }
  
  rownames(df_final_mean) <- pct_range
  rownames(df_final_min) <- pct_range
  
  return(df_final_min)
}




get_switchers <- function(path, session, pval_threshold=0.05, num_shuffles=1000) {
  
  tmp <- get_spike_train_and_stim_trace_from_path(path, session)
  
  if (grepl("CA1", path)) {
    ext = str_split(str_split(path, "CA1\\\\")[[1]][[2]], "\\\\Matrices")[[1]][1]
  } else {
    ext = str_split(str_split(path, "CA3\\\\")[[1]][[2]], "\\\\Matrices")[[1]][1]
  }
  
  if (grepl("Right", path)) {
    dir = "right"
  } else {
    dir = "left"
  }
  
  ext = sprintf("%s_%s_session_%d", ext, dir, session)
  dir.create(sprintf("~\\switchers\\%s\\", ext))
  
  spike_train <- tmp[[1]]
  stim_trace <- tmp[[2]]
  
  ### Generate a spike train
  processed_generated <- preprocess_spike_train(spike_train, stim_trace, verbose=T)
  generated_cells_active_st <- processed_generated$working_cells_spike_train
  generated_firing_rate <- processed_generated$working_firing_rate
  
  
  tmp <- compute_tuning(generated_cells_active_st, stim_trace)
  gen_stim_prob <- tmp[[1]]
  gen_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
  
  #gen_shuff <- function(spike_train, idx,st="c") { return(t(array_reshape(generate_shuffled_spike_train(spike_train, idx, 200, st, 1), c(ncol(spike_train),200))))}
  simulated_place_cell_pval <- 
    sapply(c(1:nrow(generated_cells_active_st)),
           function(i, var2)
           {compute_place_signif(i,
                                 generated_cells_active_st,
                                 stim_trace,
                                 generated_firing_rate,
                                 generated_SI,
                                 shuffle_type=var2,
                                 verbose=T,
                                 num_shuffles = num_shuffles)},
           var2="c")
  
  simulated_place_cell_percentage <- sum(simulated_place_cell_pval < pval_threshold) / len(simulated_place_cell_pval)
  
  
  simulated_place_cell_pval_r <- 
    sapply(c(1:nrow(generated_cells_active_st)),
           function(i, var2)
           {compute_place_signif(i,
                                 generated_cells_active_st,
                                 stim_trace,
                                 generated_firing_rate,
                                 generated_SI,
                                 shuffle_type=var2,
                                 verbose=T,
                                 num_shuffles = num_shuffles)},
           var2="r")
  
  simulated_place_cell_percentage_r <- sum(simulated_place_cell_pval_r < pval_threshold) / len(simulated_place_cell_pval_r)
  
  switchers <- which((simulated_place_cell_pval_r < pval_threshold) & (simulated_place_cell_pval > pval_threshold))
  
  
  switchers_tuning <- gen_tuning_curve[switchers,];
  
  switcher_rank <- order(simulated_place_cell_pval_r[switchers] + 1 - simulated_place_cell_pval[switchers])
  
  
  # 
  # idx <- switchers[switcher_rank[30]]
  # 
  # switch_c <- gen_shuff(processed_generated$working_cells_spike_train,
  #                       idx, 
  #                       "c")
  # 
  # switch_r <- gen_shuff(processed_generated$working_cells_spike_train,
  #                       idx, 
  #                       "r")
  # 
  # spike_trains <- list(processed_generated$working_cells_spike_train[idx,],
  #                      switch_r[sample(1:200,1),],
  #                      switch_c[sample(1:200,1),])
  # 
  # plot_list <- list()
  # spike_trains <- lapply(sample(1:nrow(processed_generated$working_cells_spike_train), 20), function(nr) {processed_generated$working_cells_spike_train[nr,]})
  # for (st in spike_trains) {
  #   
  #   firing_ind <- which(st > 0)
  #   cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
  #                         Times=run_df$Time[firing_ind])
  #   
  #   print(len(firing_ind))
  #   grundf <-
  #     ggplot(run_df,aes(x=Time, y=Position)) + 
  #     
  #     geom_line(col="gray30") + 
  #     theme_light() +  
  #     base_plot_theme +
  #     xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
  #     ylim(0,96) +
  #     ylab("") +
  #     xlab("Time") +
  #     geom_point(data=cell_df, aes(x=Times,y=Positions),
  #                fill="red",
  #                color="red",
  #                size=3) + 
  #     coord_flip()
  #   
  #   pdf(sprintf("%s\\switcher_%d.pdf",write_path, idx), height=1.8 * 1.6, width=1.9 * 3); plot(gf); dev.off()
  #   
  #   #plot_list <- append(plot_list, list(grundf))
  #   
  # }
  # 
  # plot_list$nrow = 1
  # gf <- do.call(grid.arrange, plot_list)
  run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * 20), 
                       Position=stim_trace * 4) 
  switchers_st <- generated_cells_active_st[switchers,]
  
  for (idx in 1:nrow(switchers_st)) {
    
    firing_ind <- which(switchers_st[idx,] > 0)
    cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
                          Times=run_df$Time[firing_ind])
    
    print(len(firing_ind))
    grundf <-
      ggplot(run_df,aes(x=Time, y=Position)) + 
      
      geom_line(col="gray30") + 
      theme_light() +  
      base_plot_theme +
      xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
      ylim(0,96) +
      ylab("") +
      xlab("Time") +
      geom_point(data=cell_df, aes(x=Times,y=Positions),
                 fill="red",
                 color="red",
                 size=3) + 
      coord_flip()
    
    png(sprintf("~\\switchers\\%s\\switcher_%d.png", ext, processed_generated$ind[switchers][idx]), height=1.8 * 1.6, width=1.9, units="in", res=100); plot(grundf); dev.off()
    
    #plot_list <- append(plot_list, list(grundf))
    
  }
}


equalize_prior <- function(spike_train, stim_trace, verbose=F) {
  equalize_bins <- c(1:3,22:24)
  indices <- 1:len(stim_trace)
  
  
  final_ind <- rep(T, times=len(stim_trace))
  
  for (ebin in equalize_bins) {
    print(ebin)
    
    num_to_draw <- gen_num_frames(stim_trace)
    if (num_to_draw > len(which(stim_trace == ebin))) {
      if (verbose) {
        print(sprintf("No need to equalize bin %d", ebin))
      }
      next
    }
    
    sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
    tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin       
    final_ind <- final_ind & tmp
  }
  
  spike_train <- spike_train[,final_ind]
  stim_trace <- stim_trace[final_ind]
  
  return(list(spike_train=spike_train, 
              stim_trace=stim_trace))
}

get_place_cells <- function(spike_train, stim_trace, pval_threshold=0.05, num_shuffles=500, cyclic_only=F, active_bins_threshold=5, verbose=F) {
  ### Generate a spike train
  processed_generated <- preprocess_spike_train(spike_train, stim_trace, active_bins_threshold=active_bins_threshold, verbose=verbose)
  active_st <- processed_generated$working_cells_spike_train
  firing_rate <- processed_generated$working_firing_rate
  
  
  tmp <- compute_tuning(active_st, stim_trace)
  stim_prob <- tmp[[1]]
  tuning_curve <- tmp[[2]]
  rm(tmp)
  
  SI <- compute_SI(stim_prob, tuning_curve, firing_rate)
  
  cyclic_place_cell_pval <- 
    sapply(c(1:nrow(active_st)),
           function(i, var2)
           {compute_place_signif(i,
                                 active_st,
                                 stim_trace,
                                 firing_rate,
                                 SI,
                                 shuffle_type=var2,
                                 verbose=verbose,
                                 num_shuffles = num_shuffles)},
           var2="c")
  
  cyclic_place_cell_percentage <- sum(cyclic_place_cell_pval < pval_threshold) / len(cyclic_place_cell_pval)
  
  if (cyclic_only) {
    
    cyclic_non_pcs <- which((cyclic_place_cell_pval > pval_threshold))
    cyclic_pcs <- which((cyclic_place_cell_pval < pval_threshold))
  
    
    return(list(ind=processed_generated$ind,
                cyclic_pval=cyclic_place_cell_pval,
                cyclic_non_pcs=cyclic_non_pcs,
                cyclic_pcs=cyclic_pcs,
                cyclic_fraction=cyclic_place_cell_percentage,
                SI=SI,
                FR=firing_rate))
  }
  
  random_place_cell_pval <- 
    sapply(c(1:nrow(active_st)),
           function(i, var2)
           {compute_place_signif(i,
                                 active_st,
                                 stim_trace,
                                 firing_rate,
                                 SI,
                                 shuffle_type=var2,
                                 verbose=verbose,
                                 num_shuffles = num_shuffles)},
           var2="r")
  
  random_place_cell_percentage <- sum(random_place_cell_pval < pval_threshold) / len(cyclic_place_cell_pval)
  
  switchers <- which((random_place_cell_pval < pval_threshold) & (cyclic_place_cell_pval > pval_threshold))
  cyclic_non_pcs <- which((cyclic_place_cell_pval > pval_threshold))
  cyclic_pcs <- which((cyclic_place_cell_pval < pval_threshold))
  random_non_pcs <- which((random_place_cell_pval > pval_threshold))
  random_pcs <- which((random_place_cell_pval < pval_threshold))
  
  return(list(ind=processed_generated$ind,
              random_pval=random_place_cell_pval,
              cyclic_pval=cyclic_place_cell_pval,
              switchers=switchers,
              random_non_pcs=random_non_pcs,
              random_pcs=random_pcs,
              cyclic_non_pcs=cyclic_non_pcs,
              cyclic_pcs=cyclic_pcs,
              random_fraction=random_place_cell_percentage,
              cyclic_fraction=cyclic_place_cell_percentage,
              SI=SI,
              FR=firing_rate))
  
}


preprocess_spike_train_no_threshold <- function(spike_train, stim_trace, verbose=T) {
  time_bins_activity_per_cell <- 
    sapply(1:nrow(spike_train), function(i) {sum(as.numeric(spike_train[i,] > 0))})
  
  activity_per_cell <- rowMeans(spike_train) / dt
  
  # Filter out sufficiently active cells 
  working_cells_spike_train <- spike_train
  working_firing_rate <- activity_per_cell
  working_time_bins_per_cells <- time_bins_activity_per_cell
  
  return(list(working_cells_spike_train=working_cells_spike_train,
              working_firing_rate=working_firing_rate,
              working_time_bins_per_cells=working_time_bins_per_cells,
              ind=1:nrow(spike_train)))
}

euclidean <- function(a, b) sqrt(sum((a - b)^2))

tuning_similarity <- function(tuning_mat_a, 
                              tuning_mat_b,
                              similarity_measure="cor") {
  
  if (similarity_measure == "cor") {
    
    sim_vec <- unlist(lapply(1:nrow(tuning_mat_a),
                            function(index) 
                              {
                              cor(tuning_mat_a[index,],
                                  tuning_mat_b[index,])
                            }))
  } else if (similarity_measure == "peaks") {
    sim_vec <- apply(tuning_mat_a, 1, which.max) - 
                  apply(tuning_mat_b, 1, which.max)
  } else if (similarity_measure == "eucl") {
    sim_vec <- unlist(lapply(1:nrow(tuning_mat_a),
                             function(index) 
                             {
                               euclidean(tuning_mat_a[index,],
                                   tuning_mat_b[index,])
                             }))
  } else if (similarity_measure == "cosine") {
    sim_vec <- unlist(lapply(1:nrow(tuning_mat_a),
                             function(index) 
                             {
                               cosine(tuning_mat_a[index,],
                                   tuning_mat_b[index,])
                             }))
  } else if (similarity_measure == "wasserstein") {
    normed_tuning_mat_a <- tuning_mat_a / (rowSums(tuning_mat_a) + 10^-30)
    normed_tuning_mat_b <- tuning_mat_b / (rowSums(tuning_mat_b) + 10^-30)
    
    sim_vec <- unlist(lapply(1:nrow(normed_tuning_mat_b),
                             function(index) 
                             {
                              vec_a <- rep(1:ncol(normed_tuning_mat_a), 
                                           times=(normed_tuning_mat_a[index,] * 1000))
                              vec_b <- rep(1:ncol(normed_tuning_mat_b), 
                                           times=(normed_tuning_mat_b[index,] * 1000))
                              
                              if (len(vec_a) == 0 || len(vec_b) == 0) {
                                return(NA)
                              }
                               return(wasserstein_dist(vec_a,vec_b))
                             }))
  }
  return(sim_vec)
}

get_num_traversals <- function(working_path, old=F, equalize_frames=F, simulated_spike_trains=F, sessions_to_use=c(1:16)) {
  
  path <- working_path
  
  if (old) {
    m <- readMat(paste(path, "\\stimulus_trace.mat", sep=""))
    stim_trace <- m$stimulus.trace
    
    
    m <- readMat(paste(path, spike_train_filename, sep=""))
    spike_train <- m$spike.train
  } else {  
    m <- readMat(paste(path, stimulus_trace_filename, sep=""))
    
    if(len(grep("Left", path)) > 0) {
      stim_trace <- m$position.left.trials
    } else {
      stim_trace <- m$position.right.trials
    }
    
    
    m <- readMat(paste(path, spike_train_filename, sep=""))
    
    if(len(grep("Left", path)) > 0) {
      spike_train <- m$spikes.left.trials
    } else {
      spike_train <- m$spikes.right.trials
    }
  }
  
  tmp_path <- path
  tmp_spike_train <- spike_train
  tmp_stimulus_trace <- stim_trace
  
  l = len(tmp_spike_train)
  
  traversals_mat <- list()
  sliced_ind <- list()
  
  for (i in (sessions_to_use)) {
    # Restore spike train / stim trace list
    spike_train <- tmp_spike_train
    stim_trace <- tmp_stimulus_trace
    spike_train <- spike_train[[i]][[1]]
    spike_train <- t(as.matrix(spike_train))
    stim_trace <- stim_trace[[i]][[1]]
    stim_trace <-  as.vector(stim_trace)
    
    path <- tmp_path
    
    
    if (!old && equalize_frames ) {
      equalize_bins <- c(1:3,22:24)
      indices <- 1:len(stim_trace)
      
      
      final_ind <- rep(T, times=len(stim_trace))
      
      for (ebin in equalize_bins) {
        print(ebin)
        
        num_to_draw <- gen_num_frames(stim_trace)
        if (num_to_draw > len(which(stim_trace == ebin))) {
          print(sprintf("No need to equalize bin %d", ebin))
          next
        }
        
        sampled_ind <- sample(which(stim_trace==ebin), num_to_draw)
        tmp <- indices %in% sampled_ind | stim_trace[indices] != ebin      
        final_ind <- final_ind & tmp
      }
      
      spike_train <- spike_train[,final_ind]
      stim_trace <- stim_trace[final_ind]
      sliced_ind <- append(sliced_ind, list(final_ind))
    }
    
    trav <- get_traversales(stim_trace)
    traversals_mat <- append(traversals_mat, list(trav))
    
  }
  
  return(list(traversals_ind=traversals_mat,
              num_of_trav=sapply(traversals_mat, length),
              sliced_ind=sliced_ind))
}

get_traversales <- function(stim_trace, nbins=5) {
  max_bin <- unique(sort(stim_trace))[(24-nbins):24]
  min_bin <- unique(sort(stim_trace))[1:nbins]
  min_ind <- which(stim_trace %in% min_bin)
  max_ind <- which(stim_trace %in% max_bin)
  reward_ind <- c(min_ind, max_ind)
  
  if (min_ind[1] < max_ind[1]) {
    names(reward_ind) <- c(rep(1, times=len(min_ind)),
                           rep(0, times=len(max_ind)))
  } else {
    names(reward_ind) <- c(rep(0, times=len(min_ind)),
                           rep(1, times=len(max_ind)))        
  }
  
  reward_ind <- sort(reward_ind)
  diff_reward <- diff(as.numeric(names(reward_ind)))
  
  nth_traversal_frame_ind <- which(diff_reward == -1)
  return(reward_ind[nth_traversal_frame_ind])
}
