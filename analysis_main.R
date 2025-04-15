library("reticulate")
library("testit")
library("dplyr")
library("plyr")
library("OneR")
library("stringr")
library("R.matlab")
library("ggplot2")
library("ggridges")
library('latex2exp')
library("gridExtra")
library("reshape2")
library("RColorBrewer")
library("lsa")

# Enter your designated path for the data
remote_path = "Y:/livneh/itayta/AllActiveCells/new_f/raw_data/"
#base_path = "C:\\Users\\itayta\\Desktop\\raw_data\\raw_data_07_02_23"
base_path = "/Users/itayta/Downloads/yaniv_place_cells_project/raw_data/raw_data_07_02_23"



mice_ids <- c("C5M1", "C6M3", "C6M4", "C8M2", "C14M4", "C15M2", "C23M4", "C24M3", "C24M4")

data_path_1 <- sprintf("%s/CA1/C5M1/Matrices/Right", base_path)
data_path_2 <- sprintf("%s/CA1/C5M1/Matrices/Left", base_path)
data_path_3 <- sprintf("%s/CA1/C6M3/Matrices/Right", base_path)
data_path_4 <- sprintf("%s/CA1/C6M3/Matrices/Left", base_path)
data_path_5 <- sprintf("%s/CA1/C6M4/Matrices/Right", base_path)
data_path_6 <- sprintf("%s/CA1/C6M4/Matrices/Left", base_path)
data_path_7 <- sprintf("%s/CA1/C8M2/Matrices/Right", base_path)
data_path_8 <- sprintf("%s/CA1/C8M2/Matrices/Left", base_path)
data_path_9 <- sprintf("%s/CA3/C14M4/Matrices/Right", base_path)
data_path_10 <- sprintf("%s/CA3/C14M4/Matrices/Left", base_path)
data_path_11 <- sprintf("%s/CA3/C15M2/Matrices/Right", base_path)
data_path_12 <- sprintf("%s/CA3/C15M2/Matrices/Left", base_path)
data_path_13 <- sprintf("%s/CA3/C23M4/Matrices/Right", base_path)
data_path_14 <- sprintf("%s/CA3/C23M4/Matrices/Left", base_path)
data_path_15 <- sprintf("%s/CA3/C24M3/Matrices/Right", base_path)
data_path_16 <- sprintf("%s/CA3/C24M3/Matrices/Left", base_path)
data_path_17 <- sprintf("%s/CA3/C24M4/Matrices/Right", base_path)
data_path_18 <- sprintf("%s/CA3/C24M4/Matrices/Left", base_path)

all_data_paths <- c(data_path_1,  data_path_2,  data_path_3,  data_path_4,  data_path_5,  data_path_6,
                    data_path_7,  data_path_8,  data_path_9,  data_path_10, data_path_11, data_path_12,
                    data_path_13, data_path_14, data_path_15, data_path_16, data_path_17, data_path_18)

SUBSAMPLE_DURATION = 1
BY_ACTIVITY = 2
SIMULATION_PARAM_ESTIMATION = 3
LIKELIHOOD_PLACE_CELL_PCT = 4
SIMULATED_VS_REAL_BY_SAMPLE = 5
PVAL_BY_SUBSAMPLE_DURATION = 6
TUNING_CORR_WITHIN = 7
SLIDING_TUNING_PROP=8
PROB_BY_POS_ANALYSIS = 9
DILUTE_PCS_ANALYSIS=10
PEAK_DIFF_BY_RUNS=11
SLIDING_TUNING_PROP_LIRON=12
SUBSAMPLE_DURATION_COR=13
BY_ACTIVITY_COR=14

stimulus_trace_filename = "/stim_trace.mat"
spike_train_filename = "/spike_train.mat"
registration_mat_filename = "registration.mat"

#####
main <- function(working_path, 
                 old=F, 
                 equalize_prior=T, 
                 fit_simulation=F,
                 naive_simulation=F, 
                 dual_simulations_original=F,
                 remove_edges=F,
                 analysis_type=SUBSAMPLE_DURATION, 
                 sessions_to_use=1:16) {
  
  path <- working_path
  
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
  
  l = len(tmp_spike_train)
  
  for (session in sessions_to_use) {
    # Restore spike train / stim trace list
    spike_train <- tmp_spike_train
    stim_trace <- tmp_stimulus_trace
    spike_train <- spike_train[[session]][[1]]
    spike_train <- t(as.matrix(spike_train))
    stim_trace <- stim_trace[[session]][[1]]
    stim_trace <-  as.vector(stim_trace)
    
    path <- tmp_path
    
    if (naive_simulation) {
      print("Simulating spike trains! (NAIVE)")
      
      fr <- rowMeans(spike_train) / dt
      processed_real <- preprocess_spike_train(spike_train, stim_trace)
      
      true_cells_spike_train <- processed_real$working_cells_spike_train
      true_firing_rate <- processed_real$working_firing_rate
      true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
      
      simulated_tuning_curve <-
        generate_tuning_curves_cost(n = nrow(spike_train),
                                    percentage = 0.8,
                                    average_width = 0.1, 
                                    sd_width = 0.005,
                                    fixed_fr=fr,
                                    noise=0.02,
                                    double_peak_pct = 0.4,
                                    plot=T)
      
      pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                              true_time_bins_per_cells,
                                              stim_trace,
                                              verbose = T,
                                              jump=0.05)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)
      spike_train <- generated_spike_train
      
      path <- sprintf("%s/naive_simulation/", path)
      dir.create(path)
      
    } 
    
    if (!old && equalize_prior) {
      dir.create(sprintf("%s/equalized/", path))
      path <- sprintf("%s/equalized/session_%d", path, session)  
      tmp <- equalize_prior(spike_train, stim_trace, verbose=T)
      
      spike_train <- tmp$spike_train
      stim_trace <- tmp$stim_trace
      
      rm(tmp)
      
      
    } else {
      path <- sprintf("%s/session_%d", path, session)
    }
    
    if (!old && remove_edges) {
      print("REMOVING EDGES!")
      path <- sprintf("%s/edge_free", path)
      non_edges_ind <- which(stim_trace %in% c(3:22))
      stim_trace <- stim_trace[non_edges_ind]
      spike_train <- spike_train[,non_edges_ind]
      stim_trace <- stim_trace - 2 # Refer to bin 3 as bin 1
    }
    
    if (fit_simulation) {
      print("Using fit parameters for simulation")
      fr <- rowMeans(spike_train) / dt
      processed_real <- preprocess_spike_train(spike_train, stim_trace)
      
      true_cells_spike_train <- processed_real$working_cells_spike_train
      true_firing_rate <- processed_real$working_firing_rate
      true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
      
      #params <- get_fit_params(path, estimation_path="simulations_2_new_new", likelihood_path = "simulations_3_new", pct_range=seq(0.5,1,by=0.1))
      params <- get_fit_params(path )
      
      print(params)
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
                                              verbose = T,
                                              jump=0.05)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)
      
      orig_st <- spike_train
      orig_path <- path
      spike_train <- generated_spike_train
      
      path <- sprintf("%s/fit_simulation/", path)
    }
    
    dir.create(path)
    dir.create(sprintf("%s/subsampled_matrices", path))
    
    if (dual_simulations_original) {
      dir.create(sprintf("%s/subsampled_matrices", orig_path))
    }
    
    print(sprintf("#### runnig session %d ####", session))
    if (analysis_type == SUBSAMPLE_DURATION) {

      main_calc_place_fraction_by_subsamples(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_place_fraction_by_subsamples(orig_st, stim_trace, orig_path)
      }
    } else if (analysis_type == BY_ACTIVITY){

      main_calc_place_fraction_by_activity(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_place_fraction_by_activity(orig_st, stim_trace, orig_path)
      }
    } else if (analysis_type == PVAL_BY_SUBSAMPLE_DURATION) {
      main_calc_pvalue_by_subsamples(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_pvalue_by_subsamples(orig_st, stim_trace, orig_path)
      }
      
    } else if (analysis_type == TUNING_CORR_WITHIN) {
      main_calc_tuning_cor_by_subsamples(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_tuning_cor_by_subsamples(orig_st, stim_trace, orig_path)
      }
      
    } else if (analysis_type == SIMULATION_PARAM_ESTIMATION) {
      main_estimate_simulation_parameters(spike_train, stim_trace, remove_edges = remove_edges, path)
    } else if (analysis_type == LIKELIHOOD_PLACE_CELL_PCT) {
      main_calculate_place_cell_fraction_likelihood(spike_train, stim_trace, remove_edges = remove_edges, path)
    } else if (analysis_type == SIMULATED_VS_REAL_BY_SAMPLE) {
      main_overlaid_simulation_real_by_sample_size(spike_train, stim_trace, path)
    }  else if (analysis_type == SLIDING_TUNING_PROP) {
      main_calc_sliding_tuning_similarity(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_sliding_tuning_similarity(orig_st, stim_trace, orig_path)
      }
      
    }  else if (analysis_type == PROB_BY_POS_ANALYSIS) {
      main_calc_prob_by_position_firing(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_prob_by_position_firing(orig_st, stim_trace, orig_path)
      }
    } else if (analysis_type == DILUTE_PCS_ANALYSIS) {
      main_calc_diluted_place_cell_fraction(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_diluted_place_cell_fraction(orig_st, stim_trace, orig_path)
      }
    } else if (analysis_type == PEAK_DIFF_BY_RUNS) {
      main_calc_peak_movement_runs(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_peak_movement_runs(orig_st, stim_trace, orig_path)
      }
    } else if (analysis_type == SLIDING_TUNING_PROP_LIRON) {

      main_calc_sliding_tuning_similarity_liron(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_sliding_tuning_similarity_liron(orig_st, stim_trace, orig_path)
      }
      
      
    } else if (analysis_type == SUBSAMPLE_DURATION_COR) {
      
      main_calc_place_fraction_by_subsamples_tuning_cor(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_place_fraction_by_subsamples_tuning_cor(orig_st, stim_trace, orig_path)
      }
    } else if (analysis_type == BY_ACTIVITY_COR) {
      
      main_calc_place_fraction_by_activity_cor(spike_train, stim_trace, path)
      
      if (fit_simulation && dual_simulations_original){
        print(orig_path)
        main_calc_place_fraction_by_activity_cor(orig_st, stim_trace, orig_path)
      }
    } 
    
    
  }
}



main_both_directions <-  function(mice_id, 
                                  nreps=5, 
                                  pval_t=0.05,
                                  pvals=c(0.001, 0.01,0.025,0.05,0.1,0.15),
                                  old=F, 
                                  equalize_prior=T, 
                                  num_shuffles=500,
                                  use_MI=F,
                                  use_PP=F,
                                  ext="",
                                  equalize_traversals=F,
                                  sessions_to_use = 1:16) {
  
  paths_both <- all_data_paths[grep(mice_id, all_data_paths)]
  path_left <- paths_both[1]
  path_right <- paths_both[2]
  
  m <- readMat(paste(path_left, stimulus_trace_filename, sep=""))
  
  if(len(grep("Left", path_left)) > 0) {
    stim_trace_left <- m$position.left.trials
  } else {
    stim_trace_left <- m$position.right.trials
  }
  
  
  m <- readMat(paste(path_left, spike_train_filename, sep=""))
  
  if(len(grep("Left", path_left)) > 0) {
    spike_train_left <- m$spikes.left.trials
  } else {
    spike_train_left <- m$spikes.right.trials
  }
  
  m <- readMat(paste(path_right, stimulus_trace_filename, sep=""))
  
  if(len(grep("Left", path_right)) > 0) {
    stim_trace_right <- m$position.left.trials
  } else {
    stim_trace_right <- m$position.right.trials
  }
  
  
  m <- readMat(paste(path_right, spike_train_filename, sep=""))
  
  if(len(grep("Left", path_right)) > 0) {
    spike_train_right <- m$spikes.left.trials
  } else {
    spike_train_right <- m$spikes.right.trials
  }
  
  tmp_path_left <- path_left
  tmp_spike_train_left <- spike_train_left
  tmp_stimulus_trace_left <- stim_trace_left
  
  tmp_path_right <- path_right
  tmp_spike_train_right <- spike_train_right
  tmp_stimulus_trace_right <- stim_trace_right
  
  final_df_list_cyclic <- list()
  final_df_list_rand <- list()
  
  for (repetition in 1:nreps) {
    if (!old && equalize_traversals) {
      
      ntraversals_left <- get_num_traversals(path_left, old = F, equalize_frames = T)
      ntraversals_right <- get_num_traversals(path_right, old = F, equalize_frames = T)
      
      nleft_min <- min(ntraversals_left$num_of_trav) - 1
      nright_min <- min(ntraversals_right$num_of_trav) - 1
    }
  
    result_df_cyclic <- list()
    result_df_rand <- list()
    
    
    for (session in (sessions_to_use)) {
      print(sprintf("### Session %d Rep %d", session, repetition))
      
      # Restore spike train / stim trace list
      spike_train_right <- tmp_spike_train_right
      stim_trace_right <- tmp_stimulus_trace_right
      spike_train_right <- spike_train_right[[session]][[1]]
      spike_train_right <- t(as.matrix(spike_train_right))
      stim_trace_right <- stim_trace_right[[session]][[1]]
      stim_trace_right <-  as.vector(stim_trace_right)
      
      spike_train_left <- tmp_spike_train_left
      stim_trace_left <- tmp_stimulus_trace_left
      spike_train_left <- spike_train_left[[session]][[1]]
      spike_train_left <- t(as.matrix(spike_train_left))
      stim_trace_left <- stim_trace_left[[session]][[1]]
      stim_trace_left <-  as.vector(stim_trace_left)    
      
      path_left <- tmp_path_left
      path_right <- tmp_path_right
      
      if (!old && equalize_prior) {
        equalize_bins <- c(1:3,22:24)
        indices <- 1:len(stim_trace_left)
        
        
        final_ind <- rep(T, times=len(stim_trace_left))
        
        for (ebin in equalize_bins) {
          print(ebin)
          
          num_to_draw <- gen_num_frames(stim_trace_left)
          if (num_to_draw > len(which(stim_trace_left == ebin))) {
            print(sprintf("No need to equalize bin %d", ebin))
            next
          }
          
          sampled_ind <- sample(which(stim_trace_left==ebin), num_to_draw)
          tmp <- indices %in% sampled_ind | stim_trace_left[indices] != ebin      
          final_ind <- final_ind & tmp
        }
        
        spike_train_left <- spike_train_left[,final_ind]
        stim_trace_left <- stim_trace_left[final_ind]
        
        print(dim(spike_train_left))
        print(len(stim_trace_left))
        
        indices <- 1:len(stim_trace_right)
        
        final_ind <- rep(T, times=len(stim_trace_right))
        
        for (ebin in equalize_bins) {
          print(ebin)
          
          num_to_draw <- gen_num_frames(stim_trace_right)
          if (num_to_draw > len(which(stim_trace_right == ebin))) {
            print(sprintf("No need to equalize bin %d", ebin))
            next
          }
          
          sampled_ind <- sample(which(stim_trace_right==ebin), num_to_draw)
          tmp <- indices %in% sampled_ind | stim_trace_right[indices] != ebin      
          final_ind <- final_ind & tmp
        }
        
        spike_train_right <- spike_train_right[,final_ind]
        stim_trace_right <- stim_trace_right[final_ind]  
        
        print(dim(spike_train_right))
        print(len(stim_trace_right))
        
      } else if (!old && equalize_traversals) {
        if (ntraversals_left$num_of_trav[[session]] > nleft_min) {
          
          sampled_ind <- sample(1:(len(ntraversals_left$traversals_ind[[session]]) - 1), nleft_min - 1)
          sampled_ind <- sort(sampled_ind)
          sampled_ind_2 <-
            lapply(sampled_ind,
                   function(si) { ntraversals_left$traversals_ind[[session]][si:(si+1)] })
          sampled_ind_f <- lapply(sampled_ind_2, function(sc) { return((sc[1]):(sc[2])) })
          
          final_ind <- unlist(sampled_ind_f)
          print(sprintf("Got %d runs", len(sampled_ind_f)))
          
          spike_train_left <- spike_train_left[,which(ntraversals_left$sliced_ind[[session]])[final_ind]]
          stim_trace_left <- stim_trace_left[which(ntraversals_left$sliced_ind[[session]])[final_ind]]  
          
          print(dim(spike_train_left))
          print(len(stim_trace_left))
        }
        
        if (ntraversals_right$num_of_trav[[session]] > nright_min) {
          
          sampled_ind <- sample(1:(len(ntraversals_right$traversals_ind[[session]]) - 1), nright_min - 1)
          sampled_ind <- sort(sampled_ind)
          sampled_ind_2 <-
            lapply(sampled_ind,function(si) {ntraversals_right$traversals_ind[[session]][si:(si+1)]})
          sampled_ind_f <- lapply(sampled_ind_2, function(sc) { return((sc[1]):(sc[2]))})
          
          final_ind <- unlist(sampled_ind_f)
          
          spike_train_right <- spike_train_right[,which(ntraversals_right$sliced_ind[[session]])[final_ind]]
          stim_trace_right <- stim_trace_right[which(ntraversals_right$sliced_ind[[session]])[final_ind]]  
          
          print(dim(spike_train_right))
          print(len(stim_trace_right))
        }
      } else {
        print("NOT EQUALIZING anything!")
      }
      
      fr <- rowMeans(spike_train_left) / dt
      processed_real_left <- preprocess_spike_train(spike_train_left, stim_trace_left)
      
      true_cells_spike_train <- processed_real_left$working_cells_spike_train
      true_firing_rate <- processed_real_left$working_firing_rate
      true_time_bins_per_cells <- processed_real_left$working_time_bins_per_cells
      
      tmp <- compute_tuning(true_cells_spike_train, stim_trace_left)
      stim_prob <- tmp[[1]]
      true_tuning_curve <- tmp[[2]]
      rm(tmp)
      
      if (use_MI) {
        true_MI <- compute_MI(stim_trace_left, true_cells_spike_train)
        print("##########")
        print(sum(is.nan(true_MI)))
        print(session)
        print(which(is.nan(true_MI)))
        print("##########")
        place_cell_pval_left_cyclic <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif_MI(i,
                                          true_cells_spike_train,
                                          stim_trace_left,
                                          true_MI,
                                          shuffle_type=var2,
                                          num_shuffles=num_shuffles)},
                 var2="c")
        
        place_cell_pval_left_rand <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif_MI(i,
                                          true_cells_spike_train,
                                          stim_trace_left,
                                          true_MI,
                                          shuffle_type=var2,
                                          num_shuffles=num_shuffles)},
                 var2="r")        
        
      } else {
        true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
        
        place_cell_pval_left_cyclic <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif(i,
                                       true_cells_spike_train,
                                       stim_trace_left,
                                       true_firing_rate,
                                       true_SI,
                                       shuffle_type=var2,
                                       num_shuffles=num_shuffles,
                                       verbose=F)},
                 var2="c")
        
        place_cell_pval_left_rand <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif(i,
                                       true_cells_spike_train,
                                       stim_trace_left,
                                       true_firing_rate,
                                       true_SI,
                                       shuffle_type=var2,
                                       num_shuffles=num_shuffles,
                                       verbose=F)},
                 var2="r")        
      }
      
      
      fr <- rowMeans(spike_train_right) / dt
      processed_real_right <- preprocess_spike_train(spike_train_right, stim_trace_right)
      
      true_cells_spike_train <- processed_real_right$working_cells_spike_train
      true_firing_rate <- processed_real_right$working_firing_rate
      true_time_bins_per_cells <- processed_real_right$working_time_bins_per_cells
      
      tmp <- compute_tuning(true_cells_spike_train, stim_trace_right)
      stim_prob <- tmp[[1]]
      true_tuning_curve <- tmp[[2]]
      rm(tmp)
      
      if (use_MI) {
        true_MI <- compute_MI(stim_trace_right, true_cells_spike_train)
        print("##########")
        print(sum(is.nan(true_MI)))
        print(session)
        print(which(is.nan(true_MI)))
        print("##########")
        
        place_cell_pval_right_cyclic <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif_MI(i,
                                          true_cells_spike_train,
                                          stim_trace_right,
                                          true_MI,
                                          shuffle_type=var2,
                                          num_shuffles=500)},
                 var2="c")
        
        place_cell_pval_right_rand <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif_MI(i,
                                          true_cells_spike_train,
                                          stim_trace_right,
                                          true_MI,
                                          shuffle_type=var2,
                                          num_shuffles=500)},
                 var2="r")
      } else {
        true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
        
        place_cell_pval_right_cyclic <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif(i,
                                       true_cells_spike_train,
                                       stim_trace_right,
                                       true_firing_rate,
                                       true_SI,
                                       shuffle_type=var2,
                                       num_shuffles=500,
                                       verbose=F)},
                 var2="c")
        
        
        place_cell_pval_right_rand <-
          sapply(c(1:nrow(true_cells_spike_train)),
                 function(i, var2)
                 {compute_place_signif(i,
                                       true_cells_spike_train,
                                       stim_trace_right,
                                       true_firing_rate,
                                       true_SI,
                                       shuffle_type=var2,
                                       num_shuffles=500,
                                       verbose=F)},
                 var2="r")
        
      }
      
      cell_indices <- unique(c(processed_real_left$ind, processed_real_right))
      
      pct_df_cyclic <- data.frame(matrix(rep(NA, times=nrow(spike_train_left) * 2), nrow=nrow(spike_train_left)))
      colnames(pct_df_cyclic) <- c("Left", "Right")
      
      pct_df_cyclic$Left[processed_real_left$ind] <- place_cell_pval_left_cyclic
      pct_df_cyclic$Right[processed_real_right$ind] <- place_cell_pval_right_cyclic
      
      pct_df_rand <- data.frame(matrix(rep(NA, times=nrow(spike_train_left) * 2), nrow=nrow(spike_train_left)))
      colnames(pct_df_rand) <- c("Left", "Right")
      
      pct_df_rand$Left[processed_real_left$ind] <- place_cell_pval_left_rand
      pct_df_rand$Right[processed_real_right$ind] <- place_cell_pval_right_rand
      
      for (pval_t in pvals) {
      percentage_vector_cyclic <-
        apply(pct_df_cyclic,
              1,
              function(r) {
                if(sum(is.na(r)) == 2) {
                  return(NA)
                }
                
                row_ind <- !is.na(r)
                return(sum(r[row_ind] < pval_t))  
                
              })  
      
      percentage_vector_rand <-
        apply(pct_df_rand,
              1,
              function(r) {
                if(sum(is.na(r)) == 2) {
                  return(NA)
                }
                
                row_ind <- !is.na(r)
                return(sum(r[row_ind] < pval_t))  
                
              })
      
      non_na_vec_cyclic <- percentage_vector_cyclic[!is.na(percentage_vector_cyclic)]
      non_na_vec_rand <- percentage_vector_rand[!is.na(percentage_vector_rand)]
      
      print("____________________________")
      print(sprintf("_____ pval = %f ____", pval_t))
      print(place_cell_pval_right_cyclic)
      print(place_cell_pval_left_cyclic)
      print(non_na_vec_cyclic)
      print(place_cell_pval_right_rand)
      print(place_cell_pval_left_rand)
      print(non_na_vec_rand)
      print("____________________________")
      
      result_df_cyclic[[as.character(pval_t)]] <-  
                  rbind(result_df_cyclic[[as.character(pval_t)]], 
                                                    c(sum(non_na_vec_cyclic > 0) / len(non_na_vec_cyclic),
                                                      sum(place_cell_pval_right_cyclic < pval_t) / len(place_cell_pval_right_cyclic),
                                                      sum(place_cell_pval_left_cyclic < pval_t) / len(place_cell_pval_left_cyclic),
                                                      
                                                      sum(non_na_vec_cyclic > 0) / len(percentage_vector_cyclic),
                                                      sum(place_cell_pval_right_cyclic < pval_t) / nrow(spike_train_right),
                                                      sum(place_cell_pval_left_cyclic < pval_t) / nrow(spike_train_left),
                                  
                                                      len(cell_indices) / nrow(spike_train_left),
                                                      len(place_cell_pval_right_cyclic) / nrow(spike_train_right),
                                                      len(place_cell_pval_left_cyclic) / nrow(spike_train_left)))
      
      
      
      result_df_rand[[as.character(pval_t)]] <-   
                                          rbind(result_df_rand[[as.character(pval_t)]], 
                                                c(sum(non_na_vec_rand > 0) / len(non_na_vec_rand),
                                                  sum(place_cell_pval_right_rand < pval_t) / len(place_cell_pval_right_rand),
                                                  sum(place_cell_pval_left_rand < pval_t) / len(place_cell_pval_left_rand),
                                                      
                                                  sum(non_na_vec_rand > 0) / len(non_na_vec_rand),
                                                  sum(place_cell_pval_right_rand < pval_t) / nrow(spike_train_right),
                                                  sum(place_cell_pval_left_rand < pval_t) / nrow(spike_train_left),
                                                      
                                                  len(cell_indices) / nrow(spike_train_left),
                                                  len(place_cell_pval_right_rand) / nrow(spike_train_right),
                                                  len(place_cell_pval_left_rand) / nrow(spike_train_left)))
      
      
      print(sum(non_na_vec_cyclic > 0) / len(non_na_vec_cyclic))
      print(sum(non_na_vec_rand > 0) / len(non_na_vec_rand))
      
      colnames(result_df_cyclic[[as.character(pval_t)]]) <- 
                                    c("All", "Right", "Left",
                                      "All_out_of_all", "Right_out_of_all", "Left_out_of_all",
                                      "All active", "Right active", "Left active")
      
      colnames(result_df_rand[[as.character(pval_t)]]) <-
                                  c("All", "Right", "Left",
                                    "All_out_of_all", "Right_out_of_all", "Left_out_of_all",
                                    "All active", "Right active", "Left active")
    }

      
    if (use_MI) {
      print(names(result_df_cyclic))
      final_df_list_cyclic <- append(final_df_list_cyclic, list(result_df_cyclic[["0.05"]]))
      final_df_list_rand <- append(final_df_list_rand, list(result_df_rand[["0.05"]]))
    } else {
      final_df_list_cyclic <- append(final_df_list_cyclic, list(result_df_cyclic))
      final_df_list_rand <- append(final_df_list_rand, list(result_df_rand))
    }
    
    }
  }
  
  
  
  if (use_MI) {
    cyclic_file_name = sprintf("percent_df_%s_%s_%s.R",
                               "cyclic",
                               ifelse(use_MI, "MI", "SI"),
                               ifelse(equalize_prior, "eq_prior", ifelse(equalize_traversals, "eq_trav", "reg")))
    
    rand_file_name = sprintf("percent_df_%s_%s_%s.R",
                             "rand",
                             ifelse(use_MI, "MI", "SI"),
                             ifelse(equalize_prior, "eq_prior", ifelse(equalize_traversals, "eq_trav", "reg")))
  } else {
  cyclic_file_name = sprintf("multiple_pvals_percent_df_%s_%s_%s.R",
                             "cyclic",
                             ifelse(use_MI, "MI", "SI"),
                             ifelse(equalize_prior, "eq_prior", ifelse(equalize_traversals, "eq_trav", "reg")))
  
  rand_file_name = sprintf("multiple_pvals_percent_df_%s_%s_%s.R",
                           "rand",
                           ifelse(use_MI, "MI", "SI"),
                           ifelse(equalize_prior, "eq_prior", ifelse(equalize_traversals, "eq_trav", "reg")))
  }
  
  save(file=sprintf("%s/%s", 
                    str_split(path_left, "Matrices")[[1]][[1]],
                    cyclic_file_name),
       final_df_list_cyclic)
  
  save(file=sprintf("%s/%s", 
                    str_split(path_left, "Matrices")[[1]][[1]],
                    rand_file_name),
       final_df_list_rand)
  
}



main_across_learning <- function (path, result_file_name="registered_cells_categories.Rda", old=F) {
  
  # Groups
  DORMENT = 0;
  NON_RANDOM_NON_CYCLIC = 1;
  RANDOM_NON_CYLIC = 2; # Swithcers
  NON_RANDOM_CYCLIC = 3;
  RANDOM_CYCLIC = 4;
  
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
  
  
  registration_path <- sprintf("%s%s",
                               str_split(path, "Matrices//")[[1]][1], registration_mat_filename)
  
  registration_mat <- readMat(registration_path)
  registration_mat <- registration_mat$cell.to.index.map
  
  print("Registered: ")
  print(apply(registration_mat, 2, function(r) {sum(r>0)}))
  print("Detected: ")
  print(unlist(lapply(spike_train, function(st){ncol(st[[1]])})))
  
  cells_mat <- matrix(rep(-1, times=len(registration_mat)), nrow=nrow(registration_mat))
  

  for (session_idx in 1:ncol(cells_mat)) { 
  print(sprintf("Running session %d ", session_idx))
  session_spike_train <- t(spike_train[[session_idx]][[1]])
  session_stim_trace <- as.vector(stim_trace[[session_idx]][[1]])

  tmp <- equalize_prior(session_spike_train, session_stim_trace, verbose=T)


  session_spike_train <- tmp$spike_train
  session_stim_trace <- tmp$stim_trace
  rm(tmp)
  
  
  ind <- unlist(lapply(1:nrow(session_spike_train), 
                       function(neur_idx) {which(registration_mat[,session_idx] == neur_idx)}))
  
  spike_train_ind <- registration_mat[ind,session_idx]
  
  
  cells_groups <- rep(DORMENT, len(spike_train_ind))
  
  place_cells_metadata <- get_place_cells(session_spike_train, session_stim_trace, verbose=F)
  active_ind <- place_cells_metadata$ind
  
  
  non_random_non_cyclic_ind <- 
    which(spike_train_ind %in% active_ind[place_cells_metadata$random_pval >= 0.05 & 
                                          place_cells_metadata$cyclic_pval >= 0.05])

  random_non_cyclic_ind <- 
    which(spike_train_ind %in% active_ind[place_cells_metadata$random_pval < 0.05 & 
                                            place_cells_metadata$cyclic_pval >= 0.05])  
  
  non_random_cyclic_ind <- 
    which(spike_train_ind %in% active_ind[place_cells_metadata$random_pval >= 0.05 & 
                                            place_cells_metadata$cyclic_pval < 0.05])    
  
  random_cyclic_ind <- 
    which(spike_train_ind %in% active_ind[place_cells_metadata$random_pval < 0.05 & 
                                            place_cells_metadata$cyclic_pval < 0.05])    
  
  cells_groups[non_random_non_cyclic_ind] <- NON_RANDOM_NON_CYCLIC
  cells_groups[random_non_cyclic_ind] <- RANDOM_NON_CYLIC
  cells_groups[non_random_cyclic_ind] <- NON_RANDOM_CYCLIC
  cells_groups[random_cyclic_ind] <- RANDOM_CYCLIC
  cells_mat[ind,session_idx] <- cells_groups

  }
  
  fpath <- sprintf("%s/%s", path, result_file_name)
  print(sprintf("Done, saving to %s", fpath))
  save(file=fpath, cells_mat)
}


main_place_cells_pooled_across_days <- function (path, result_file_name="pooled_st.Rda", old=F) {
  
  # Groups
  DORMENT = 0;
  NON_RANDOM_NON_CYCLIC = 1;
  RANDOM_NON_CYLIC = 2; # Swithcers
  NON_RANDOM_CYCLIC = 3;
  RANDOM_CYCLIC = 4;
  
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
  
  
  registration_path <- sprintf("%s%s",
                               str_split(path, "Matrices//")[[1]][1], registration_mat_filename)
  
  registration_mat <- readMat(registration_path)
  registration_mat <- registration_mat$cell.to.index.map
  
  print("Registered: ")
  print(apply(registration_mat, 2, function(r) {sum(r>0)}))
  print("Detected: ")
  print(unlist(lapply(spike_train, function(st){ncol(st[[1]])})))
  

  
  final_mat <- c()
  final_stim_trace <- c()
  session_sizes <- c()
  
  for (session_idx in 1:ncol(cells_mat)) { 
    print(sprintf("Running session %d ", session_idx))
    session_spike_train <- t(spike_train[[session_idx]][[1]])
    session_stim_trace <- as.vector(stim_trace[[session_idx]][[1]])
    
    tmp <- equalize_prior(session_spike_train, session_stim_trace, verbose=T)
    
    
    session_spike_train <- tmp$spike_train
    session_stim_trace <- tmp$stim_trace
    rm(tmp)
    
    
    ind <- unlist(lapply(1:nrow(session_spike_train), 
                         function(neur_idx) {which(registration_mat[,session_idx] == neur_idx)}))
    
    spike_train_ind <- registration_mat[ind,session_idx]
    
    
    cells_groups <- rep(DORMENT, len(spike_train_ind))
    
    
    
    session_pooled_spike_train <- matrix(rep(0, times=nrow(registration_mat) * ncol(session_spike_train)),
                                         nrow=nrow(registration_mat))
     
    session_pooled_spike_train[ind,] <- session_spike_train[spike_train_ind,]
    
    final_mat <- cbind(final_mat,
                       session_pooled_spike_train)
    
    session_sizes <- c(session_sizes, ncol(session_spike_train))
    
    final_stim_trace <- c(final_stim_trace, session_stim_trace)

  }
  
  env_a_pooled_st <- final_mat[,1:sum(session_sizes[1:8])]
  env_b_pooled_st <- final_mat[,(sum(session_sizes[1:8]) + 1):sum(session_sizes)]
  env_a_pooled_stim_trace <- final_stim_trace[1:sum(session_sizes[1:8])]
  env_b_pooled_stim_trace <- final_stim_trace[(sum(session_sizes[1:8]) + 1):sum(session_sizes)]
  
  place_cells_metadata_env_A <- get_place_cells(env_a_pooled_st, env_a_pooled_stim_trace, active_bins_threshold = 20, cyclic_only=T, verbose=T)
  place_cells_metadata_env_B <- get_place_cells(env_b_pooled_st, env_b_pooled_stim_trace, verbose=T)
  # active_ind <- place_cells_metadata$ind
  # 
  # fpath <- sprintf("%s/%s", path, result_file_name)
  # print(sprintf("Done, saving to %s", fpath))
  # save(file=fpath, cells_mat)
}
  # 
  # first_session_spike_train <- t(spike_train[[1]][[1]])
  # first_session_stim_trace <- as.vector(stim_trace[[1]][[1]])
  # 
  # tmp <- equalize_prior(first_session_spike_train, first_session_stim_trace, verbose=T)
  # 
  # 
  # first_session_spike_train <- tmp$spike_train
  # first_session_stim_trace <- tmp$stim_trace
  # rm(tmp)
  # 
  # first_session_place_cells <- get_place_cells(first_session_spike_train,
  #                                              first_session_stim_trace,
  #                                              verbose=T)
  # 
  # first_active <- first_session_place_cells$ind
  # first_switchers <- first_session_place_cells$ind[first_session_place_cells$switchers]
  # first_cyclic_non_pcs <- first_session_place_cells$ind[first_session_place_cells$cyclic_non_pcs]
  # first_random_non_pcs <- first_session_place_cells$ind[first_session_place_cells$random_non_pcs]
  # first_cyclic_pcs <- first_session_place_cells$ind[first_session_place_cells$cyclic_pcs]
  # first_random_pcs <- first_session_place_cells$ind[first_session_place_cells$random_pcs]
  # 
  # cell_types_list <- list(first_active,
  #                         first_switchers,
  #                         first_cyclic_non_pcs,
  #                         first_random_non_pcs,
  #                         first_cyclic_pcs,
  #                         first_random_pcs)
  # 
  # 
  # for (idx in 2:len(spike_train)) {
  #   session_spike_train <- t(spike_train[[idx]][[1]])
  #   session_stim_trace <- as.vector(stim_trace[[idx]][[1]])
  #   
  #   tmp <- equalize_prior(session_spike_train, session_stim_trace, verbose=T)
  #   session_spike_train <- tmp$spike_train
  #   session_stim_trace <- tmp$stim_trace
  #   rm(tmp)
  #   
  #   session_place_cells <- get_place_cells(session_spike_train,
  #                                          session_stim_trace,
  #                                          verbose=T)
  #   
  #   session_active <- session_place_cells$ind
  #   session_switchers <- session_place_cells$ind[session_place_cells$switchers]
  #   session_cyclic_non_pcs <- session_place_cells$ind[session_place_cells$cyclic_non_pcs]
  #   session_random_non_pcs <- session_place_cells$ind[session_place_cells$random_non_pcs]
  #   session_cyclic_pcs <- session_place_cells$ind[session_place_cells$cyclic_pcs]
  #   session_random_pcs <- session_place_cells$ind[session_place_cells$random_pcs]
  #   
  #   active_cells_on_first_detected_on_second <- registration_mat[which(registration_mat[,1] %in% first_active),idx]
  #   active_cells_on_first_detected_on_second <- sort(active_cells_on_both_session[active_cells_on_both_session != 0])
  #   active_on_both <- session_active[which(session_active %in% active_cells_on_first_detected_on_second)]
  #   
  #   random_non_pc_on_first_detected_on_second <- registration_mat[which(registration_mat[,1] %in% first_random_non_pcs),idx]
  #   random_non_pc_on_first_detected_on_second <- sort(random_non_pc_on_first_detected_on_second[random_non_pc_on_first_detected_on_second != 0])
  #   random_non_pc_on_first_active_on_second <- session_active[which(session_active %in% random_non_pc_on_first_detected_on_second)]
  #   random_non_pc_on_first_switcher_on_second <- session_switchers[which(session_switchers %in% random_non_pc_on_first_detected_on_second)]
  #   random_non_pc_on_first_cyclic_pc_on_second <- session_cyclic_pcs[which(session_cyclic_pcs %in% random_non_pc_on_first_detected_on_second)]
  #   random_non_pc_on_first_random_non_pc_on_second <- session_random_non_pcs[which(session_random_non_pcs %in% random_non_pc_on_first_detected_on_second)]
  #   
  #   
  #   switcher_on_first_detected_on_second <- registration_mat[which(registration_mat[,1] %in% first_switchers),idx]
  #   switcher_on_first_detected_on_second <- sort(switcher_on_first_detected_on_second[switcher_on_first_detected_on_second != 0])
  #   switcher_on_first_active_on_second <- session_active[which(session_active %in% switcher_on_first_detected_on_second)]
  #   switcher_on_first_switcher_on_second <- session_switchers[which(session_switchers %in% switcher_on_first_detected_on_second)]
  #   switcher_on_first_cyclic_pc_on_second <- session_cyclic_pcs[which(session_cyclic_pcs %in% switcher_on_first_detected_on_second)]
  #   switcher_on_first_random_non_pc_on_second <- session_random_non_pcs[which(session_random_non_pcs %in% switcher_on_first_detected_on_second)]
  #   
  #   
  #   cyclic_pc_on_first_detected_on_second <- registration_mat[which(registration_mat[,1] %in% first_cyclic_pcs),idx]
  #   cyclic_pc_on_first_detected_on_second <- sort(cyclic_pc_on_first_detected_on_second[cyclic_pc_on_first_detected_on_second != 0])
  #   cyclic_pc_on_first_active_on_second <- session_active[which(session_active %in% cyclic_pc_on_first_detected_on_second)]
  #   cyclic_pc_on_first_switcher_on_second <- session_switchers[which(session_switchers %in% cyclic_pc_on_first_detected_on_second)]
  #   cyclic_pc_on_first_cyclic_pc_on_second <- session_cyclic_pcs[which(session_cyclic_pcs %in% cyclic_pc_on_first_detected_on_second)]
  #   cyclic_pc_on_first_random_non_pc_on_second <- session_random_non_pcs[which(session_random_non_pcs %in% cyclic_pc_on_first_detected_on_second)]
  #   
  #   
  # }
#   
# }





