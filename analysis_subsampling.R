library(pheatmap)
generate_sampled_ind <- function(sample_duration, 
                                 total_frames,
                                 num_of_subsample_repetitions=500) {
  nframes <- sample_duration / dt 
  ind <- lapply(c(1:num_of_subsample_repetitions), function(i) {sort(sample(1:total_frames, nframes))})
  
  return(do.call(rbind,ind))
}

calculate_subsample_place_cell_fraction <- function(ind, spike_train, stim_trace, return_pval=T){
  
  sliced_spike_train <- spike_train[,ind]
  sliced_stim_trace <- stim_trace[ind]
  
  processed <- preprocess_spike_train(sliced_spike_train, 
                                      sliced_stim_trace)
  
  if (len(processed$ind) == 0) { return(licomst(empty=T))}
  
  working_cells_spike_train <- processed$working_cells_spike_train
  working_firing_rate <- processed$working_firing_rate
  working_time_bins_per_cells <- processed$working_time_bins_per_cells
  
  if (len(processed$ind) == 1) {
    working_cells_spike_train <- t(as.matrix(working_cells_spike_train))
  }
  
  tmp <- compute_tuning(working_cells_spike_train, 
                        sliced_stim_trace)
  stim_prob <- tmp[[1]]
  working_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  sliced_SI <- compute_SI(stim_prob, working_tuning_curve, working_firing_rate)
  
  
  random_SI <- lapply(c(1:nrow(working_cells_spike_train)), 
                      
                      function(i, var2) { 
                        compute_place_signif(i,
                                             working_cells_spike_train,
                                             sliced_stim_trace,
                                             working_firing_rate,
                                             sliced_SI,
                                             shuffle_type=var2,
                                             num_shuffles = 500,
                                             return_pval = return_pval,
                                             rate_to_use=2)},
                      var2="r")
  
  cyclic_SI <- lapply(c(1:nrow(working_cells_spike_train)), 
                      
                      function(i, var2) { 
                        compute_place_signif(i,
                                             working_cells_spike_train,
                                             sliced_stim_trace,
                                             working_firing_rate,
                                             sliced_SI,
                                             shuffle_type=var2,
                                             num_shuffles = 500,
                                             return_pval = return_pval,
                                             rate_to_use=2)},
                      var2="c")
  
  
  if (return_pval) {
    random_SI = unlist(random_SI)
    cyclic_SI = unlist(cyclic_SI)
  } else {
    names(random_SI) <- processed$ind
    names(cyclic_SI) <- processed$ind
  }
  res = list(empty=F,
             SI=sliced_SI,
             random=random_SI,
             cyclic=cyclic_SI,
             ind=processed$ind)
  
  return(res)
  
}

main_calc_place_fraction_by_subsamples <- function(spike_train, 
                                                   stim_trace, 
                                                   path, 
                                                   ext="", 
                                                   dur_bin=25, 
                                                   number_of_duration_samples=8,
                                                   return_df=F, 
                                                   save_res=T,
                                                   plot_result=T,
                                                   result_file_name="place_cells_by_sample_duration.Rda",
                                                   subsamples_folder="subsampled_matrices") {
  
  prcsd <- preprocess_spike_train(spike_train, stim_trace)
  
  pct_by_duration <- c()
  dur <- floor(seq(10,ncol(spike_train) * dt, length.out=number_of_duration_samples))

  n_of_reps <- floor(1 / ((dur) / max(dur)) * 3)
  
  for (i in 1:len(dur)) {
    duration <- dur[i]
    ind_mat <- generate_sampled_ind(duration, 
                                    ncol(spike_train),
                                    num_of_subsample_repetitions = n_of_reps[i])
    random_mat <- c()
    cyclic_mat <- c()
    
    for (i in 1:nrow(ind_mat)) {
      
      res = calculate_subsample_place_cell_fraction(ind_mat[i,],
                                                    spike_train,
                                                    stim_trace)
      
      random_SI = rep(NA, nrow(spike_train))
      cyclic_SI = rep(NA, nrow(spike_train))
      
      if (!res$empty) { 
        random_SI[res$ind] <- res$random
        cyclic_SI[res$ind] <- res$cyclic
      }
      
      random_mat <- cbind(random_mat, random_SI)
      cyclic_mat <- cbind(cyclic_mat, cyclic_SI)
      
    }
    
    print(path)
    
    if (save_res) {
      print("saving!")
      
      save(cyclic_mat,
           file=sprintf("%s/%s/cyclic_mat_%d.Rda",
                        path,
                        subsamples_folder,
                        duration))
      
      save(random_mat,
           file=sprintf("%s/%s/rand_mat_%d.Rda",
                        path,
                        subsamples_folder,
                        duration))      
      
    }
    
    relevant_random_cells <- apply(random_mat, 1, function(r) {sum(!is.na(r))}) > 0
    
    random_pct <- NA
    if(sum(relevant_random_cells) > 0) {
      print(sum(relevant_random_cells))
      random_mat_t <- random_mat[relevant_random_cells,]    
      
      if (sum(relevant_random_cells == 1)) {
        random_mat_t <- t(as.matrix(random_mat_t))
      }
      
      random_mat_t <- apply(random_mat_t, 2, function(r) {as.numeric(r < 0.05)})
      
      if (sum(relevant_random_cells == 1)) {
        random_mat_t <- t(as.matrix(random_mat_t))
      }
      
      random_pct <- apply(random_mat_t, 2, function(col) {sum(col, na.rm = T) / len(col)})
      random_pct <- mean(random_pct, na.rm=T)
    }
    
    cyclic_pct <- NA
    
    relevant_cyclic_cells <- apply(cyclic_mat, 1, function(r) {sum(!is.na(r))}) > 0
    if(sum(relevant_cyclic_cells) > 0) {
      
      cyclic_mat_t <- cyclic_mat[relevant_cyclic_cells,]    
      
      if (sum(relevant_cyclic_cells == 1)) {
        cyclic_mat_t <- t(as.matrix(cyclic_mat_t))
      }
      
      cyclic_mat_t <- apply(cyclic_mat_t, 2, function(r) {as.numeric(r < 0.05)})
      
      if (sum(relevant_cyclic_cells == 1)) {
        cyclic_mat_t <- t(as.matrix(cyclic_mat_t))
      }      
      
      cyclic_pct <- apply(cyclic_mat_t, 2, function(col) {sum(col, na.rm = T) / len(col)})
      cyclic_pct <- mean(cyclic_pct)
    }
    
    
    
    pct_by_duration <- rbind(pct_by_duration, c(random_pct, cyclic_pct))
    
  }
  
  tmp <- as.data.frame(pct_by_duration)
  rownames(tmp) <- dur
  tmp$total_cells <- rep(len(prcsd$ind), times=nrow(tmp))
  tmp$total_frames <- rep(ncol(spike_train), times=nrow(tmp))
  
  if (save_res) {
    save(tmp,
         file=sprintf("%s/%s%s",
                      path,                    
                      ext,
                      result_file_name))  
  }
  
  if (plot_result) {
    png(sprintf("%s/pct_by_subsample_%s.png", path, ext))
    
    
    plot(dur, tmp[,1], type = "o", pch = 19, cex=1, 
         col = "black", 
         xlab = "Sample duration (s)", 
         ylab = "Percentage of place cells (%)", 
         main=str_split(path, "set")[[1]][2])
    
    lines(dur, tmp[,2], type = "o", pch = 19, cex=1,
          col = "red")
    
    legend("bottomright", legend=c("Random shuffles", "Cyclic shuffles"),
           col=c("Black", "Red"), lty=1:1, cex=0.8, pch=19:19)
    
    dev.off()
  }
  
  if (return_df) {
    return(pct_by_duration)
  }
}


main_calc_place_fraction_by_activity <- function(spike_train, stim_trace, path, result_file_name="tuning_df.Rda") {
  
  processed <- preprocess_spike_train(spike_train, stim_trace)
  
  working_cells_spike_train <- processed$working_cells_spike_train
  working_firing_rate <- processed$working_firing_rate
  working_time_bins_per_cells <- processed$working_time_bins_per_cells
  
  tmp <- compute_tuning(working_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  working_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  working_SI <- compute_SI(stim_prob, 
                           working_tuning_curve, 
                           working_firing_rate)
  
  
  # Compute place tuning cells for cyclic shuffles
  place_tuning_cyclic <- sapply(c(1:nrow(working_cells_spike_train)), 
                                function(i, var2) 
                                {compute_place_signif(i,
                                                      working_cells_spike_train,
                                                      stim_trace,
                                                      working_firing_rate,
                                                      working_SI,
                                                      shuffle_type=var2,
                                                      num_shuffles = 500,
                                                      verbose = T)},
                                var2="c")
  
  
  # Cyclic binning and percentage calculation
  bins_cyclic <- bin(working_time_bins_per_cells, bin_size_active_time_bins)
  is_sig_cyclic <- as.numeric(place_tuning_cyclic < tuning_significance_threshold)
  
  place_tuned_by_time_bins_cyclic_df  <-
    data.frame(pval=place_tuning_cyclic,
               bins=bins_cyclic,
               sig=is_sig_cyclic) 
  
  percentage_cyclic <- c()
  
  for (lv in levels(place_tuned_by_time_bins_cyclic_df$bins)) {
    ind <- place_tuned_by_time_bins_cyclic_df$bins == lv
    subdf <- place_tuned_by_time_bins_cyclic_df[ind,]
    pct <- 0
    if (nrow(subdf) > 0) {
      pct = sum(subdf$sig == 1) / nrow(subdf)
    }
    
    percentage_cyclic <- c(percentage_cyclic, pct)
    
  }
  
  active_cyclic <- sapply(levels(bins_cyclic), function(b) 
  {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
  
  
  ind_to_rem <- which(percentage_cyclic == 0)[-1]
  
  if (len(ind_to_rem) > 0) {
    percentage_cyclic <- percentage_cyclic[-ind_to_rem]
    active_cyclic <- active_cyclic[-ind_to_rem]
  }
  
  # Compute place tuning cells for random shuffles
  place_tuning <- sapply(c(1:nrow(working_cells_spike_train)), 
                         function(i, var2) 
                         {compute_place_signif(i,
                                               working_cells_spike_train,
                                               stim_trace,
                                               working_firing_rate,
                                               working_SI,
                                               shuffle_type=var2,
                                               num_shuffles = 500,
                                               verbose = T)},
                         var2="r")
  
  bins <- bin(working_time_bins_per_cells, bin_size_active_time_bins)
  is_sig <- as.numeric(place_tuning < tuning_significance_threshold)
  
  place_tuned_by_time_bins_df  <- data.frame(pval=place_tuning,
                                             bins=bins,
                                             sig=is_sig) 
  percentage <- c()
  for (lv in levels(place_tuned_by_time_bins_df$bins)) {
    ind <- place_tuned_by_time_bins_df$bins == lv
    subdf <- place_tuned_by_time_bins_df[ind,]
    pct <- 0
    if (nrow(subdf) > 0) {
      pct = sum(subdf$sig == 1) / nrow(subdf)
    }
    percentage <- c(percentage, pct)
    
  }
  
  active <- sapply(levels(bins), function(b) 
  {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
  
  
  ind_to_rem <- which(percentage == 0)[-1]
  if (len(ind_to_rem) > 0) {
    percentage <- percentage[-ind_to_rem]
    active <- active[-ind_to_rem]
  }
  
  
  if (plot_res) { 
    
    if (save_res) {
      png(sprintf("%s\\%s", path, "pct_by_time_bins_single_session.png"))
    }
    
    plot(active, percentage, type = "o", frame = FALSE, pch = 19, cex=1, 
         col = "red", 
         xlab = "Active time bins", 
         ylab = "Percentage of place cells", 
         main=str_split(path, "set")[[1]][2])
    
    lines(active_cyclic, percentage_cyclic, type = "o", frame = FALSE, pch = 19, cex=1,
          col = "black")
    legend("bottomright", legend=c("Random shuffles", "Cyclic shuffles"),
           col=c("red", "black"), lty=1:1, cex=0.8, pch=19:19)
    
    if (save_res) {
      dev.off()
    }
  }
  
  
  tuning_df <- data.frame(active_bins=working_time_bins_per_cells,
                          random_significance=place_tuning,
                          cyclic_significance=place_tuning_cyclic)
  
  random_df <- data.frame(bins=active,
                          percentage=percentage)
  
  rownames(random_df) = c()
  
  cyclic_df <- data.frame(bins=active_cyclic,
                          percentage=percentage_cyclic)
  
  rownames(cyclic_df) = c()
  
  save(tuning_df, file=sprintf("%s\\%s",
                               path,
                               result_file_name))
  
  
  save(random_df, file=sprintf("%s\\%s",
                               path,
                               "random_df.Rda"))
  
  save(cyclic_df, file=sprintf("%s\\%s",
                               path,
                               "cyclic_df.Rda"))
}

calculate_pval_by_ind <- function(ind, spike_train, stim_trace, return_pval=T, cyclic_only=T){
  
  sliced_spike_train <- spike_train[,ind]
  sliced_stim_trace <- stim_trace[ind]
  
  processed <- preprocess_spike_train_no_threshold(sliced_spike_train, sliced_stim_trace)
  
  if (len(processed$ind) == 0) { return(list(empty=T))}
  
  working_cells_spike_train <- processed$working_cells_spike_train
  working_firing_rate <- processed$working_firing_rate
  working_time_bins_per_cells <- processed$working_time_bins_per_cells
  
  if (len(processed$ind) == 1) {
    working_cells_spike_train <- t(as.matrix(working_cells_spike_train))
  }
  
  tmp <- compute_tuning(working_cells_spike_train, 
                        sliced_stim_trace)
  stim_prob <- tmp[[1]]
  working_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  
  sliced_SI <- compute_SI(stim_prob, working_tuning_curve, working_firing_rate)
  
  
  if (!cyclic_only) {
    random_SI <- lapply(c(1:nrow(working_cells_spike_train)), 
                        
                        function(i, var2) { 
                          compute_place_signif(i,
                                               working_cells_spike_train,
                                               sliced_stim_trace,
                                               working_firing_rate,
                                               sliced_SI,
                                               shuffle_type=var2,
                                               num_shuffles = 500,
                                               return_pval = return_pval,
                                               rate_to_use=2)},
                        var2="r")
    
  } else {
    random_SI <- rep(NA, times=nrow(working_cells_spike_train))
  }
  
  cyclic_SI <- lapply(c(1:nrow(working_cells_spike_train)), 
                      
                      function(i, var2) { 
                        compute_place_signif(i,
                                             working_cells_spike_train,
                                             sliced_stim_trace,
                                             working_firing_rate,
                                             sliced_SI,
                                             shuffle_type=var2,
                                             num_shuffles = 500,
                                             return_pval = return_pval,
                                             rate_to_use=2, 
                                             verbose=F)},
                      var2="c")
  
  
  if (return_pval) {
    random_SI = unlist(random_SI)
    cyclic_SI = unlist(cyclic_SI)
  } else {
    names(random_SI) <- processed$ind
    names(cyclic_SI) <- processed$ind
  }
  res = list(empty=F,
             SI=sliced_SI,
             random=random_SI,
             cyclic=cyclic_SI,
             ind=processed$ind)
  
  return(res)
}

main_calc_pvalue_by_subsamples <- function(spike_train, 
                                           stim_trace, 
                                           path, 
                                           ext="", 
                                           return_df=F,
                                           save_res=T, 
                                           plot_res=T, 
                                           result_file_name="cyclic_pval_by_dur.Rda",
                                           number_of_duration_samples=6) {
  
  
  pct_by_duration <- c()
  dur <- seq(10,ncol(spike_train) * dt, length.out=number_of_duration_samples)
  n_of_reps <- c(rep(10, times=2), rep(5, times=2), rep(3, times=2))
  
  random_list <- list()
  cyclic_list <- list()
  
  for (i in 1:len(dur)) {
    
    duration <- dur[i]
    ind_mat <- generate_sampled_ind(duration, 
                                    ncol(spike_train),
                                    num_of_subsample_repetitions = n_of_reps[i])
    random_mat <- c()
    cyclic_mat <- c()
    
    for (i in 1:nrow(ind_mat)) {
      print(sprintf("Duration (%.1f), rep (%d)", duration, i))
      res = calculate_pval_by_ind(ind_mat[i,], spike_train, stim_trace)
      
      random_SI = rep(NA, nrow(spike_train))
      cyclic_SI = rep(NA, nrow(spike_train))
      
      random_mat <- cbind(random_mat, res$random)
      cyclic_mat <- cbind(cyclic_mat, res$cyclic)
      
    }
    
    
    random_list <- append(random_list, list(random_mat))
    cyclic_list <- append(cyclic_list, list(cyclic_mat))
    
    
    
  }
  
  cyclic_pval_mat <- do.call(cbind, lapply(cyclic_list, rowMeans))
  
  cyclic_list$spike_train_size <- ncol(spike_train)
  save(cyclic_list,
       file=sprintf("%s\\%s%s",                    
                    path,
                    result_file_name,
                    ext))  
  
  
  if (plot_res) { 
    png(sprintf("%s\\cyclic_pval_by_dur%s.png", path, ext))
    mtmp <- apply(cyclic_pval_mat, 2, 
                  function(c) {return((c + 10^-30) / (cyclic_pval_mat[,ncol(cyclic_pval_mat)] + 10^-30))})
    
    mtmp <- colMeans(mtmp)
    mtmp <- mtmp / sum(mtmp)
    
    plot(dur, mtmp, type = "o", pch = 19, cex=1, 
         col = "black", 
         xlab = "Sample duration (s)", 
         ylab = "Normalized delta Pval (AU)", 
         main=str_split(path, "set")[[1]][2])
    
    
    dev.off()
  }
  
  if (return_df) {
    return(pct_by_duration)
  }
}


main_calc_tuning_cor_by_subsamples <- function(spike_train, 
                                               stim_trace, 
                                               path, 
                                               ext="", 
                                               save_res=T, 
                                               plot_res=F, 
                                               result_file_name="properties.Rda",
                                               number_of_duration_samples=c(3:6)) {
  
  
  # Calculate whose a place cell and who isn't
  place_cells_metadata <- get_place_cells(spike_train,
                                          stim_trace,
                                          verbose=T)
  
  tuning_properties_all <- get_tuning_properties(spike_train, stim_trace)
  
  active_ind <- place_cells_metadata$ind
  
  groups_indices <- 
    list(all=1:nrow(spike_train),
         active=active_ind,
         switchers=active_ind[place_cells_metadata$switchers],
         cyclic_non_pcs=active_ind[place_cells_metadata$cyclic_non_pcs],
         random_non_pcs=active_ind[place_cells_metadata$random_non_pcs],
         cyclic_pcs=active_ind[place_cells_metadata$cyclic_pcs],
         random_pcs=active_ind[place_cells_metadata$random_pcs])
  
  tuning_properties_by_group <- 
    lapply(groups_indices, 
           function(ind) {
             tp <-  tuning_properties_all$tuning_properties[ind,]
             tp <- cbind(tp, ind)
             
             return(tp)
           })
  
  across_subsamples <- list()
  
  # Divide data into consecutive parts (2, 3, 4 ....)
  for (nd in number_of_duration_samples) {  
    indices <- floor(seq(0, ncol(spike_train), length.out=nd))
    
    # Calculate tuning properties for each sub sample
    tuning_list <- 
      lapply(1:(len(indices) - 1), 
             function(idx) {
               ind <- (indices[idx] + 1):(indices[idx + 1])
               tuning_chunk <- get_tuning_properties(as.matrix(spike_train[,ind]), stim_trace[ind])
               return(tuning_chunk)})
    
    
    # Calculate tuning correlation of subsamples
    tuning_corr_vec <- c()
    
    for (s_idx in 1:nrow(spike_train)) { 
      neur_tuning <- c()
      
      for (i in 1:(len(indices) - 1)){
        neur_tuning <- rbind(neur_tuning,
                             tuning_list[[i]]$tuning_curves[s_idx,])  
      }
      
      
      print("#####")
      print(dim(neur_tuning))
      print(neur_tuning)
      
      if (all(neur_tuning == 0)) {
        tuning_corr_vec <- c(tuning_corr_vec, NA)
        print("ALL are zero!!!")
      } else {
        cmt  <- cor(t(neur_tuning), t(neur_tuning))
        
        print(dim(cmt))
        print(cmt)
        print(1:nrow(cmt))
        print("#>>>>>>>>>>>>>>>")
        for(i in 1:nrow(cmt)) {cmt[i,i] <- NA}
        tuning_corr_vec <- c(tuning_corr_vec, mean(c(cmt), na.rm=T))
      }
      
    }
    
    # Calcualte tuning properties (Number peaks, SI ...) for each subsample
    tuning_properties_mats <- as.list(rep(0, times=ncol(tuning_list[[1]]$tuning_properties)))
    names(tuning_properties_mats) <- colnames(tuning_list[[1]]$tuning_properties)
    
    for (property in names(tuning_properties_mats)) {
      property_mat <- c()
      for (i in 1:(len(indices) - 1)) {
        property_mat <- cbind(property_mat,
                              tuning_list[[i]]$tuning_properties[,property])
      }
      
      tuning_properties_mats[[property]]  <- property_mat
    }
    
    # For each group of cells (place cells, all cells, active cells... ) slice tuning correlationm, and tuning properties (with their correlations)
    tuning_properties_by_groups_across_subsamples <- 
      lapply(groups_indices,
             function(cells_ind) {
               
               property_mats <- 
                 lapply(tuning_properties_mats,
                        function(property) { return(property[cells_ind,])})
               
               
               property_corrs <- 
                 lapply(property_mats,
                        function(mat) {
                          if (len(cells_ind == 1)) {
                            print("JUST ONE CELL!!!!")
                            return(NA)
                          }
                          cmt <- cor(mat,mat)
                          print("######@@@@@####")
                          print(mat)
                          print(dim(mat))
                          print(cmt)
                          print(1:nrow(cmt))
                          print(len(cells_ind))
                          print("######@@@@@####")
                          for(i in 1:nrow(cmt)) {cmt[i,i] <- NA}
                          
                          return(mean(cmt, na.rm=T))
                        })
               
               tuning_corr_vec
               
               return(list(property_mats=property_mats,
                           property_corrs=property_corrs,
                           tuning_corr=tuning_corr_vec[cells_ind]))
               
             })
    
    
    melted_df <- data.frame()
    
    for(group in names(tuning_properties_by_groups_across_subsamples)) {
      for(property in names(tuning_properties_by_groups_across_subsamples[[group]]$property_corrs)) {
        melted_df <- rbind(melted_df,
                           c(group,
                             property,
                             tuning_properties_by_groups_across_subsamples[[group]]$property_corrs[[property]])) 
      }
    }
    
    colnames(melted_df) <- c("Group", "Property", "Corr")
    melted_df$Corr <- as.numeric(melted_df$Corr)
    tuning_properties_by_groups_across_subsamples$melted_df <- melted_df  
    
    subsample_name <- sprintf("subsamples_%d", (nd - 1))
    across_subsamples[[subsample_name]] <- tuning_properties_by_groups_across_subsamples
  }
  
  
  result <- list(place_cells_metadata = place_cells_metadata,
                 groups_indices = groups_indices,
                 tuning_properties_by_group=tuning_properties_by_group,
                 across_subsamples=across_subsamples)    
  
  fname <- sprintf("%s\\%s", path, result_file_name)
  print(sprintf("Done! Saving tuning properties to %s", fname))
  save(file=fname, result)
}

main_overlaid_simulation_real_by_sample_size <- function(spike_train, 
                                                         stim_trace, 
                                                         path, 
                                                         pval_threshold=0.05, 
                                                         sim_reps=10,
                                                         n_bins=6) {
  params_df <- extract_simulation_params(path)
  
  fr <- rowMeans(spike_train) / dt
  processed_real <- preprocess_spike_train(spike_train, stim_trace)
  
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  real_spike_train_df <- main_calc_place_fraction_by_subsamples(path = path, spike_train, stim_trace, 
                                                                ext="real_dur_15", 
                                                                number_of_duration_samples = n_bins, 
                                                                return_df=T, 
                                                                save_res=F,
                                                                plot_res=F)
  
  pcts <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
  final_result <- list()
  for (pct in pcts) {
    params <- params_df[which(pct ==  as.numeric(rownames(params_df))),]
    
    simulated_df <- c()
    for (sim_round in 1:sim_reps) {
      
      print(sprintf("Simulating for %.3f", pct))
      simulated_tuning_curve <-
        generate_tuning_curves_cost(n = nrow(spike_train),
                                    percentage = pct,
                                    average_width = params["average_width"], 
                                    sd_width = params["sd_width"],
                                    fixed_fr=fr,
                                    noise=params["noise"],
                                    double_peak_pct = params["double_peak_pct"],
                                    plot=F)
      
      pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                              true_time_bins_per_cells,
                                              stim_trace,
                                              jump = 0.05,
                                              verbose = T)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)
      
      simulated_pct_by_sample <- 
        main_calc_place_fraction_by_subsamples(path = path, 
                                               generated_spike_train, 
                                               stim_trace, 
                                               ext=sprintf("sim_dur_15_%.3f_r%d", pct, sim_round), 
                                               number_of_duration_samples = n_bins, 
                                               return_df=T, 
                                               save_res=F,
                                               plot_res=F)
      
      par(mfrow=c(1,2))
      plot(real_spike_train_df[,2], col="black", type="l", lwd=2)
      lines(simulated_pct_by_sample[,2], col="red", lty=2, lwd=2)
      
      plot(real_spike_train_df[,1], col="black", type="l", lwd=2)
      lines(simulated_pct_by_sample[,1], col="red", lty=2, lwd=2)
      simulated_df <- cbind(simulated_df, simulated_pct_by_sample)
      
    }
    
    colnames(simulated_df) <- paste(c("Random", "Cyclic"), rep(1:sim_reps, each=2))
    final_result <- append(final_result, list(simulated_df))
  }
  time_vec <- seq(10, ncol(spike_train), length.out=n_bins) / (ncol(spike_train) / n_bins)  * 20
  real_df <- as.data.frame(real_spike_train_df)
  real_df$time <- time_vec
  colnames(real_df) <- c("Rand", "Cyclic", "Time")
  
  
  
  names(final_result) <- pcts
  
  
  
  
  
  
  g_real_cyc <-  ggplot(real_df, aes(x=Time, y=Cyclic)) + geom_line(size=2) + 
    theme_light() +
    ylim(c(0, 1)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlim(0, max(time_vec)) + 
    labs(x="Sample duration (minutes)", y = "Fraction of place cells (%)")
  
  
  g_real_rand <-  ggplot(real_df, aes(x=Time, y=Rand)) + geom_line(size=2) +
    theme_light() +
    ylim(c(0, 1)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlim(0, max(time_vec)) + 
    labs(x="Sample duration (minutes)", y = "Fraction of place cells (%)")
  
  plots_sim_cyc <- list()
  plots_sim_rand <- list()
  
  dir.create(sprintf("%s\\sim_vs_real\\", path))
  
  sim_pct_colors = brewer.pal(6, "RdYlBu")
  names(sim_pct_colors) <- c("0.5","0.6","0.7","0.8","0.9","1")
  # sim_pct_colors <- c("0.5" = "darkolivegreen3", 
  #                     "0.6" = "dodgerblue1",
  #                     "0.7" = "deeppink2",
  #                     "0.8" = "goldenrod2",
  #                     "0.9" = "orangered2",
  #                     "1" = "darkmagenta")
  
  for (sim_pct in names(final_result)) {
    
    simulation_df <- final_result[[sim_pct]]
    g_sim_cyc <- g_real_cyc #+ ggtitle(sprintf("%d%% place cells simulation", as.numeric(sim_pct) * 100))
    g_sim_rand <- g_real_rand#  + ggtitle(sprintf("%d%% place cells simulation", as.numeric(sim_pct) * 100))
    
    sim_pct_colors[sim_pct]
    
    for (i in 1:sim_reps) {
      tmp_df_cyc <- data.frame(Fraction=simulation_df[, (i * 2)],
                               Time=time_vec)
      
      tmp_df_rand <- data.frame(Fraction=simulation_df[, (i * 2) -1],
                                Time=time_vec)
      
      g_sim_cyc <- 
        g_sim_cyc + geom_line(data=tmp_df_cyc,
                              aes(x=Time, y=Fraction),
                              color=sim_pct_colors[sim_pct],
                              size=2,
                              alpha=0.5)
      
      g_sim_rand <- 
        g_sim_rand + geom_line(data=tmp_df_rand,
                               aes(x=Time, y=Fraction),
                               color=sim_pct_colors[sim_pct],
                               size=2,
                               alpha=0.5)
    }
    
    plots_sim_cyc <- append(plots_sim_cyc, list(g_sim_cyc))
    plots_sim_rand <- append(plots_sim_rand, list(g_sim_rand))
    
    pdf(file=sprintf("%s\\sim_vs_real\\cyclic_sim_vs_real_%s.pdf", path,  sim_pct), 
        height=5,
        width=5)
    
    plot(g_sim_cyc)
    
    dev.off()
    
    pdf(file=sprintf("%s\\sim_vs_real\\random_sim_vs_real_%s.pdf", path,  sim_pct), 
        height=5,
        width=5)
    
    plot(g_sim_rand)
    
    dev.off()   
  }
  
  plots_sim_cyc$nrow = 1
  plots_sim_rand$nrow = 1
  
  pdf(file=sprintf("%s\\sim_vs_real\\random_all_1.5.pdf", path), 
      height=1.5,
      width=1.5 * (len(plots_sim_cyc) - 1))
  
  plot(do.call(arrangeGrob, plots_sim_rand))
  
  dev.off()
  
  pdf(file=sprintf("%s\\sim_vs_real\\random_all_1.75.pdf", path), 
      height=1.75,
      width=1.75 * (len(plots_sim_cyc) - 1))
  
  plot(do.call(arrangeGrob, plots_sim_rand))
  
  dev.off()
  
  pdf(file=sprintf("%s\\sim_vs_real\\random_all_2.pdf", path), 
      height=2,
      width=2 * (len(plots_sim_cyc) - 1))
  
  plot(do.call(arrangeGrob, plots_sim_rand))
  
  dev.off()
  
  pdf(file=sprintf("%s\\sim_vs_real\\random_all_3.pdf", path), 
      height=3,
      width=3 * (len(plots_sim_cyc) - 1))
  
  plot(do.call(arrangeGrob, plots_sim_rand))
  
  dev.off()
  
  
  
  pdf(file=sprintf("%s\\sim_vs_real\\cyclic_all_1.5.pdf", path), 
      height=1.5,
      width=1.5 * (len(plots_sim_cyc) - 1))
  
  
  plot(do.call(arrangeGrob, plots_sim_cyc))
  
  dev.off()
  
  pdf(file=sprintf("%s\\sim_vs_real\\cyclic_all_1.75.pdf", path), 
      height=1.75,
      width=1.75 * (len(plots_sim_cyc) - 1))
  
  
  plot(do.call(arrangeGrob, plots_sim_cyc))
  
  dev.off()  
  
  pdf(file=sprintf("%s\\sim_vs_real\\cyclic_all_2.pdf", path), 
      height=2,
      width=2 * (len(plots_sim_cyc) - 1))
  
  
  plot(do.call(arrangeGrob, plots_sim_cyc))
  
  dev.off()
  
  pdf(file=sprintf("%s\\sim_vs_real\\cyclic_all_3.pdf", path), 
      height=3,
      width=3 * (len(plots_sim_cyc) - 1))
  
  
  plot(do.call(arrangeGrob, plots_sim_cyc))
  
  dev.off()
  
  final_result_2 <- final_result
  final_result_2$real <- real_df
  
  save(file=sprintf("%s\\sim_vs_real\\%s", path, "sim_vs_real_by_sample.R"), final_result_2)
  
}




main_calc_sliding_tuning_similarity <- function(spike_train, 
                                                stim_trace, 
                                                path, 
                                                ext="", 
                                                save_res=T, 
                                                plot_res=F, 
                                                result_file_name="partial_peaks_sliding_tuning_3.Rda",
                                                verbose=F,
                                                use_fractions=F) {
  
  
  # Calculate whose a place cell and who isn't
  place_cells_metadata <- get_place_cells(spike_train,
                                          stim_trace,
                                          verbose=verbose)
  
  
  if (use_fractions) {
    window_sizes=seq(0.1, 1, by=0.1)
    jumps <- seq(0.05,.7,by=0.05)
    
  } else {
    window_sizes=c(200, 215, 225, 235, 250, 275, 300, 350, 400, 450, 500, 600, 650, 750,800,900,1000,1100, 1200, 1300)
    jumps <- c(50)#,100,200,300)
  }
  
  similarity_results <- list()
  
  # <- <- <- <- metrics <- c("cor", "peaks", "eucl", "cosine", "wasserstein")
  metrics <- c("peaks")
  
  for (window_frac in window_sizes) {
    for (jump_frac in jumps) {
      
      if (use_fractions) {
        window = round(ncol(spike_train) * window_frac)
        jump = round(ncol(spike_train) * jump_frac) 
      } else {
        window = window_frac
        jump = jump_frac
      }
      initial_centroid_index = window
      final_centroid_index = (ncol(spike_train) - 1 - window)
      
      if (final_centroid_index < initial_centroid_index) {
        next
      }
      
      if (use_fractions) {
        print(sprintf("fracs - (%f:%f) %d - %d - %d", window_frac, jump_frac, initial_centroid_index,
                      final_centroid_index,
                      jump))
      } else {
        print(sprintf("(%d:%d) %d - %d - %d", window, jump, initial_centroid_index,
                      final_centroid_index,
                      jump))        
      }
      
      centroid_indices <- round(seq(initial_centroid_index, final_centroid_index, by=jump))
      similarity_list <- list() 
      
      for (similarity_metric in metrics) {
        similarity_list[[similarity_metric]]  <- c()
      }
      
      for (cent_idx in centroid_indices) {
        
        first_half <- c((cent_idx - window + 1):cent_idx)
        second_half <- c(cent_idx + 1):(cent_idx + window)
        
        if (verbose) {
          print(sprintf("Comparing (%d:%d) to (%d:%d)",
                        min(first_half), max(first_half), min(second_half), max(second_half)))
        }
        first_tuning_curve <- compute_tuning(spike_train[,first_half], stim_trace[first_half])[[2]]
        second_tuning_curve <- compute_tuning(spike_train[,second_half], stim_trace[second_half])[[2]]
        
        
        for (similarity_metric in metrics) {
          similarity_list[[similarity_metric]] <- 
            cbind(similarity_list[[similarity_metric]],
                  tuning_similarity(first_tuning_curve,
                                    second_tuning_curve,
                                    similarity_metric))
        }
      }        
      
      sliding_window <- sprintf("W%f_J%f", window_frac, jump_frac)
      
      print(sprintf("Done Calculating for %s, (%d metrics)", sliding_window,
                    len(similarity_list)))
      print(dim(similarity_list[[1]]))
      
      similarity_results[[sliding_window]] <- similarity_list
      
    }
  }
  
  final_result <- list()
  final_result$similarity <- similarity_results
  final_result$place_cells_metadata <- place_cells_metadata
  
  fname <- sprintf("%s\\%s", path, result_file_name)
  print(sprintf("Done! Saving sliding tuning properties to %s", fname))
  save(file=fname, final_result)
}


main_calc_prob_by_position_firing <- function(spike_train, 
                                              stim_trace, 
                                              path, 
                                              ext="", 
                                              save_res=T, 
                                              plot_res=T, 
                                              result_file_name="new_probability_to_fire_by_pos.Rda",
                                              verbose=F) {
  
  
  # Calculate whose a place cell and who isn't
  place_cells_metadata <- get_place_cells(spike_train,
                                          stim_trace,
                                          verbose=verbose, 
                                          cyclic_only = T)
  
  
  probability_belt_range <- c(-1,2,3,6,10,15,20,25)
  
  firing_pos <- 
    apply(spike_train,
          1,
          function(c)
          {firing_ind <- which(c > 0);
          stim_trace[firing_ind]})
  
  
  position_difference <- 
    lapply(firing_pos, function(pos) {abs(diff(pos))})
  
  position_diff_prob_mat <- 
    lapply(position_difference, function(diffr)
    {
      cont_vec <- rep(0, times=len(probability_belt_range) - 1)
      
      tabulated <- 
        table(as.numeric(cut(diffr, breaks = probability_belt_range)))
      
      cont_vec[as.numeric(names(tabulated))] <- tabulated
      return(cont_vec)
    })
  
  
  traversals <- get_traversales(stim_trace, nbins=2)
  
  
  traversals_firing <- list()
  traversals_tuning_curve <- c()
  trvs <- seq(1, len(traversals), by=2)
  trvs[len(trvs)] <- len(traversals)
  
  for (trv_idx in 1:(len(trvs)-1)) {
    
    if (trv_idx == 1){
      index_i = 1
    } else {
      index_i <- traversals[trvs[trv_idx]]
    }
    
    
    
    if(trv_idx == len(traversals)) {
      index_f = len(stim_trace)
    } else {
      index_f = traversals[trvs[trv_idx]  + 1]
    }
    
    
    traversal_tuning <- 
    compute_tuning(spike_train = spike_train[,index_i:index_f],
                   stim_trace = stim_trace[index_i:index_f])
    
    
    traversals_tuning_curve <- rbind(traversals_tuning_curve,
                                     apply(traversal_tuning[[2]], 1, function(tc) {if(all(tc==0)){NA} else{which.max(tc)}}))
    
    trav_firing_pos <- 
      apply(spike_train,
            1,
            function(c)
            {
              firing_ind <- which(c[index_i:index_f] > 0);
              stim_trace[firing_ind]
              tabulated <- table(stim_trace[firing_ind])
              if (len(tabulated) == 0) {
                return(NA)
              }
              
              return(as.numeric(names(which.max(tabulated))))
            })
    
    traversals_firing <- append(traversals_firing,
                                list(unlist(trav_firing_pos)))
    
    
  }
  
  trav_firing_mat <- do.call(cbind, traversals_firing)
  traversals_prob_mat <- 
    apply(trav_firing_mat,
          1,
          function(cell) {
            #runs_firing_ind <- which(cell != -1)
            firing_pos <- cell
            #firing_pos <- cell[runs_firing_ind[which(diff(runs_firing_ind) == 1)]]
            diffr <- diff(firing_pos)
            print(diffr)
            diffr <- abs(diffr[!is.na(diffr)])
            
            cont_vec <- rep(0, times=len(probability_belt_range) - 1)
            
            tabulated <- 
              table(as.numeric(cut(diffr, breaks = probability_belt_range)))
            
            cont_vec[as.numeric(names(tabulated))] <- tabulated
            return(cont_vec)
          })
  
  traversals_position_prob_mat <- t(traversals_prob_mat)
  
  position_prob_mat <- do.call(rbind, position_diff_prob_mat)
  
  titles = list(switchers="SW",
                random_non_pcs="R_NPC",
                random_pcs="R_PC",
                cyclic_non_pcs="C_NPC",
                cyclic_pcs="C_PC")
  
  groups = list(switchers=place_cells_metadata$switchers,
                random_non_pcs=place_cells_metadata$random_non_pcs,
                random_pcs=place_cells_metadata$random_pcs,
                cyclic_non_pcs=place_cells_metadata$cyclic_non_pcs,
                cyclic_pcs=place_cells_metadata$cyclic_pcs)
  
  final_result <- list()
  plot_list <- list()
  all_pos <- unique(stim_trace)
  
  for (group_name in names(groups)) {
    
    
    group_ind <- place_cells_metadata$ind[groups[[group_name]]]
    
    if (len(group_ind) <= 1 || min(dim(position_prob_mat[group_ind,])) < 2){
      print(len(group_ind))
      print(dim(position_prob_mat[group_ind,]))
      plot_list <- append(plot_list,
                          list(ggplot()))
      next
    }
    
    pos_matrix <- matrix(rep(0, times=len(all_pos) ** 2),
                         nrow=len(all_pos))
    
    rownames(pos_matrix) <- all_pos
    colnames(pos_matrix) <- all_pos
    
    
    for (fp in firing_pos[group_ind]) {
      if (len(fp) == 1) {
        next
      }
      
      for (idx in 1:(len(fp) - 1)) {
        pos_matrix[fp[idx],
                   fp[idx + 1]] <- 
          pos_matrix[fp[idx],
                     fp[idx + 1]] + 1
      }
    }
    
    ph_prob_pos_mat <- pheatmap(pos_matrix / sum(pos_matrix), cluster_rows=F, cluster_cols=F, border_col=NA, legend=F,
                                main=titles[[group_name]])
    ph_prob_diff_mat <- pheatmap(position_prob_mat[group_ind,]/sum(position_prob_mat[group_ind,]), cluster_rows=F, cluster_cols=F, border_col=NA, legend=F)
    
    
    
    pos_diff_df <- 
      data.frame(prob=colSums(position_prob_mat[group_ind,])/sum(position_prob_mat[group_ind,]),
                 pos_diff=probability_belt_range[-1])
    
    gpos_diff <- 
      ggplot(pos_diff_df) + 
      geom_bar(aes(x=pos_diff, y=prob), stat="identity") + 
      base_plot_theme + scale_y_continuous(expand=c(0,0)) + 
      ylab("") + xlab("Position jump")
    
    
    trav_pos_diff_df <- 
      data.frame(prob=colSums(traversals_position_prob_mat[group_ind,])/sum(traversals_position_prob_mat[group_ind,]),
                 pos_diff=probability_belt_range[-1])
    
    gpos_trav_diff <- 
      ggplot(trav_pos_diff_df) + 
      geom_bar(aes(x=pos_diff, y=prob), stat="identity", fill="royalblue4") + 
      base_plot_theme + scale_y_continuous(expand=c(0,0)) + 
      ylab("") + xlab("Position jump")
    
    gf <- 
      plot_grid(ph_prob_pos_mat[[4]],
                ph_prob_diff_mat[[4]],
                gpos_diff,
                gpos_trav_diff,
                nrow=4,
                rel_heights=c(2,2,1.5, 1.5),
                align="hv")
    
    
    plot_list <- append(plot_list,
                        list(gf))
    
    final_result[[group_name]] <-
      list(pos_matrix = pos_matrix,
           pos_diff = position_prob_mat[group_ind,],
           trav_pos_diff = traversals_position_prob_mat[group_ind,])
  }
  
  plot_list$nrow <- 1
  plot_f_all <- do.call(plot_grid, plot_list)
  
  pos_matrix <- matrix(rep(0, times=len(all_pos) ** 2),
                       nrow=len(all_pos))
  
  rownames(pos_matrix) <- all_pos
  colnames(pos_matrix) <- all_pos
  
  
  for (fp in firing_pos) {
    if (len(fp) == 1) {
      next
    }
    
    for (idx in 1:(len(fp) - 1)) {
      pos_matrix[fp[idx],
                 fp[idx + 1]] <- 
        pos_matrix[fp[idx],
                   fp[idx + 1]] + 1
    }
  }
  
  
  final_result$all_pos_matrix <- pos_matrix
  final_result$all_firing_pos <- firing_pos
  final_result$all_pos_prob_mat <- position_prob_mat
  final_result$place_cells_metadata <- place_cells_metadata
  final_result$trav_firing_mat  <- trav_firing_mat
  final_result$trav_tuning_curve <- t(traversals_tuning_curve)
  final_result$trav_prob <- traversals_position_prob_mat
  
  
  
  if (plot_res) {
    png(sprintf("%s/probability_of_firing_by_pos.png",
                path),
        res = 200,
        units="in",
        height=3 * 1.75,
        width=5 * 1.75)
    plot(plot_f_all)
    dev.off()
  }
  
  
  
  fname <- sprintf("%s\\%s", path, result_file_name)
  print(sprintf("Done! Saving sliding tuning properties to %s", fname))
  save(file=fname, final_result)
}


main_calc_diluted_place_cell_fraction <- function(spike_train, stim_trace, path, verbose=F, result_file_name="diluted_pcs.Rda") {
  
  
  n_iter=10
  processed_st <- preprocess_spike_train(spike_train, stim_trace)
  working_st <- processed_st$working_cells_spike_train
  
  events_per_cell <- apply(working_st, 1, function(r) {sum(r>0)})
  n_20 <- floor(len(events_per_cell) * .2)
  
  
  
  bottom_20_cells <- working_st[order(events_per_cell)[1:n_20],]
  top_20_cells <- working_st[order(events_per_cell, decreasing=T)[1:n_20],]
  
  lower_20_events <- apply(bottom_20_cells, 1, function(r) {sum(r>0)})
  
  pc_frac_bottom_20 <- get_place_cells(bottom_20_cells, stim_trace, cyclic_only = T)
  pc_frac_top_20 <- get_place_cells(top_20_cells, stim_trace, cyclic_only = T)
  
  
  diluted_pval <- c()
  diluted_frac <- c()
  for(itr in 1:n_iter) {
    
    sampled_events <- sample(lower_20_events, n_20)
    
    diluted_spike_train <- 
      lapply(1:nrow(top_20_cells),
             function (cell_idx) {
                  firing_ev <- which(top_20_cells[cell_idx,] > 0)
                  new_diluted_st <- rep(0, times=ncol(top_20_cells))
                  sampled_ev <- sample(firing_ev, sampled_events[cell_idx])
                  new_diluted_st[sampled_ev] <- top_20_cells[cell_idx,sampled_ev]
                  return(new_diluted_st)
                })
    
    diluted_spike_train <- do.call(rbind, diluted_spike_train)
    diluted_pc_frac <- get_place_cells(diluted_spike_train, stim_trace, cyclic_only=T)
    diluted_frac <- c(diluted_frac, diluted_pc_frac$cyclic_fraction)
    diluted_pval <- rbind(diluted_pval, diluted_pc_frac$cyclic_pval)
    
    if (verbose) {
      print(diluted_pc_frac$cyclic_fraction)
    }
    
  }
  
  
  final_result <- list(diluted_frac=diluted_frac,
                       diluted_pval=diluted_pval,
                       lower_20_frac=pc_frac_bottom_20$cyclic_fraction,
                       lower_20_pval=pc_frac_bottom_20$cyclic_pval,
                       top_20_frac=pc_frac_top_20$cyclic_fraction,
                       top_20_pval=pc_frac_top_20$cyclic_pval,
                       events=lower_20_events,
                       top_st=top_20_cells,
                       bottom_st=bottom_20_cells,
                       stim_trace=stim_trace)
  
  print(sprintf("Top: %.3f, Bottom: %.3f, Diluted: %.3f",
                final_result$top_20_frac,
                final_result$lower_20_frac,
                mean(final_result$diluted_frac)))
  
  fname <- sprintf("%s\\%s", path, result_file_name)
  print(sprintf("Done! Saving %s (diluted pcs) to %s",result_file_name,fname))
  save(file=fname, final_result)
}



main_calc_peak_movement_runs <- function(spike_train, 
                                         stim_trace, 
                                         path, 
                                         ext="", 
                                         save_res=T, 
                                         plot_res=F, 
                                         result_file_name="peak_movement_by_runs.Rda",
                                         verbose=F,
                                         use_fractions=F) {
  
  
  # Calculate whose a place cell and who isn't
  place_cells_metadata <- get_place_cells(spike_train,
                                          stim_trace,
                                          verbose=verbose,
                                          cyclic_only = T)
  
  events <- apply(spike_train, 1, function(st) {sum(st>0)})
  npc_ind <- place_cells_metadata$ind[place_cells_metadata$cyclic_non_pcs]
  
  traversals <- get_traversales(stim_trace, nbins=2)
  
  all_traversals_tuning_matrices <- list()
  
  for (i in 1:len(traversals)) {
    
    index_i <- traversals[i]
    
    if (i == len(traversals)) {
      index_f <- len(stim_trace)
    } else {
      index_f <- traversals[i + 1]
    }
    
    
    tuning_mt <- compute_tuning(spike_train[,index_i:index_f], stim_trace[index_i:index_f])
    all_traversals_tuning_matrices <- append(all_traversals_tuning_matrices, list(tuning_mt[[2]]))
  }
  
  
  all_neurons_tuning_all_traversals <- 
  lapply(1:nrow(spike_train),
         function(neur_idx) {

           neur_traversals_tuning <- 
           do.call(rbind, lapply(all_traversals_tuning_matrices, 
                                 function(tuning_mt) {tuning_mt[neur_idx, ]}))
           
           return(neur_traversals_tuning)
           
         })
  
  
  all_data <- list()
  
  used_ind <- c()
  all_neur_idx <- c()
  for (neur_idx in 1:len(all_neurons_tuning_all_traversals)) {
    #print(neur_idx)
    neur_mt <- all_neurons_tuning_all_traversals[[neur_idx]]
  
    
  fired_traversals <- which(apply(neur_mt, 1, function(r) {!all(r == 0)}))
  emt <- neur_mt[fired_traversals,]
  
  if (len(fired_traversals) <= 1) {
    used_ind <- c(used_ind, FALSE)
    next
  } else {
    used_ind <- c(used_ind, TRUE) 
    all_neur_idx <- c(all_neur_idx, neur_idx)
  }
  
  diff_mt <- do.call(rbind, lapply(fired_traversals, function(d) {fired_traversals - d}))
  colnames(diff_mt) <- fired_traversals
  rownames(diff_mt) <- fired_traversals
  
  peak_diff <- function(a,b) {abs(which.max(a) - which.max(b))}
  
  pairwise_cor_mt <- cor(t(emt), t(emt))
  pairwise_peak_mt <- apply(emt, 1, function(a) {apply(emt, 1, function(b) {peak_diff(a,b)})})
  pairwise_cosine_mt <- cosine(t(emt))
  pairwise_euc_mt <- apply(emt, 1, function(a) {apply(emt, 1, function(b) {euclidean(a,b)})})
  
  res_mt <- cbind(diff_mt[upper.tri(diff_mt)],
                  pairwise_cor_mt[upper.tri(pairwise_cor_mt)],
                  pairwise_peak_mt[upper.tri(pairwise_peak_mt)],
                  pairwise_cosine_mt[upper.tri(pairwise_cosine_mt)],
                  pairwise_euc_mt[upper.tri(pairwise_euc_mt)])
  
  
  colnames(res_mt) <- c("TravDiff",
                          "Correlation",
                          "PeakDiff",
                          "Cosine",
                          "Euclid") 
  
  all_data <- append(all_data, list(res_mt))
  
  }
  

  
  final_result <- list(spike_train=spike_train,
                       stim_trace=stim_trace,
                       sliding_df=all_data,
                       place_cells_md=place_cells_metadata,
                       binary_events=events,
                       tuning_mt=all_neurons_tuning_all_traversals,
                       traversals=traversals,
                       cells_with_more_than_one_traversal_fired=used_ind,
                       used_neur_indices=all_neur_idx)

  
  fname <- sprintf("%s\\%s", path, result_file_name)
  print(sprintf("--Done! Saving sliding tuning properties to %s", fname))
  save(file=fname, final_result)
}


main_calc_sliding_tuning_similarity_liron <- function(spike_train, 
                                                      stim_trace, 
                                                      path, 
                                                      ext="", 
                                                      save_res=T, 
                                                      plot_res=F, 
                                                      result_file_name="partial_peaks_sliding_tuning_liron.Rda",
                                                      verbose=F,
                                                      use_fractions=F) {
  
  
  # Calculate whose a place cell and who isn't
  place_cells_metadata <- get_place_cells(spike_train,
                                          stim_trace,
                                          verbose=verbose,
                                          cyclic_only=T)
  
  
  
    window_sizes <- c(200,300,400)
    jumps  <- c(100, 200, 215, 225, 235, 250, 275, 300, 350, 400, 450, 
                500, 600, 650, 750,800,900,1000,1100, 1200, 1300)
  
  
  similarity_results <- list()
  
  # <- <- <- <- metrics <- c("cor", "peaks", "eucl", "cosine", "wasserstein")
  metrics <- c("peaks")
  
  for (window_frac in window_sizes) {
    for (jump_frac in jumps) {
      
      
      window = window_frac
      jump = jump_frac
    
      initial_centroid_index = window
      final_centroid_index = (ncol(spike_train) - 1 - window)
      
      if (final_centroid_index < initial_centroid_index) {
        next
      }

        print(sprintf("(%d:%d) %d - %d - %d", window, jump, initial_centroid_index,
                      final_centroid_index,
                      jump))        
      
      
      centroid_indices <- round(seq(initial_centroid_index, final_centroid_index, by=jump))
      similarity_list <- list() 
      
      if (len(centroid_indices) <= 1) {
        next
      }
      
      for (similarity_metric in metrics) {
        similarity_list[[similarity_metric]]  <- c()
      }
      
      for (cent_idx in 1:(len(centroid_indices) - 1)) {
        
        first_half <- (centroid_indices[cent_idx] - window + 1):centroid_indices[cent_idx]
        second_half <- (centroid_indices[cent_idx + 1] - window + 1):centroid_indices[cent_idx + 1]
        
        if (verbose) {
          print(sprintf("Comparing (%d:%d) to (%d:%d)",
                        min(first_half), max(first_half), min(second_half), max(second_half)))
        }
        first_tuning_curve <- compute_tuning(spike_train[,first_half], stim_trace[first_half])[[2]]
        second_tuning_curve <- compute_tuning(spike_train[,second_half], stim_trace[second_half])[[2]]
        
        
        for (similarity_metric in metrics) {
          similarity_list[[similarity_metric]] <- 
            cbind(similarity_list[[similarity_metric]],
                  tuning_similarity(first_tuning_curve,
                                    second_tuning_curve,
                                    similarity_metric))
        }
      }        
      
      sliding_window <- sprintf("W%f_J%f", window_frac, jump_frac)
      
      print(sprintf("Done Calculating for %s, (%d metrics)", sliding_window,
                    len(similarity_list)))
      print(dim(similarity_list[[1]]))
      
      similarity_results[[sliding_window]] <- similarity_list
      
    }
  }
  
  final_result <- list()
  final_result$similarity <- similarity_results
  final_result$place_cells_metadata <- place_cells_metadata
  final_result$frames <- ncol(spike_train)
  
  fname <- sprintf("%s\\%s", path, result_file_name)
  print(sprintf("Done! Saving sliding tuning properties to %s", fname))
  save(file=fname, final_result)
}



main_calc_place_fraction_by_subsamples_tuning_cor <- function(spike_train, 
                                                   stim_trace, 
                                                   path, 
                                                   ext="", 
                                                   dur_bin=25, 
                                                   number_of_duration_samples=4,
                                                   return_df=F, 
                                                   save_res=T,
                                                   plot_result=T,
                                                   result_file_name="place_cells_by_sample_duration_cor.Rda",
                                                   subsamples_folder="subsampled_matrices") {
  
  prcsd <- preprocess_spike_train(spike_train, stim_trace)
  
  pct_by_duration <- c()
  dur <- floor(seq(20,ncol(spike_train) * dt, length.out=15))
  
  n_of_reps <- floor(1 / ((dur) / max(dur)) * 3)
  
  for (i in 1:len(dur)) {
    duration <- dur[i]
    ind_mat <- generate_sampled_ind(duration, 
                                    ncol(spike_train),
                                    num_of_subsample_repetitions = n_of_reps[i])
    pct <- c()
    
    
    for (i in 1:nrow(ind_mat)) {
      
      
      ind = ind_mat[i,]
      ind <- sample(ind, len(ind))
      
      prcsd <- preprocess_spike_train(spike_train[,ind], stim_trace, verbose = F)
      
      half_size = (len(ind) - len(ind) %% 2) / 2
      
      
      sliced_spike_train_half_1 <- prcsd$working_cells_spike_train[,1:half_size]
      #sliced_spike_train_half_1 <- spike_train[,ind[1:half_size]]
      sliced_stim_trace_half_1 <- stim_trace[ind[1:half_size]]
      
      sliced_spike_train_half_2 <- prcsd$working_cells_spike_train[,(half_size + 1):(half_size * 2)]
      #sliced_spike_train_half_2 <- spike_train[,ind[(half_size + 1):(half_size * 2)]]
      sliced_stim_trace_half_2 <- stim_trace[ind[(half_size + 1):(half_size * 2)]]
      
      #processed_half_1 <- preprocess_spike_train(sliced_spike_train_half_1, sliced_stim_trace_half_1)
      #processed_half_2 <- preprocess_spike_train(sliced_spike_train_half_2, sliced_stim_trace_half_2)
      
      #joint_ind_1 = which(processed_half_1$ind %in% processed_half_2$ind)
      #joint_ind_2 = which(processed_half_2$ind %in% processed_half_1$ind)
      
      # just to make sure no mistakes
      
      joint_ind_1 = joint_ind_1[order(processed_half_1$ind[processed_half_1$ind %in% processed_half_2$ind])]
      joint_ind_2 = joint_ind_2[order(processed_half_2$ind[processed_half_2$ind %in% processed_half_1$ind])]
      
      if (len(joint_ind_1) <= 1) {  next }
      
      #working_cells_spike_train_half_1 <- processed_half_1$working_cells_spike_train[joint_ind_1,]
      #working_cells_spike_train_half_2 <- processed_half_2$working_cells_spike_train[joint_ind_2,]
      
      working_cells_spike_train_half_1 <- sliced_spike_train_half_1 
      working_cells_spike_train_half_2 <- sliced_spike_train_half_2
      
      if (len(joint_ind_1) == 1) {
        
        #working_cells_spike_train_half_1 <- t(as.matrix(processed_half_1$working_cells_spike_train[joint_ind_1,]))
        #working_cells_spike_train_half_2 <- t(as.matrix(processed_half_2$working_cells_spike_train[joint_ind_2,]))
        
        working_cells_spike_train_half_1 <- t(as.matrix(sliced_spike_train_half_1))
        working_cells_spike_train_half_2 <- t(as.matrix(sliced_spike_train_half_2))
      }
      
      tmp <- compute_tuning(working_cells_spike_train_half_1,  sliced_stim_trace_half_1)
      working_tuning_curve_half_1 <- tmp[[2]]
      rm(tmp)
      
      tmp <- compute_tuning(working_cells_spike_train_half_2,  sliced_stim_trace_half_2)
      working_tuning_curve_half_2 <- tmp[[2]]
      rm(tmp)
      
      
     cor_mt <- 
     lapply(1:nrow(working_tuning_curve_half_1),
            function(idx)
            {
              c(cor.test(working_tuning_curve_half_1[idx,],
                         working_tuning_curve_half_2[idx,])$p.value,
                cor(working_tuning_curve_half_1[idx,],
                    working_tuning_curve_half_2[idx,]))
            })
     
     cor_mt <- do.call(rbind, cor_mt)
     cor_mt <- cor_mt[!is.na(cor_mt)[,1] & !is.na(cor_mt)[,2],]
     
     ncor_mt <- cor_mt[cor_mt[,2] > 0,]
     
     if (nrow(ncor_mt) == 1) {
       next
     }
     
     
     pct <- c(pct, sum(ncor_mt[,1] < .05) / nrow(cor_mt))
  
    }
    
    pct_by_duration <- rbind(pct_by_duration, median(pct))
  }
  
  tmp <- as.data.frame(pct_by_duration)
  rownames(tmp) <- dur
  tmp$total_cells <- rep(len(prcsd$ind), times=nrow(tmp))
  tmp$total_frames <- rep(ncol(spike_train), times=nrow(tmp))
  
  if (save_res) {
    
    res_file=sprintf("%s/%s%s",
                       path,                    
                       ext,
                       result_file_name)
    
    print(res_file)
    save(tmp,
         file=res_file)  
  }
  
  
  if (return_df) {
    return(pct_by_duration)
  }
}




main_calc_place_fraction_by_activity_cor <- function(spike_train, stim_trace, path, result_file_name="activity_tuning_df_cor.Rda") {
  
  processed <- preprocess_spike_train(spike_train, stim_trace)
  
  working_cells_spike_train <- processed$working_cells_spike_train
  working_firing_rate <- processed$working_firing_rate
  working_time_bins_per_cells <- processed$working_time_bins_per_cells
  
  
  half_size = (ncol(working_cells_spike_train) - ncol(working_cells_spike_train) %% 2) / 2
  
  
  sliced_spike_train_half_1 <- working_cells_spike_train[,1:half_size]
  #sliced_spike_train_half_1 <- spike_train[,ind[1:half_size]]
  sliced_stim_trace_half_1 <- stim_trace[1:half_size]
  
  sliced_spike_train_half_2 <- working_cells_spike_train[,(half_size + 1):(half_size * 2)]
  #sliced_spike_train_half_2 <- spike_train[,ind[(half_size + 1):(half_size * 2)]]
  sliced_stim_trace_half_2 <- stim_trace[(half_size + 1):(half_size * 2)]
  
  
  
  working_cells_spike_train_half_1 <- sliced_spike_train_half_1 
  working_cells_spike_train_half_2 <- sliced_spike_train_half_2
  
  if (len(joint_ind_1) == 1) {
    
    #working_cells_spike_train_half_1 <- t(as.matrix(processed_half_1$working_cells_spike_train[joint_ind_1,]))
    #working_cells_spike_train_half_2 <- t(as.matrix(processed_half_2$working_cells_spike_train[joint_ind_2,]))
    
    working_cells_spike_train_half_1 <- t(as.matrix(sliced_spike_train_half_1))
    working_cells_spike_train_half_2 <- t(as.matrix(sliced_spike_train_half_2))
  }
  
  tmp <- compute_tuning(working_cells_spike_train_half_1,  sliced_stim_trace_half_1)
  working_tuning_curve_half_1 <- tmp[[2]]
  rm(tmp)
  
  tmp <- compute_tuning(working_cells_spike_train_half_2,  sliced_stim_trace_half_2)
  working_tuning_curve_half_2 <- tmp[[2]]
  rm(tmp)
  
  
  cor_mt <- 
    lapply(1:nrow(working_tuning_curve_half_1),
           function(idx)
           {
             c(cor.test(working_tuning_curve_half_1[idx,],
                        working_tuning_curve_half_2[idx,])$p.value,
               cor(working_tuning_curve_half_1[idx,],
                   working_tuning_curve_half_2[idx,]))
           })
  
  cor_mt <- do.call(rbind, cor_mt)
  
  #cor_mt <- cor_mt[!is.na(cor_mt)[,1] & !is.na(cor_mt)[,2],]
  
  

  bins <- bin(working_time_bins_per_cells, bin_size_active_time_bins)
  pct <- c()
  for (b in levels(bins)) {
    work_cor_mt <- cor_mt[which(bins == b),]
    
    print(dim(work_cor_mt))
    
    if (len(which(bins == b)) == 1) {
      work_cor_mt <- t(as.matrix(work_cor_mt))
    }
    
    pct <- c(pct,
              nrow(work_cor_mt[!is.na(work_cor_mt[,1]) & 
                               !is.na(work_cor_mt[,2]) & 
                               work_cor_mt[,2] > 0 & 
                               work_cor_mt[,1] < .05,]) / nrow(work_cor_mt))
    
  }
  
  
  
  active <- sapply(levels(bins), function(b) 
  {as.numeric(str_split(str_split(b, "\\(")[[1]][2], ",")[[1]][1])})
  
  active <- active[1:len(pct)]
  
  
  if (plot_res) { 
    
    if (save_res) {
      png(sprintf("%s/%s", path, "pct_by_time_bins_single_session_cor.png"))
    }
    
    plot(active, pct, type = "o", frame = FALSE, pch = 19, cex=1, 
         col = "red", 
         xlab = "Active time bins", 
         ylab = "Percentage of place cells", 
         main=str_split(path, "set")[[1]][2])

    
    if (save_res) {
      dev.off()
    }
  }
  
  
  tuning_df <- data.frame(bins=active,
                          percentage=pct)
 
  
  save(tuning_df, file=sprintf("%s/%s",
                               path,
                               result_file_name))
}
