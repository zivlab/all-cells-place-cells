main_calculate_place_cell_fraction_likelihood <- function(spike_train, 
                                                          stim_trace, 
                                                          path, 
                                                          pval_threshold=0.05, 
                                                          sim_reps=100,
                                                          n_shuffles=500,
                                                          plot_example_dists=F,
                                                          binarize=T,
                                                          save_folder="KS_simulations_likelihood_c",
                                                          estimation_folder="JSD_simulations_estimate",
                                                          result_file_name="pct_estimate_df.R",
                                                          shuffle_type="c",
                                                          remove_edges=F,
                                                          save_PNG=F,
                                                          save_PDF=F,
							  override=F) {
  
  print(sprintf("Likelihood analysis on spike_train(%d X %d) stim_trace(%d)", 
                nrow(spike_train),
                ncol(spike_train),
                len(stim_trace)))
  
  print(max(stim_trace))
  
  print(sprintf("path: %s", path))
  
  path_to_save <- sprintf("%s\\%s", path, save_folder) 
  dir.create(path_to_save)
  print(path_to_save)

  if (!override) {
    if ("pct_estimate_df.R" %in% list.files(path_to_save, full.names=F)) {
	print("Continuing")
	return()
    }
  }
  
  
  fr <- rowMeans(spike_train) / dt
  processed_real <- preprocess_spike_train(spike_train, stim_trace)
  
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(true_cells_spike_train, stim_trace, remove_edges = remove_edges)
  stim_prob <- tmp[[1]]
  true_tuning_curve <- tmp[[2]]
  print(dim(true_tuning_curve));
  
  rm(tmp)
  
  true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
  #Compute place tuning cells for cyclic shuffles
  
  if (binarize) {
    binnarized_st <- t(apply(spike_train, 1, function(N) {as.numeric(N > 0)}))
    
    bin_processed_real <- preprocess_spike_train(binnarized_st, stim_trace)
    
    bin_true_cells_spike_train <- bin_processed_real$working_cells_spike_train
    bin_true_firing_rate <- bin_processed_real$working_firing_rate
    bin_true_time_bins_per_cells <- bin_processed_real$working_time_bins_per_cells
    
    tmp <- compute_tuning(bin_true_cells_spike_train, stim_trace)
    bin_stim_prob <- tmp[[1]]
    bin_true_tuning_curve <- tmp[[2]]
    rm(tmp)
    
    bin_true_SI <- compute_SI(bin_stim_prob, bin_true_tuning_curve, bin_true_firing_rate)
    place_cell_pval <- 
      sapply(c(1:nrow(bin_true_cells_spike_train)),
             function(i, var2)
             {compute_place_signif(i,
                                   bin_true_cells_spike_train,
                                   stim_trace,
                                   bin_true_firing_rate,
                                   bin_true_SI,
                                   shuffle_type=var2,
                                   num_shuffles = n_shuffles,
				   verbose=T)},
             var2=shuffle_type)
    
    bin_percent <- sum(place_cell_pval < pval_threshold) / len(place_cell_pval)
  }
  
  place_cell_pval <- 
    sapply(c(1:nrow(true_cells_spike_train)),
           function(i, var2)
           {compute_place_signif(i,
                                 true_cells_spike_train,
                                 stim_trace,
                                 true_firing_rate,
                                 true_SI,
                                 shuffle_type="c",
                                 num_shuffles = n_shuffles,
                                 verbose=T)},
           var2=shuffle_type)
  # 
  
  place_cell_percentage <- sum(place_cell_pval < pval_threshold) / len(place_cell_pval)
  print(place_cell_percentage)
  
  fit_params_df <- extract_simulation_params(path, 
                                             plot=F,  
                                             absolute_min = T, 
                                             estimation_path = estimation_folder, 
                                             pct_range =  seq(0.5, 1, by=0.1),
                                             file_name="pct_estimate_df.R") 
  

  result_df <- c()
 
  # Iterate over all fit parameters and simulate Nreps 
  for (simulated_p in as.character((c(.9, 1.0, .8, .7, .6, .5)))) {
    
    params <- fit_params_df[simulated_p,]
    simulated_percents <- c()
    
    if (plot_example_dists) { 
      to_plot_example = T  
    }
    
    for (i in 1:sim_reps) {
      
      tuning_curve <-
        generate_tuning_curves_cost(n = nrow(spike_train),
                                    percentage = as.numeric(simulated_p),
                                    average_width = params["average_width"], 
                                    sd_width = params["sd_width"],
                                    fixed_fr=fr,
                                    noise=params["noise"],
                                    n_bins=ifelse(remove_edges, 20, 24),
                                    double_peak_pct = params["double_peak_pct"],
                                    plot=F)
      
      # pois_factor <- currate_spike_train_cost(tuning_curve, 
      #                                         true_time_bins_per_cells,
      #                                         stim_trace,
      #                                         verbose = F)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=1.0,
                                   fs=1)
      
      
      
      
      ### Generate a spike train
      processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=T)
      generated_cells_active_st <- processed_generated$working_cells_spike_train
      generated_firing_rate <- processed_generated$working_firing_rate
      
      
      tmp <- compute_tuning(generated_cells_active_st, stim_trace, remove_edges = remove_edges)
      gen_stim_prob <- tmp[[1]]
      gen_tuning_curve <- tmp[[2]]
      print(dim(gen_tuning_curve))
      print(remove_edges)
      rm(tmp)
      
      generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
      
      
      if (F && to_plot_example) {
        graphs <- plot_gen_place_cell_histograms(average_width = params["average_width"],
                                                 sd_width = params["sd_width"],
                                                 noise_level = params["noise"],
                                                 double_peak_percent = params["double_peak_pct"],
                                                 true_spike_train = spike_train,
                                                 stim_trace = stim_trace,
                                                 place_cell_percentage = as.numeric(simulated_p),
                                                 pval_func = wilcox.test,
                                                 hist_density = T)
        
        fg <- arrangeGrob(graphs[[1]], graphs[[2]], graphs[[3]], graphs[[4]], nrow=2)
        
        png(sprintf("%s/example_dist_hists_%s.png", path_to_save, simulated_p),
            unit="in",
            res=500,
            width=8,
            height=8)
        plot(fg)
        dev.off()
        
        pdf(file=sprintf("%s/example_dist_hists_%s.pdf", path_to_save, simulated_p),
           width=8,
           height=8)
         plot(fg)
         dev.off()       
        to_plot_example = F
      }
      
      simulated_place_cell_pval <- 
        sapply(c(1:nrow(generated_cells_active_st)),
               function(i, var2)
               {compute_place_signif(i,
                                     generated_cells_active_st,
                                     stim_trace,
                                     generated_firing_rate,
                                     generated_SI,
                                     shuffle_type="c",
                                     num_shuffles = n_shuffles)},
               var2=shuffle_type)
      # 
      
      simulated_place_cell_percentage <- sum(simulated_place_cell_pval < pval_threshold) / len(simulated_place_cell_pval)
      simulated_percents <- c(simulated_percents,
                              simulated_place_cell_percentage)
      print(sprintf("Rep %d, got percent %f", i, simulated_percents))
    }
    
    result_df <- rbind(result_df, simulated_percents)
    print(result_df)
  } 
  
  
  rownames(result_df) <- rev(rownames(fit_params_df))
  
  likelihood <- 
    unlist(apply(result_df, 1, 
                 function(sim_vec)
                 {
                   kdf_func <- approxfun(density(sim_vec, width=0.08, from=0, to=1))
                   
                   if (place_cell_percentage > mean(sim_vec)) { 
                     pval <- integrate(kdf_func, place_cell_percentage, 1)$value
                   } else {
                     pval <- integrate(kdf_func, 0, place_cell_percentage)$value
                   }
                   return(list(pval))
                 }))
  
  mdf <- melt(result_df)
  mdf$Var1 <- factor(mdf$Var1)
  
  bar_df <- data.frame(Likelihood=likelihood, 
                       Percentage=factor(sprintf("%d%%", as.numeric(rownames(fit_params_df)) * 100),
                                         levels=sprintf("%d%%", seq(40,100,by=10))))
  
  limits_of_plots <- c(min(mdf$value) - 0.05, max(mdf$value) + 0.1)
  
  gbox <- 
    ggplot(mdf) +
    geom_boxplot(aes(y=value, group=Var1, fill=Var1, x=Var1), size=1, alpha=0.9) +
    
    scale_fill_brewer(palette="Spectral") +
    theme_light() +
    xlab("True simulation percentage (%)") +
    ylab("Measured simulation percentage (%)") +
    base_plot_theme +
    geom_hline(yintercept=place_cell_percentage, linetype="dashed", col="gray60", size=2) +
    geom_text(y=place_cell_percentage + 0.01, x=2.5, label="Measured place cell percentage (actual data)",col="gray60",size=3) 
  
  gdot <- 
    ggplot(mdf) +
    geom_boxplot(aes(y=value, group=Var1, fill=Var1, x=Var1), size=1, alpha=0.4, fill="white") +
    geom_dotplot(aes(y=value, group=Var1, fill=Var1,  x=Var1), binaxis='y', 
                 stackdir='center', dotsize = 0.4) +
    
    scale_fill_brewer(palette="Spectral") +
    theme_light() +
    xlab("True simulation percentage (%)") +
    ylab("Measured simulation percentage (%)") +
    base_plot_theme +
    geom_hline(yintercept=place_cell_percentage, linetype="dashed", col="gray60", size=2) +
    geom_text(y=place_cell_percentage + 0.01, x=2.5, label="Measured place cell percentage (actual data)",col="gray60",size=3)
  
  gbar <- ggplot(bar_df, aes(x=Percentage, y=Likelihood)) +
    geom_bar(stat="identity", fill="gray50", color="black") + 
    theme_light() +
    base_plot_theme
  
  gdensity <- 
    ggplot(mdf,  aes(y = Var1, x = value)) + 
    geom_density_ridges(bandwidth=0.02, aes(fill=Var1), size=1, alpha=0.8, scale=0.9) + 
    coord_flip() + 
    scale_fill_brewer(palette="Spectral") +
    theme_light() +
    ylab("True simulation percentage (%)") +
    xlab("Measured simulation percentage (%)") +
    base_plot_theme +
    geom_vline(xintercept=place_cell_percentage, linetype="dashed", col="gray60", size=2) +
    geom_text(y=place_cell_percentage + 0.01, x=2.5, label="Measure place cell percentage (actual data)",col="gray60") 
  
  
  if (binarize) {
    g_bin_box <- 
      ggplot(mdf) +
      # geom_dotplot(aes(y=value, group=Var1, fill=Var1,  x=Var1), binaxis='y', 
      #              stackdir='center', dotsize = 0.3) +
      geom_boxplot(aes(y=value, group=Var1, fill=Var1, x=Var1), size=1, alpha=0.9) +
      
      scale_fill_brewer(palette="Spectral") +
      theme_light() +
      xlab("True simulation percentage (%)") +
      ylab("Measured simulation percentage (%)") +
      base_plot_theme+
      geom_hline(yintercept=place_cell_percentage, linetype="dashed", col="gray60", size=2) +
      geom_text(y=place_cell_percentage + 0.01, x=2.5, label="Measured place cell percentage (actual data)",col="gray60",size=3) +
      geom_hline(yintercept=bin_percent, linetype="dashed", col="gray60", size=2)
    
    
    g_bin_density <- 
      ggplot(mdf,  aes(y = Var1, x = value)) + 
      geom_density_ridges(bandwidth=0.02, aes(fill=Var1), size=1, alpha=0.8, scale=0.9) + 
      coord_flip() + 
      scale_fill_brewer(palette="Spectral") +
      theme_light() +
      ylab("True simulation percentage (%)") +
      xlab("Measured simulation percentage (%)") +
      base_plot_theme +
      geom_vline(xintercept=place_cell_percentage, linetype="dashed", col="gray60", size=2) +
      geom_text(y=place_cell_percentage + 0.01, x=2.5, label="Measure place cell percentage (actual data)",col="gray60")  +
      geom_vline(xintercept=bin_percent, linetype="dashed", col="gray60", size=2)
    
    gcomb_bin <-  arrangeGrob(g_bin_box + ylim(limits_of_plots[1], limits_of_plots[2]), 
                              g_bin_density + scale_x_continuous(limits=limits_of_plots), 
                              gbar, 
                              nrow=1)
    
  }
  
  final_df <- cbind(result_df, likelihood, rep(place_cell_percentage, times=len(likelihood)))
  print(dim(final_df))
  colnames(final_df) <- c(sprintf("rep_%d", 1:sim_reps), "likelihood", "measured_pct")
  
  gcomb <- arrangeGrob(gbox + ylim(limits_of_plots[1], limits_of_plots[2]), 
                       gdensity + scale_x_continuous(limits=limits_of_plots), 
                       gbar, 
                       nrow=1)
  
  save(final_df, file=sprintf("%s\\%s", path_to_save, result_file_name))
  
  if (save_PDF) {
    if (binarize) {
      pdf(file=sprintf("%s\\bin_combined_plot.pdf", path_to_save),
          width=24,
          height=8)
      plot(gcomb_bin)
      dev.off()
    }
    
    pdf(file=sprintf("%s\\combined_plot.pdf", path_to_save),
        width=24,
        height=8)
    plot(gcomb)
    dev.off()
    
    pdf(file=sprintf("%s\\barplot.pdf", path_to_save),
        width=8,
        height=8)
    plot(gbar)
    dev.off()
    
    pdf(file=sprintf("%s\\boxplot.pdf", path_to_save),
        width=8,
        height=8)
    plot(gbox)
    dev.off()
    
    pdf(file=sprintf("%s\\density.pdf", path_to_save),
        width=8,
        height=8)
    plot(gdensity)
    dev.off()
    
    pdf(file=sprintf("%s\\dotplot.pdf", path_to_save),
        width=8,
        height=8)
    plot(gdot)
    dev.off()
  }
  
  if (save_PNG) {
    if (binarize) {
      png(sprintf("%s\\bin_combined_plot.png", path_to_save),
          unit="in",
          res=500,
          width=24,
          height=8)
      plot(gcomb_bin)
      dev.off()
    }
    
    png(sprintf("%s\\combined_plot.png", path_to_save),
        unit="in",
        res=500,
        width=24,
        height=8)
    plot(gcomb)
    dev.off()
    
    # png(sprintf("%s\\barplot.png", path_to_save),
    #     unit="in",
    #     res=500,
    #     width=8,
    #     height=8)
    # plot(gbar)
    # dev.off()
    # 
    # png(sprintf("%s\\boxplot.png", path_to_save),
    #     unit="in",
    #     res=500,
    #     width=8,
    #     height=8)
    # plot(gbox)
    # dev.off()
    # 
    # png(sprintf("%s\\density.png", path_to_save),
    #     unit="in",
    #     res=500,
    #     width=8,
    #     height=8)
    # plot(gdensity)
    # dev.off()
    # 
    # png(sprintf("%s\\dotplot.png", path_to_save),
    #     unit="in",
    #     res=500,
    #     width=8,
    #     height=8)
    # plot(gdot)
    # dev.off()
  }
  
}

main_estimate_simulation_parameters <- function(spike_train, 
                                                stim_trace, 
                                                path, 
                                                computed_kernel_path="/Users/itayta/Downloads/yaniv_place_cells_project/NewCleanCodeBase/param_space_cov_mat.R", 
                                                out_folder = "JSD_simulations_estimate",
                                                result_file_name="pred_df.R",
                                                tweak=F,
						                                    remove_edges=F,
                                                compute_kernel=F,
                                                plot_contours=F,
                                                use_KS=F) {
  
  
  
  dir.create(sprintf("%s/%s", path, out_folder))
  
  param_space <- list(double_peak_pct=seq(0, 0.5,length.out=7),
                      average_width=seq(0.01, 0.17, by=0.01),
                      noise=seq(0.01, 0.08, by=0.01),
                      sd_width=seq(0.0001, 0.015, length.out=10))
  
  print(param_space)
  
  
  parameter_set <- as.matrix(cross_df(param_space))
  colnames(parameter_set) <- c("double_peak_pct", "average_width", "noise", "sd_width")
  
  if (use_KS) {
    objective_func=place_cell_params_sim_objective
  } else {
    objective_func=place_cell_params_sim_objective_new
  }
  
  if (!compute_kernel) {
    print("Reading computed kernel for params covariance mat")
    load(computed_kernel_path)
  } else {
    print("Computing kernel for params covariance mat")
    cov_mat_pred <- rbf_kernel(parameter_set, parameter_set, sigma=0.35)
    save(cov_mat_pred, file = computed_kernel_path)
  }
  
  true_percentage_range <- rev(seq(0.5, 1, by=0.1))
  # true_percentage_range <- c(0.7,.6,.5,.8,.9,1)
  #true_percentage_range <- c(.9,1)
  
  for (true_p in true_percentage_range) {
    work_folder <- sprintf("%s/%s/param_fit_%.3f\\", path, out_folder, true_p)    
    res_file_name <- sprintf("%s/pred_df.R", work_folder)

    print("pred_df.R" %in% list.files(work_folder, full.names=F))

    # if ("pred_df.R" %in% list.files(work_folder, full.names=F)) {
    #     print("Continuing")
    #     next
    # }
    
    
    dir.create(work_folder)
    mres <- 
      bayesian_optimizer(objective = objective_func,
                         parameter_set = parameter_set,
                         cov_mat_params = cov_mat_pred,
                         n0 = 30,
                         max_iterations = 35,
                         improval_step = 10,
                         sigma_n = 0.35,
                         aquisition = "expected",
                         true_spike_train = spike_train,
                         stim_trace = stim_trace,
                         use_SI=T,
                         remove_edges=remove_edges,
                         percentage = true_p,
                         dt=dt,
                         dist_fun=ifelse(use_KS, "KS", "JSD"))
    
    predicted_df <- cbind(parameter_set, t(mres))
    colnames(predicted_df)[5:24] <- sprintf("predicted_%d", 1:20)
    save(predicted_df, file=sprintf("%s\\%s", work_folder, result_file_name))
    
    if (plot_contours) {
      contour_maps(mres, parameter_set, work_folder)
    }
  }
  

  }




tmp_generation = function() {
  
  fr <- rowMeans(spike_train) / dt
  tuning_curve <-
    generate_tuning_curves_cost(n = nrow(spike_train),
                                percentage = 1,
                                average_width =.12, 
                                sd_width = .02,
                                fixed_fr=unlist(lapply(1:1, function(i) {sample(fr, nrow(spike_train))})),
                                noise=.02,
                                n_bins=ifelse(F, 20, 24),
                                double_peak_pct = .1,
                                plot=F)
  
  
  
  norm_tc <- t(apply(tuning_curve, 1, function(r) {(r - min(r))/(max(r) - min(r))}))
  ph(norm_tc[order(unlist(apply(norm_tc, 1, which.max))),], color=viridisLite::inferno(50))
  
  
  old_st = spike_train
  processed_real <- preprocess_spike_train(old_st, stim_trace)
  
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  pois_factor <- currate_spike_train_cost(tuning_curve, 
                                          true_time_bins_per_cells,
                                          stim_trace,
                                          verbose = T)
  
  generated_spike_train <- 
    generate_spike_trains_cost(tuning_curves = tuning_curve,
                               stim_trace = stim_trace,
                               factor=pois_factor,
                               fs=1)
  
  
  path = "/Users/itayta/Downloads/yaniv_place_cells_project//simulated_dataset09/"
  dir.create(path)
  main_estimate_simulation_parameters()
  
  a= get_spike_train_and_stim_trace_from_path(all_data_paths[5], 7)
}