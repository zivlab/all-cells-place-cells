supp_figure_non_place_cell_identification <- function() 
{
  write_path <- sprintf("%s\\supp_figure_classification_sensitivity\\",figures_path)
  dir.create(write_path)

  write_path <- sprintf("%s\\supp_figure_classification_sensitivity\\sensitivity_plots\\",figures_path)
  dir.create(write_path)
  
  directions_to_use = c("Both", "Left", "Right")
  sizes=c(big=3,
          big_1=2.5,
          medium=2,
          medium_1=1.75,
          medium_2=1.5)
  
  
  classification_df_path <- sprintf("%s\\samples\\classification_df\\classification_df.Rda", base_path)
  create_df = F
  
  if (!create_df) {
    load(classification_df_path, verbose=T)
  } else {
    tmp <- get_spike_train_and_stim_trace_from_path(all_data_paths[[2]], 1)
    
    spike_train <- tmp[[1]]
    stim_trace <- tmp[[2]]
    
    fr <- rowMeans(spike_train) / dt
    processed_real <- preprocess_spike_train(spike_train, stim_trace)
    
    true_cells_spike_train <- processed_real$working_cells_spike_train
    true_firing_rate <- processed_real$working_firing_rate
    true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
    classification_df <- c()
    
    stim_trace_orig <- stim_trace
    
    for (p in seq(0.2,0.8,by=.2)) { 
      
      tmp <-
        generate_tuning_curves_cost(n = nrow(spike_train),
                                    percentage = p,
                                    average_width = 0.1, 
                                    sd_width = 0.001,
                                    fixed_fr=fr,
                                    noise=0.02,
                                    double_peak_pct = 0.4,
                                    plot=F,
                                    return_non_place = T)
      
      simulated_tuning_curve <- tmp[[1]]
      non_place_ind <- tmp[[2]]
      
      pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                              true_time_bins_per_cells,
                                              stim_trace_orig,
                                              verbose = T,
                                              jump=.05)    
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace_orig,
                                   factor=pois_factor,
                                   fs=1)
      
      processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace_orig)
      ind_to_use <- processed_generated$ind
      non_place_ind_f <- which(ind_to_use %in% non_place_ind)
      
      
      for (nrep in 1:4) {
        
        stim_trace <- rep(stim_trace_orig, nrep)
        
        for (sim_rep in 1:8) {
          
          print(sprintf("Running simulation rep %d, on length %d, on percent %.2f", sim_rep, nrep, p))
          
          generated_spike_train <- 
            generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                       stim_trace = stim_trace,
                                       factor=pois_factor,
                                       fs=1)
          
          
          generated_firing_rate <- rowMeans(generated_spike_train[ind_to_use,]) / dt
          
          
          tmp <- compute_tuning(generated_spike_train[ind_to_use,], stim_trace)
          gen_stim_prob <- tmp[[1]]
          gen_tuning_curve <- tmp[[2]]
          rm(tmp)
          
          generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
          
          simulated_place_cell_pval <- 
            sapply(c(1:nrow(generated_spike_train[ind_to_use,])),
                   function(i, var2)
                   {compute_place_signif(i,
                                         generated_spike_train[ind_to_use,],
                                         stim_trace,
                                         generated_firing_rate,
                                         generated_SI,
                                         shuffle_type=var2,
                                         num_shuffles = 500,
                                         verbose=F)},
                   var2="c")
          
          
          place_cells <- simulated_place_cell_pval < 0.05
          pct <- sum(place_cells / len(simulated_place_cell_pval))
          correctly_classified <- sum(!place_cells[non_place_ind_f]) / len(place_cells[non_place_ind_f])
          
          classification_df <- rbind(classification_df, c(pct, correctly_classified, sim_rep, nrep, p))
          
          print(classification_df)
          
        }
      }
    }
    
    colnames(classification_df) <- c("MP", "C", "Rep", "Time", "SP")
    classification_df[,"Time"] <- classification_df[,"Time"] * 20
    
    save(file=classification_df_path, classification_df)
  }
  
  cdf <- as.data.frame(classification_df)
  classification_df <- cdf

  
  
  sum_cdf <- 
    ddply(cdf, .(SP), function(sub_df) {
      res <- 
        ddply(sub_df, .(Time), function(subdf_i) {
          
          res_i <- c(colMeans(subdf_i[,c("MP", "C")]), apply(subdf_i[,c("MP", "C")], 2, sd))
          names(res_i) <- c("MP", "C", "sd_MP", "sd_C")
          #      print(res_i)
          return(res_i)
        })
      
      return(res)})
  
  sum_cdf$SP <- sprintf("%d%%", rep(seq(20,80,by=20), each=4))
  
  summary_df <- melt(sum_cdf, id.vars = c("MP", "C", "sd_MP", "sd_C", "Time"))
  
  
  measured_plot <- 
    ggplot(summary_df, aes(y=MP, x=Time, group=value)) +
    geom_line(aes(color=value)) + 
    #geom_errorbar(aes(ymin=MP - sd_MP, ymax=MP + sd_MP, color=value),
    #              width=2,alpha=0.8)  + 
    geom_point(data=cdf, aes(x=Time, y=MP, group=SP), color="gray20", alpha=0.55) +
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title=element_text(size=11),
          text=element_text(size=15, color="black"),
          axis.text=element_text(color="black"),
          axis.ticks = element_line(color="black")) +
    ylab("Naive fraction") + 
    xlab ("Simulation length (Minutes)") +
    ggtitle("Place cells (of active)") + 
    scale_color_brewer("Simulation",palette = "RdYlBu") +
    geom_hline(yintercept = c(0.8, 0.6, 0.4, 0.2), linetype="longdash") +
    ylim(c(0.1,0.9))
  
  
  classified_plot <- 
    ggplot(summary_df, aes(y=C, x=Time, group=value)) +
    geom_line(aes(color=value)) + 
    #geom_errorbar(aes(ymin=C - sd_C, ymax=C + sd_C), width=0.5)  + 
    geom_point(data=cdf, aes(x=Time, y=C, group=SP), color="gray20", alpha=0.55, size=0.3) +
    ylim(c(0,1)) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title=element_text(size=11),
          text=element_text(size=15, color="black"),
          axis.text=element_text(color="black"),
          axis.ticks = element_line(color="black")) +
    scale_color_brewer("Simulation", palette = "RdYlBu") +
    ylab("Correctly classified") +
    ggtitle("Non place cells")
    xlab ("Simulation length (Minutes)")
  
  
  for (size_name in names(sizes)) {
    size = sizes[size_name]
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
  
    pdf(sprintf("%s\\%s\\classification_no_legend.pdf", write_path, size_name),
        height = size,
        width = size)
    plot(classified_plot + theme(legend.position="NA"))
    dev.off()
    
    pdf(sprintf("%s\\%s\\classification.pdf", write_path, size_name),
        height = size,
        width = size)
    plot(classified_plot)
    dev.off()    
    
    pdf(sprintf("%s\\%s\\measured_no_legend.pdf", write_path, size_name),
        height = size,
        width = size)
    plot(measured_plot + theme(legend.position="NA"))
    dev.off()
    
    pdf(sprintf("%s\\%s\\measured.pdf", write_path, size_name),
        height = size,
        width = size)
    plot(measured_plot)
    dev.off()    
  
  }
}

supp_figure_num_shuffles <- function()   
{
  write_path <- sprintf("%s\\supp_figure_classification_sensitivity\\",figures_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\supp_figure_classification_sensitivity\\num_shuffles_plots\\",figures_path)
  dir.create(write_path)
  
  datasets <- list(list(path=3,
                        session=5,
                        name="dataset_3_5"),
                   list(path=2,
                        session=1,
                        name="dataset_2_1"),
                   list(path=11,
                        session=6,
                        name="dataset_11_6"),
                   list(path=8,
                        session=3,
                        name="dataset_8_3"))
  
  num_shuffles_df <- c()
  num_shuffles_df_path <- sprintf("%s\\samples\\classification_df\\num_shuffles_df.Rda", base_path)
  create_df = F
  
  sizes=c(big=2,
          medium_1=1.5,
          medium_2=1.25,
          small=1)
  
  if (!create_df) {
    load(num_shuffles_df_path, verbose=T)
  } else {
    for (ds in datasets[3:4]) {
      tmp <- get_spike_train_and_stim_trace_from_path(all_data_paths[ds$path], ds$session)
      
      spike_train <- tmp[[1]]
      stim_trace <- tmp[[2]]
      
      processed_generated <- preprocess_spike_train(spike_train, stim_trace)
      ind_to_use <- processed_generated$ind
      
      
      firing_rate <- rowMeans(spike_train[ind_to_use,]) / dt
      
      
      tmp <- compute_tuning(spike_train[ind_to_use,], stim_trace)
      stim_prob <- tmp[[1]]
      tuning_curve <- tmp[[2]]
      rm(tmp)
      
      SI <- compute_SI(stim_prob, tuning_curve, firing_rate)
      
      for (nshuffles in c(500,1000,2000,5000)) {
        for (rp in 1:5) {
          
          print(sprintf("%d - %d", rp, nshuffles))
          simulated_place_cell_pval <- 
            sapply(c(1:nrow(spike_train[ind_to_use,])),
                   function(i, var2)
                   {compute_place_signif(i,
                                         spike_train[ind_to_use,],
                                         stim_trace,
                                         firing_rate,
                                         SI,
                                         shuffle_type=var2,
                                         num_shuffles = nshuffles,
                                         verbose=F)},
                   var2="c")
          
          place_cells <- simulated_place_cell_pval < 0.05
          pct <- sum(place_cells / len(simulated_place_cell_pval))
          
          results_df <- rbind(results_df,
                              c(pct, nshuffles, rep, ds$name))
        }
      }
    }
    
    num_shuffles_df <- as.data.frame(num_shuffles_df)
    num_shuffles_df <- num_shuffles_df[,c(1:2,4)]

    save(file=num_shuffles_df_path, num_shuffles_df)
  }
  
  
  colnames(num_shuffles_df) <- c("Percent",
                                 "Number",
                                 "Dataset")
  num_shuffles_df$Percent <- unlist(num_shuffles_df$Percent)
  num_shuffles_df$Number <- unlist(num_shuffles_df$Number)
  num_shuffles_df$Dataset <- unlist(num_shuffles_df$Dataset)
  # nsplot <- 
  
  gall <- list()
  for (dataset_idx in 1:len(unique(num_shuffles_df$Dataset))) {
     dataset_name <- unique(num_shuffles_df$Dataset)[dataset_idx]
     
     df_to_use <- num_shuffles_df[num_shuffles_df$Dataset == dataset_name,]
     
     g <- 
     ggplot(df_to_use) + 
     geom_bar(aes(x=factor(Number), y=Percent, group=Number), width=0.25, position=position_dodge(0.75), stat="summary") +
     geom_point(aes(x=factor(Number), y=Percent, group=Number, color=factor(Number)), position=position_dodge(0.75),
                   size=.5) +
     ylim(c(0,1)) + xlab("") + scale_color_brewer(palette = "Spectral") + 
    theme_light() +     
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text=element_text(size=15),
          axis.text.x = element_text(angle=45, vjust=.4),
          legend.position="NA",
          plot.title=element_text(size=11)) +
    ylab("Naive fraction") +
    xlab("Number of shuffles (#)") +
    ggtitle(sprintf("Dataset %d", dataset_idx))
    #scale_y_continuous(expand=c(0,0))
    
    gall <- append(gall,
                   list(g))
  }

  gf <- do.call(arrangeGrob, gall)
  for (size_name in names(sizes)) {
    size = sizes[size_name]
    
    dir.create(sprintf("%s\\%s",
                       write_path,
                       size_name))
    
    pdf(sprintf("%s\\%s\\num_shuffles_sensitivty.pdf", write_path, size_name),
        height = size * 2,
        width = size * 2)
    plot(gf)
    dev.off()    
    
  } 
}
