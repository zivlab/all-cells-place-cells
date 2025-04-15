figure_2_raster_plots <- function() {
  dir.create(sprintf("%s\\figure_2\\",figures_path))
  write_path <- sprintf("%s\\figure_2\\rasters",figures_path)
  dir.create(write_path)
  
  for (ext in c("", "\\fit_simulation")) {
    for (p in all_data_paths[c(6,1,8,3,4,5)]) {
      for (session in c(sample(c(1:2,9:10), 2),
                        sample(c(6:8,14:16), 2),
                        sample(c(3:5,11:13), 1))) {
        
        a <- get_spike_train_and_stim_trace_from_path(p,session)
        spike_train <- a[[1]]
        stim_trace <- a[[2]]
        
        
        if (ext != "") {
          fr <- rowMeans(spike_train) / dt
          processed_real <- preprocess_spike_train(spike_train, stim_trace)
          
          true_cells_spike_train <- processed_real$working_cells_spike_train
          true_firing_rate <- processed_real$working_firing_rate
          true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
          
          #params <- get_fit_params(path, estimation_path="simulations_2_new_new", likelihood_path = "simulations_3_new", pct_range=seq(0.5,1,by=0.1))
          params <- get_fit_params(sprintf("%s\\equalized\\session_%d", p, session))
          
          print(params)
          simulated_tuning_curve <-
            generate_tuning_curves_cost(n = nrow(spike_train),
                                        percentage = params$pct,
                                        average_width = params$params["average_width"], 
                                        sd_width = params$params["sd_width"],
                                        fixed_fr=fr,
                                        noise=params$params["noise"],
                                        double_peak_pct = params$params["double_peak_pct"],
                                        n_bins=24,
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
          
          pval_by_time = pct_by_sample_size_figures(stim_trace = stim_trace,
                                                    spike_train=generated_spike_train)
        } else {
          simulated_tuning_curve <- list()
          pval_by_time = pct_by_sample_size_figures(stim_trace = stim_trace,
                                                    spike_train=spike_train)
        }
        
        load(sprintf("%s\\equalized\\session_%d%s\\%s", p,session,ext,  "properties.Rda"), verbose=T)
        
        decreasing_rand <- apply(t(apply(pval_by_time$rand, 1, diff)), 1, function(r) {sum(r < 0)})
        
        session_length = 20
        
        
        run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * session_length), 
                             Position=stim_trace * 4) 
        run_df$Position[which(run_df$Position == 4)] <- 0
        actual_times <- (pval_by_time$duration / (ncol(spike_train) * dt)) * session_length
        
        
        if (grepl("CA3", p)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', p))
          mice_str <- substr(p, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', p))
          mice_str <- substr(p, mice_str_index, mice_str_index+3) 
        }
        
        
        dir <- ifelse(grepl("Left", p), "L", "R")
        issim <- ifelse(ext == "", "", "_simulation")
        mice_path <- sprintf("%s\\mice_%s_%s_%d%s",
                             write_path,
                             mice_str,
                             dir,
                             session,
                             issim)
        dir.create(mice_path)
        
        save(file=sprintf("%s\\pval_by_time_df.Rda", mice_path),
             pval_by_time)
        
        pcs <- get_place_cells(spike_train, stim_trace)
        group_indices <- list()
        group_indices$random_non_pcs <- pcs$ind[pcs$random_non_pcs]
        group_indices$random_pcs <- pcs$ind[pcs$random_pcs]
        group_indices$cyclic_pcs <- pcs$ind[pcs$cyclic_pcs]
        group_indices$cyclic_non_pcs <- pcs$ind[pcs$cyclic_non_pcs]
        group_indices$switchers <- pcs$ind[pcs$switchers]
        groups <- names(group_indices)
        
        for (group in groups) {
          mice_group_path <- sprintf("%s\\%s", mice_path, group)
          dir.create(mice_group_path)
          dir.create(sprintf("%s\\regular_size", mice_group_path))
          dir.create(sprintf("%s\\big_size", mice_group_path))
          dir.create(sprintf("%s\\med_size", mice_group_path))
          
          for (idx in group_indices[[group]]) {
            
            #idx <- round(runif(1,1,nrow(spike_train)))
            firing_ind <- which(spike_train[idx,] > 0)
            
            cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
                                  Times=run_df$Time[firing_ind])
            
            pval_df_all <- data.frame(Pvalue=c(pval_by_time$rand[idx,], pval_by_time$cyclic[idx,]),
                                      Time=c(actual_times,actual_times),
                                      `Shuffle type`=c(rep("Random", times=length(actual_times)),
                                                       rep("Cyclic", times=length(actual_times))))
            
            pval_df <- data.frame(Pvalue=c(pval_by_time$cyclic[idx,]),
                                  Time=c(actual_times),
                                  `Shuffle type`=c(rep("Cyclic", times=length(actual_times))))
            
            tuning_curve <- (unlist(lapply(sort(unique(stim_trace)), function (sb) {mean(spike_train[idx,which(stim_trace == sb)]) / dt}))) 
            barp <- barplot(tuning_curve)
            smoothed_tuning <- smooth.spline(barp, tuning_curve, all.knots = T, lambda=1e-4)
            smoothed_tuning$y[smoothed_tuning$y < 0] <- 0
            pos <- sort(unique(stim_trace)) * 4; pos[pos == 4] <- 0
            smoothed_df <- data.frame(FR=smoothed_tuning$y,
                                      Positions=pos)
            
            if (len(simulated_tuning_curve) != 0) {
              sim_barp <- barplot(simulated_tuning_curve[idx,])
              sim_smoothed_tuning <- smooth.spline(sim_barp, simulated_tuning_curve[idx,], all.knots = T, lambda=1e-4)
              sim_smoothed_tuning$y[sim_smoothed_tuning$y < 0] <- 0
              sim_smoothed_df <- data.frame(FR=sim_smoothed_tuning$y / (dt * pois_factor),
                                            Positions=pos)
            }
            
            polygon_bin <- which(pval_df$Pvalue[2:30] < 0.05)[1]
            if(len(polygon_bin) > 0) {
              polygon_df <- polygon_df <- data.frame(x=c(pval_df$Time[polygon_bin],# - (pval_df$Time[polygon_bin] - pval_df$Time[polygon_bin-1]) * 0.5,
                                                         20, 20,
                                                         pval_df$Time[polygon_bin]),# -(pval_df$Time[polygon_bin] - pval_df$Time[polygon_bin-1]) * 0.5),
                                                     y=c(0,0,96,96))
              signif_text_pos = ifelse(polygon_bin - 5 > 0, polygon_bin -5, polygon_bin +5)
            } else {
              polygon_df(x=c(20,20,20,20), y=c(0,0,96,96))
              signif_text_pos = 4
            }
            
            grundf <-
              ggplot(polygon_df, aes(x=x,y=y)) + 
              geom_polygon(fill=adjustcolor("deepskyblue1", alpha=0.2)) +
              #ggplot(run_df) +
              geom_line(data=run_df, aes(x=Time, y=Position),
                        col="gray30") + 
              theme_light() +  
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    legend.position="NA",
                    #panel.border = element_blank(),
                    panel.background = element_blank()) +
              xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
              ylim(0,96) +
              ylab("Position (cm)") +
              xlab("Time (minutes)") +
              geom_vline(xintercept = polygon_df$x[1],
                         linetype="dashed",
                         color="black") +
              geom_point(data=cell_df, aes(x=Times,y=Positions),
                         fill="red",
                         color="red",
                         size=2) + 
              coord_flip()
            
            
            gtuning_df <- ggplot(smoothed_df) + 
              geom_line(aes(x=Positions, y=FR), col="red", size=2) +
              theme_light() +
              ylab("Firing rate (spikes/sec)") + 
              xlab("Position (cm)") +
              xlim(0,96) +
              ylim(0,ifelse(max(smoothed_df$FR) > 1.5, max(smoothed_df$FR), 1.5)) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position="NA",
                    panel.border = element_blank(),
                    panel.background = element_blank()) 
            
            if (len(simulated_tuning_curve) != 0) { 
              
              scale_factor = 0.5
              combined_tuning_df <- rbind(sim_smoothed_df, smoothed_df)
              combined_tuning_df$Type <- rep(c("Defined tuning curve", "Measured"), each=nrow(smoothed_df))
              
              gtuning_df <- ggplot(combined_tuning_df) + 
                geom_line(aes(x=Positions, y=FR, group=Type, color=Type), size=2, alpha=0.8) +
                theme_light() +
                ylab("Firing rate (spikes/sec)") + 
                xlab("Position (cm)") +
                scale_color_manual(breaks=c("Defined tuning curve", "Measured"), values=c("royalblue4", "red")) +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="NA",
                      panel.border = element_blank(),
                      panel.background = element_blank()) +
                ylim(c(0,ifelse(max(smoothed_df$FR) > 1.5, max(smoothed_df$FR), 1.5))) +
                xlim(0,96)
              
            }
            
            gpvaldf <- ggplot(pval_df) + 
              geom_line(aes(x=Time, y=Pvalue, color=`Shuffle.type`)) + 
              geom_point(aes(x=Time, y=Pvalue, color=`Shuffle.type`)) + 
              theme_light() +  
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position="NA",
                    panel.border = element_blank(),
                    panel.background = element_blank()) + 
              ylab(TeX(r'($P_{value}$)')) +
              ylim(0,1) + 
              xlim(0, 20) +
              geom_hline(yintercept = 0.05, linetype="dashed", col="gray60") +
              scale_color_manual(name = "Shuffle type", 
                                 breaks = c("Cyclic", "Random"),
                                 values = c("mediumorchid2", "lightsalmon2"))
            
            thm_big_text <- theme(text=element_text(size=15))
            thm_med_text <- theme(text=element_text(size=14))
            
            grundf_bigger_text <- grundf + thm_big_text
            gtuning_df_bigger_text <- gtuning_df + thm_big_text
            gpvaldf_bigger_text <- gpvaldf + thm_big_text
            
            grundf_med_text <- grundf + thm_med_text
            gtuning_df_med_text <- gtuning_df + thm_med_text
            gpvaldf_med_text <- gpvaldf + thm_med_text
            
            res <- align_plots(grundf, 
                               gtuning_df, 
                               gpvaldf, align="v")
            
            res_big_text <- align_plots(grundf_bigger_text, 
                                        gtuning_df_bigger_text, 
                                        gpvaldf_bigger_text, 
                                        align="v")
            
            res_med_text <- align_plots(grundf_med_text, 
                                        gtuning_df_med_text, 
                                        gpvaldf_med_text, 
                                        align="v")
            
            g <- grid.arrange(res[[1]], res[[2]], res[[3]], heights=c(1.8,1,1))
            g_big <- grid.arrange(res_big_text[[1]], 
                                  res_big_text[[2]], 
                                  res_big_text[[3]], 
                                  heights=c(1.8,1,1))
            
            g_med <- grid.arrange(res_med_text[[1]], 
                                  res_med_text[[2]], 
                                  res_med_text[[3]], 
                                  heights=c(1.8,1,1))
            

            png(sprintf("%s\\regular_size\\neur_%d.png", mice_group_path, idx), unit="in", height=5, width=1.5, res=300); 
            plot(g); dev.off()   
            
            pdf(sprintf("%s\\regular_size\\neur_%d.pdf", mice_group_path, idx), height=5, width=1.5); 
            plot(g); dev.off()   
            
            pdf(sprintf("%s\\big_size\\neur_%d.pdf", mice_group_path, idx), height=5, width=1.5); 
            plot(g_big); dev.off()   
            
            pdf(sprintf("%s\\med_size\\neur_%d.pdf", mice_group_path, idx), height=5, width=1.5); 
            plot(g_med); dev.off()

          }
        }
      }
    }
  }
}

figure_2_bias_plots <- function() {
  write_path <- sprintf("%s/figure_2/",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s/figure_2/bias_plots/",figures_path)
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s/%s", all_data_paths[c(1:8)], "equalized")
  ses_ind <- c(1:16)
  #by_activity <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T, fit_sim = "fit_simulation_new_v_fit\\")
  
  sizes =c(big=3.5,
           medium=2.5,
           medium_1=2.2,
           small=2)
  
  
    by_pval_real <- plot_pval_by_sample_duration_figure(true_data_paths,
                                                        ext = "", 
                                                        fit_sim = "", 
                                                        verbose=F, type1 = F, 
                                                        ses_ind = ses_ind)
      
    by_duration_real <- plot_place_cell_pct_by_sample_duration_figure(true_data_paths, bin_size = 2.5, 
                                                                      ext = "", 
                                                                      cyclic_only=T, 
                                                                      fit_sim = "", 
                                                                      verbose=F, 
                                                                      ses_ind = ses_ind)
    
    big_theme_text = theme(text=element_text(size=15))
    
    for (size_name in names(sizes)) {
      # dir.create(sprintf("%s\\%s", write_path, size_name))
      # pdf(file=sprintf("%s\\%s\\duration.pdf",
      #                  write_path,
      #                  size_name),
      #     height=sizes[[size_name]],
      #     width=sizes[[size_name]])
      # 
      # plot(by_duration_real[[1]] + big_theme_text)
      # dev.off()
      # 
      # pdf(file=sprintf("%s\\%s\\pval.pdf",
      #                  write_path,
      #                  size_name),
      #     height=sizes[[size_name]],
      #     width=sizes[[size_name]])
      # 
      # plot(by_pval_real[[1]] + big_theme_text)
      # dev.off()
      
      pdf(file=sprintf("%s/%s/activity_by_cor_lines.pdf",
                       write_path,
                       size_name),
          height=sizes[[size_name]],
          width=sizes[[size_name]])

      plot(gcor3 + big_theme_text)
      dev.off()
  }
}
