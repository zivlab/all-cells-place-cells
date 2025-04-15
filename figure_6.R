figure_6_non_place_cells_rasters <- function() {
  dir.create(sprintf("%s\\figure_6\\",figures_path))
  write_path <- sprintf("%s\\figure_6\\rasters",figures_path)
  dir.create(write_path)
  
  
  for (p in all_data_paths[c(2,4,6,8,10,12,14,16)]) {
      for (session in c(1,2,3,4,9,10)) {
        
        session_length <- 20
        a <- get_spike_train_and_stim_trace_from_path(p,session)
        spike_train <- a[[1]]
        stim_trace <- a[[2]]
        
        for (ext in c("")) {
          
          simulated_tuning_curve <- list()
          
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
            
            spike_train <- generated_spike_train
            
          }
          
          
          tmp <- compute_tuning(spike_train, stim_trace)
          SI <- compute_SI(tmp[[1]], tmp[[2]], rowMeans(spike_train) / dt)[[1]]
          
          run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * session_length), 
                               Position=stim_trace * 4) 
          run_df$Position[which(run_df$Position == 4)] <- 0
          
          
          
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
          
          
          pcs <- get_place_cells(spike_train, stim_trace)
          group_indices <- list()
          group_indices$random_non_pcs <- pcs$ind[pcs$random_non_pcs]
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
              
              
              tuning_curve <- (unlist(lapply(sort(unique(stim_trace)), function (sb) {mean(spike_train[idx,which(stim_trace == sb)]) / dt}))) 
              barp <- barplot(tuning_curve)
              smoothed_tuning <- smooth.spline(barp, tuning_curve, all.knots = T, lambda=1e-4)
              smoothed_tuning$y[smoothed_tuning$y < 0] <- 0
              pos <- sort(unique(stim_trace)) * 4; pos[pos == 4] <- 0
              smoothed_df <- data.frame(FR=smoothed_tuning$y,
                                        Positions=pos)
              
              smoothed_df$FR <- smoothed_df$FR / sum(smoothed_df$FR)
              
              if (len(simulated_tuning_curve) != 0) {
                sim_barp <- barplot(simulated_tuning_curve[idx,])
                sim_smoothed_tuning <- smooth.spline(sim_barp, simulated_tuning_curve[idx,], all.knots = T, lambda=1e-4)
                sim_smoothed_tuning$y[sim_smoothed_tuning$y < 0] <- 0
                sim_smoothed_df <- data.frame(FR=sim_smoothed_tuning$y / (dt * pois_factor),
                                              Positions=pos)
                
                sim_smoothed_df$FR <- sim_smoothed_df$FR / sum(sim_smoothed_df$FR)
              }
              
              
              title = sprintf("SI (%f) E(%d)",
                              SI[[idx]],
                              len(firing_ind))
              
              grundf <-
                ggplot(run_df, aes(x=Time, y=Position)) + 
                #geom_line(col="gray30") + 
                theme_light() +  
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.position="NA",
                      axis.ticks=element_line(color="black"),
                      axis.text=element_text(color="black"),
                      panel.border = element_rect(color="black"),
                      panel.background = element_blank(),
                      plot.title=element_text(size=7)) +
                xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
                ylim(0,96) +
                ylab("Position (cm)") +
                xlab("Time (minutes)") +
                # geom_vline(xintercept = polygon_df$x[1],
                #            linetype="dashed",
                #            color="black") +
                geom_point(data=cell_df, aes(x=Times,y=Positions),
                           fill="#DA1C5C",
                           color="#DA1C5C",
                           stroke=0,
                           size=2,
                           alpha=.8) + 
                coord_flip() +
                ggtitle(title)
              
              
              gtuning_df <- ggplot(smoothed_df) + 
                geom_line(aes(x=Positions, y=FR), col="red", size=2) +
                theme_light() +
                ylab("Firing rate (spikes/sec)") + 
                xlab("Position (cm)") +
                xlim(0,96) +
                ylim(0,ifelse(max(smoothed_df$FR) > .18, max(smoothed_df$FR), .18)) +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.ticks=element_line(color="black"),
                      axis.text=element_text(color="black"),
                      legend.position="NA",
                      panel.border = element_blank(),
                      panel.background = element_blank()) 
              
              gpvaldf <- ggplot() + 
                theme_light() +  
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="NA",
                      panel.border = element_blank(),
                      panel.background = element_blank()) + 
                ylab(TeX(r'($P_{value}$)')) +
                ylim(0,1) + 
                xlim(0, 20)
              
              
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
                  ylim(c(0,ifelse(max(smoothed_df$FR) > .18, max(smoothed_df$FR), .18))) +
                  xlim(0,96)
                
                
                
              }
              
              #### Note the pvalue dfs are empty just for alignment of plots
              
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

figure_6_corlineplots <- function() {
  
  load(sprintf("%s\\samples\\tuning_dataframes\\KS_stability_dataframes.R", 
               base_path), 
       verbose=T)
  
  write_path_base <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path_base)
  write_path_base <- sprintf("%s\\figure_6\\07_03_corline_plots\\",figures_path)
  dir.create(write_path_base)
  

  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.3,
          small=2)
  
  
  session_groups <- list(all=1:16)
    # first_2=c(1:2,9:10),
    # first_3=c(1:3,9:11),
    # first_4=c(1:4,9:12))
    # 

  
  
  for (session_group_name in names(session_groups)) {
    
    write_path <- sprintf("%s\\%s", write_path_base, session_group_name)
    dir.create(write_path)
    dir.create(sprintf("%s\\big\\", write_path))
    dir.create(sprintf("%s\\medium_1\\", write_path))
    dir.create(sprintf("%s\\medium\\", write_path))
    dir.create(sprintf("%s\\small\\", write_path))
    sessions_to_use <- session_groups[[session_group_name]]
    
    for (df_name in names(list_of_data_frames)) {
      corr_df <- list_of_data_frames[[df_name]]
      corr_df <- corr_df[corr_df$Session %in% sessions_to_use,]
      
      for (subsample in unique(corr_df$Subsamples)) {
        subsample_cor_df <- corr_df[corr_df$Subsamples == subsample,]
        
        plts <- list()
        plts_small <- list()
        
        properties_to_use <- c("peaks",
                               "time_bins",
                               "SI",
                               "spatial_bins")
        #properties_to_use <- unique(subsample_cor_df$Property)[c(1,2,4,6)]
        
        for (subfield in unique(subsample_cor_df$Subfield)) {
          for (prop in properties_to_use)  { 
            
            fdf <- data.frame()
            fdf2 <- data.frame()
            
            pcs_real <- c()
            non_pcs_real <- c()
            pcs_sim <- c()
            non_pcs_sim <- c()
            
            real <- c()
            sim <- c()
            
            pcs <- c()
            non_pcs <- c()
            
            for (mice in unique(subsample_cor_df$Mice)) {
              for (session in unique(subsample_cor_df$Session)){
                for (dir in unique(subsample_cor_df$Direction)) {
                  for (group in c("cyclic_pcs", "cyclic_non_pcs")) {
                    
                    
                    working_df <- subsample_cor_df[subsample_cor_df$Subfield == subfield &
                                                     subsample_cor_df$Mice == mice &   
                                                     subsample_cor_df$Session == session &
                                                     subsample_cor_df$Direction == dir &
                                                     subsample_cor_df$Property == prop &
                                                     subsample_cor_df$Group %in%  group,]
                    
                    if (nrow(working_df) < 2) {
                      next
                    }
                    
                    df <- 
                      data.frame(Real_corr=as.numeric(working_df$Corr[working_df$isSim=="Real"]),
                                 MOS_corr=as.numeric(working_df$Corr[working_df$isSim=="Simulation"]),
                                 Direction=dir,
                                 Session=session,
                                 Mice=mice,
                                 Propery=prop,
                                 Group=group)
                    
                    df2 <- 
                      data.frame(Corr=c(as.numeric(working_df$Corr[working_df$isSim=="Real"]),
                                        as.numeric(working_df$Corr[working_df$isSim=="Simulation"])),
                                 Direction=c(dir, dir),
                                 Session=c(session, session),
                                 Mice=c(mice, mice),
                                 Propery=c(prop, prop),
                                 Group=ifelse(group == "cyclic_pcs", "PCs", "Non-PCs"),
                                 isSim=c("Real", "Simulation"))
                    
                    
                    fdf <- rbind(fdf,df)
                    fdf2 <- rbind(fdf2, df2)
                    
                    if (group == "cyclic_pcs") {
                      pcs_real <- c(pcs_real, df2[df2$isSim == "Real","Corr"])
                      pcs_sim <- c(pcs_sim, df2[df2$isSim != "Real","Corr"])
                      
                      pcs <- c(pcs, c(df2[,"Corr"]))
                    } else {
                      non_pcs_real <- c(non_pcs_real, df2[df2$isSim == "Real","Corr"])
                      non_pcs_sim <- c(non_pcs_sim, df2[df2$isSim != "Real","Corr"]) 
                      
                      non_pcs <- c(non_pcs, c(df2[,"Corr"]))
                    }
                    
                    real <- c(real, df2[df2$isSim == "Real","Corr"])
                    sim <- c(sim, df2[df2$isSim != "Real","Corr"])
                    
                    
                  }
                }
              }
            }
            t
            all_corr_vec <- subsample_cor_df[subsample_cor_df$Property == prop,"Corr"]
            #all_corr_vec <- corr_df$Corr[corr_df$Corr]
            
            # lims <- c(min(all_corr_vec[!is.na(all_corr_vec)]),
            #           max(all_corr_vec[!is.na(all_corr_vec)]))
            
            lims <- c(min(fdf[,c(1,2)]),
                      1)
            
            ldf <- data.frame(x=c(min(lims),max(lims)),
                              y=c(min(lims),max(lims)))
            
            g <- 
              ggplot(fdf, aes(x=Real_corr, y=MOS_corr)) + 
              geom_line(data=ldf,
                        aes(x=x, y=y),
                        linetype="longdash",
                        col="gray40")+
              geom_point(aes(col=Group)) + 
              ylim(lims) + 
              xlim(lims) + 
              xlab("Stability (data)") + 
              ylab("Stability (simulation)") + 
              ggtitle(sprintf("%s - %s", subfield, prop)) +
              scale_color_manual(values=c(adjustcolor(c("#A6A4A4", "#405EAB", alpha=0.8)))) +
              base_plot_theme
            
            
            
            
            
            plts_small <- append(plts_small,
                                 list(g))
            
            plts <- append(plts,
                           list(g + theme(text=element_text(size=13.5), plot.title=element_text(size=10))))
          }
        }
        
        plts_small$nrow <- 2
        gf_small <- do.call(grid.arrange, plts_small)
        
        plts$nrow <- 2
        gf <- do.call(grid.arrange, plts)
        
        for (size_name in names(sizes)) {
          size = sizes[size_name]
          dir.create(sprintf("%s\\%s\\%s", 
                             write_path,
                             size_name,
                             df_name))
          
          pdf(file=sprintf("%s\\%s\\%s\\corline_subsamples_%s.pdf",
                           write_path,
                           size_name,
                           df_name, 
                           subsample),
              height=size * 2,
              width=size * 4)
          
          
          plot(gf)
          dev.off()
          
          
          pdf(file=sprintf("%s\\%s\\%s\\small_text_corline_subsamples_%s.pdf",
                           write_path,
                           size_name,
                           df_name, 
                           subsample),
              height=size * 2,
              width=size * 4)
          
          
          plot(gf_small)
          dev.off()
        }
      }
    }
  }
}

figure_6_stability_comparisions <- function() {  
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\stability_comparision\\",figures_path)
  dir.create(write_path)
  
  load(sprintf("%s\\samples\\tuning_dataframes\\KS_stability_dataframes.R", 
               base_path), 
       verbose=T)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  sessions_to_use=c(1:16)
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  
  metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
  metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
  colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
  
  dir.create(sprintf("%s\\big\\", write_path))
  dir.create(sprintf("%s\\medium_1\\", write_path))
  dir.create(sprintf("%s\\medium\\", write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes=c(big=3.5,
          medium=3,
          medium_1=2.75,
          small=2.5)
  
  groups_to_use = c(1:2)
  
  comp_groups <-
    c("PC_REAL_PC_SIM",
      "NPC_REAL_NPC_SIM",
      "PC_REAL_NPC_REAL",
      "PC_SIM_NPC_SIM")
  
  statistics_df <- data.frame()
  values_df <- data.frame()
  
  for (df_name in names(list_of_data_frames)) {
    corr_df <- list_of_data_frames[[df_name]]
    
    for (subsample in unique(corr_df$Subsamples)) {
      subsample_cor_df <- corr_df[corr_df$Subsamples == subsample,]
      
      metric_plot_lists <- list()
      
      properties_to_use <- c("peaks",
                             "time_bins",
                             "SI",
                             "spatial_bins")
      metric_functions <- list(KL=KL_dist,
                               JSD=JSD_dist,
                               dprime=dprime_dist,
                               wasserstein=wasserstein_dist,
                               MSE=mse)
      
      for (metric_name in names(metric_functions)) {
        metric_plot_lists[[metric_name]] <- list()
      }
      
      #properties_to_use <- unique(subsample_cor_df$Property)[c(1,2,4,6)]
      subfield_list <- list()
      npc_error_list <- list()
      npc_by_session_error_list <- list()
      
      for (subfield in unique(subsample_cor_df$Subfield)) {
        for (prop in properties_to_use)  {
          
          comparisions_df <- data.frame()
          subsample_cor_df$Session <- as.numeric(subsample_cor_df$Session)
          
          mice_all <- c()
          sessions <- c()
          pcs_real_all <- c()
          non_pcs_real_all <- c()
          pcs_sim_all <- c()
          non_pcs_sim_all <- c()
          est_all <- c()
          
          for (session in 1:16){
            pcs_real <- c()
            non_pcs_real <- c()
            pcs_sim <- c()
            non_pcs_sim <- c()
            
            for (mice in unique(subsample_cor_df$Mice)) {
              for (env in c(0)) {
                for (dir in unique(subsample_cor_df$Direction)) {
                 
                  
                  for (group in c("cyclic_pcs", "cyclic_non_pcs")) {

                    

                                        
                    working_df <- subsample_cor_df[subsample_cor_df$Subfield == subfield &
                                                     subsample_cor_df$Mice == mice &   
                                                     subsample_cor_df$Session == (session + env) &
                                                     subsample_cor_df$Direction == dir &
                                                     subsample_cor_df$Property == prop &
                                                     subsample_cor_df$Group %in%  group,]
                    
                    if (nrow(working_df) < 2) {
                      next
                    }
                    
                    
                    est  <- metadata_df[metadata_df$Session == session & 
                                          metadata_df$Mice == mice & 
                                          metadata_df$Direction == dir,]$Est
                    # Save only once
                    if (group == "cyclic_pcs") {
                      est_all = c(est_all, est)
                      sessions <- c(sessions, session)
                      mice_all <- c(mice_all, mice)
                    }
                    
                    
                    df2 <-
                      data.frame(Corr=c(as.numeric(working_df$Corr[working_df$isSim=="Real"]),
                                        as.numeric(working_df$Corr[working_df$isSim=="Simulation"])),
                                 Direction=c(dir, dir),
                                 Session=c(session, session),
                                 Mice=c(mice, mice),
                                 Propery=c(prop, prop),
                                 Group=ifelse(group == "cyclic_pcs", "PCs", "Non-PCs"),
                                 isSim=c("Real", "Simulation"),
                                 isAll=ifelse(est == 1, "100", "Not100"))
                    
                    
                    if (group == "cyclic_pcs") {
                      pcs_real <- c(pcs_real, df2[df2$isSim == "Real","Corr"])
                      pcs_real_all <- c(pcs_real_all, df2[df2$isSim == "Real","Corr"])
                      pcs_sim <- c(pcs_sim, df2[df2$isSim != "Real","Corr"])
                      pcs_sim_all <- c(pcs_sim_all, df2[df2$isSim != "Real","Corr"])
                    } else {
                      non_pcs_real <- c(non_pcs_real, df2[df2$isSim == "Real","Corr"])
                      non_pcs_real_all <- c(non_pcs_real_all, df2[df2$isSim == "Real","Corr"])
                      non_pcs_sim <- c(non_pcs_sim, df2[df2$isSim != "Real","Corr"]) 
                      non_pcs_sim_all <- c(non_pcs_sim_all, df2[df2$isSim != "Real","Corr"])
                      
                    }
                    
                  }
                }
              }  
            }
            
            
            for (metric_name in names(metric_functions)) {
              
              
              func <- metric_functions[[metric_name]]
              
              
              comparisions_vector <- 
                c(`PC_REAL_PC_SIM`=func(pcs_real, pcs_sim),
                  `NPC_REAL_NPC_SIM`=func(non_pcs_real, non_pcs_sim),
                  `PC_SIM_NPC_SIM`=func(pcs_sim, non_pcs_sim),
                  `PC_REAL_NPC_REAL`=func(pcs_real, non_pcs_real))
              
              for (comp in names(comparisions_vector))  {
                
                df <- data.frame(dist=comparisions_vector[comp],
                                 metric=metric_name,
                                 session=session,
                                 group=comp)
                
                comparisions_df <- rbind(comparisions_df,
                                         df)
              }
            }
          }
          
          
          
          
          mse_comparisions_vector_all <- 
            list(`PC_REAL_PC_SIM`=(pcs_real_all -  pcs_sim_all) ** 2,
              `NPC_REAL_NPC_SIM`=(non_pcs_real_all - non_pcs_sim_all) ** 2,
              `PC_SIM_NPC_SIM`=(pcs_sim_all - non_pcs_sim_all) ** 2,
              `PC_REAL_NPC_REAL`=(pcs_real_all - non_pcs_real_all) ** 2)
          
          tmp <- data.frame(pcs_real_all,
                            pcs_sim_all,
                            non_pcs_real_all,
                            non_pcs_sim_all,
                            mice=mice_all,
                            sessions=sessions)
          
          
          mse_df <- do.call(cbind, mse_comparisions_vector_all["NPC_REAL_NPC_SIM"])
          mse_df <- as.data.frame(mse_df)
          colnames(mse_df) <- "NPC_REAL_NPC_SIM"
          mse_df$Estimated <- est_all
          
          melted <- melt(mse_df, id.vars =  "Estimated")
          colnames(melted) <- c("Estimated", "Group", "Error")
          melted_2 <- melted
          
          #melted_2$Estimated <- as.character(melted_2$Estimated)
          melted_2$Estimated[melted_2$Estimated == "1"] <- "GT = 100%"
          melted_2$Estimated[melted_2$Estimated != "GT = 100%"] <- "GT < 100%"
          melted_2$Session <- sessions
          #melted_2$Session <- melted_2$Session %% 8
          #melted_2$Session[melted_2$Session == 0] <- 8
          
          max_dist <- max(melted_2$Error)
          g_npc_error <- 
            ggplot(melted_2, aes(x=factor(Estimated), y=Error, group=1)) +
            #geom_jitter(aes(color=Estimated), alpha=0.3, size=1) + 
            #geom_boxplot(aes(fill=Estimated, color=Estimated), size=1, alpha=0.8) + 
            #geom_line(fill="gray80", color="black", width=.7, alpha=0.8, stat="summary")  + 
            geom_errorbar(width=.3, stat="summary") +
              geom_point( stat="summary") +
              geom_line(stat="summary") +
            
            #scale_fill_manual(values=c(adjustcolor(c("gray80","gray80"),)))+
            #scale_fill_manual(values=c(adjustcolor(c("#5574BC","gray80"),)))+
            #scale_color_manual(values=c(adjustcolor(c("#3F5DAB","gray65") , alpha=1))) +           
            ylab("MSE") + 
            xlab("") +
            theme(text=element_text(size=13.5),
                  plot.title = element_text(size=11)) +
            ggtitle(sprintf("%s - %s", subfield, prop)) + 
            #scale_y_continuous(expand=c(0,0)) +
            base_plot_theme
          

          
          npc_error_list <- append(npc_error_list, list(g_npc_error))

          for (metric_name in unique(comparisions_df$metric)) {
            work_df <- comparisions_df[comparisions_df$metric == metric_name,]
            work_df <- work_df[work_df$group %in% comp_groups[groups_to_use],]
            
            grp_names <- unique(work_df$group)
            pval_mt <- matrix(rep(1, times=len(grp_names) ** 2),
                              nrow=len(grp_names))
            rownames(pval_mt) <- grp_names
            colnames(pval_mt) <- grp_names
            
            work_df$group <- factor(work_df$group, levels=comp_groups[groups_to_use])

                  
              pc_df <- work_df[work_df$group == "PC_REAL_PC_SIM",]
              npc_df <- work_df[work_df$group == "NPC_REAL_NPC_SIM",]
                  
                  wx <- 
                    wilcox.test(pc_df[order(pc_df$session),"dist"],
                                npc_df[order(pc_df$session),"dist"],
                                paired=T, correct=F,
                                alternative="less")
                  
                  pval_mt[grp, grp2] <- wx$p.value
                  pval_df <- data.frame(pval=wx$p.value,
                                        statistic=wx$statistic,
                                        method=wx$method,
                                        alternative=wx$alternative,
                                        signif=signif.num(wx$p.value),
                                        subfield=subfield,
                                        metric=metric_name,
                                        conf=df_name,
                                        subsampled=subsample,
                                        prop=prop,
                                        N=nrow(npc_df),
                                        mean_NPC=mean(npc_df$dist, na.rm=T),
                                        mean_PC=mean(pc_df$dist, na.rm=T),
                                        sd_NPC=sd(npc_df$dist, na.rm=T),
                                        Sd_PC=sd(pc_df$dist, na.rm=T),
                                        sem_NPC=sem(npc_df$dist),
                                        sem_PC=sem(pc_df$dist))
                  
                  statistics_df <- rbind(statistics_df,
                                         pval_df)
                

            
            #max_dist <- max(work_df$dist)
            max_dist_err <- max(c(mean(work_df[work_df$group == "PC_REAL_PC_SIM",]$dist, na.rm=T) + sem(work_df[work_df$group == "PC",]$dist),
                              mean(work_df[work_df$group == "NPC_REAL_NPC_SIM",]$dist, na.rm=T) + sem(work_df[work_df$group == "NPC",]$dist)))
            
            
            work_df$group <- unlist(lapply(str_split(work_df$group, "_"), function(l) {l[[1]]}))
            g <- 
              ggplot(work_df, aes(x=group, y=dist)) +
              #geom_jitter(aes(color=group), alpha=0.3, size=1) + 
              #geom_boxplot(aes(fill=group, color=group), size=1, alpha=0.8) + 
              geom_bar(aes(fill=group), color=NA, width=.7, size=1, alpha=0.8, stat="summary") + 
              geom_errorbar(width=.3, stat="summary") + 
              scale_fill_manual(values=c(adjustcolor(c("#5574BC","gray80"),)))+
              scale_color_manual(values=c(adjustcolor(c("#3F5DAB","gray65") , alpha=1))) +           
              ylab(sprintf("%s", metric_name)) + 
              xlab("") +
              theme(text=element_text(color="black", size=13.5),
                    plot.title = element_text(size=11)) +
              ggtitle(sprintf("%s - %s", subfield, prop)) +
              scale_y_continuous(expand=c(0,0)) +
              base_plot_theme  
            
              if (len(groups_to_use) == 2) {
                g <- g + 
                  #ylim(c(0, max_dist * 1.25)) +
                  geom_line(data=data.frame(x=c(1,1,2,2),
                                            y=c(max_dist * 1.05,
                                                max_dist * 1.1,
                                                max_dist * 1.1,
                                                max_dist * 1.05)),
                            aes(x=x, y=y), color="black") + 
                  geom_text(x=1.5,
                            y=max_dist_err,
                            label=signif.num(pval_mt[comp_groups[groups_to_use[1]],
                                                     comp_groups[groups_to_use[2]]]))
              }

            metric_plot_lists[[metric_name]][[prop]] <- g
          }
        }
        
        
        subfield_list <- append(subfield_list,
                                list(metric_plot_lists))
      }
      
      
      
      npc_error_list$nrow <- 2
      
      gnpc_error_f <- do.call(plot_grid, npc_error_list)
      
      
      for (size_name in names(sizes)) {
        size = sizes[size_name]
        
        dir.create(sprintf("%s\\%s", 
                           write_path,
                           size_name))
        
        dir.create(sprintf("%s\\%s\\%s", 
                           write_path,
                           size_name,
                           df_name))
        
        dir.create(sprintf("%s\\%s\\%s\\subsamples_%s",
                           write_path,
                           size_name,
                           df_name, 
                           subsample))
        
        pdf(file=sprintf("%s\\%s\\%s\\subsamples_%s\\npc_error_new.pdf",
                         write_path,
                         size_name,
                         df_name, 
                         subsample),
            height=size * 2,
            width=size * 2.2)
        plot(gnpc_error_f)
        dev.off()
        
      }
      
      for (metric_name in names(metric_functions)) {
        
        final_plt_list <- list() 
      
        for (sbf in c(1:len(subfield_list))) {
          final_plt_list <- append(final_plt_list,
                                   subfield_list[[sbf]][[metric_name]])
        }
        
        final_plt_list$nrow <- 2
        gf <- do.call(arrangeGrob, final_plt_list)
        
        for (size_name in names(sizes)) {
          size = sizes[size_name]

          pdf(file=sprintf("%s\\%s\\%s\\subsamples_%s\\%s.pdf",
                           write_path,
                           size_name,
                           df_name, 
                           subsample,
                           metric_name),
              height=size * 2,
              width=size * 2.2)
          plot(gf)
          dev.off()
          
        }
      }
    }
  }
  
  write.csv(file=sprintf("%s//statistics.csv", write_path), statistics_df)  
}

figure_6_sliding_tuning <- function() {
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\sliding_tuning\\",figures_path)
  dir.create(write_path)
  
  session_groups <- list(all=1:16)
                         #first_2=c(1,2,9,10),
                         #first_4=c(1:4,9:13))
  
  
  df_res <- data.frame()
  
  metrics <- unique(sliding_tuning_df$Metric)
  subfields <- unique(sliding_tuning_df$Subfield)
  windows <- unique(sliding_tuning_df$WindowP)
  cell_funcs <- unique(sliding_tuning_df$CellFunc)
  window_funcs <- unique(sliding_tuning_df$WindowFunc)
  
  for (session_group_name in names(session_groups)) {
    sessions_to_use <- session_groups[[session_group_name]]
    
    write_path_session <- sprintf("%s\\%s",
                                  write_path, session_group_name)
    
    dir.create(write_path_session)
    for (metric in metrics) {
      
      metric_path_session <- sprintf("%s\\%s",
                                     write_path_session, metric)
      
      dir.create(metric_path_session)
      for (window in windows) {
        for (cf in c("Median", "Mean")) {
          for (wf in c("Median", "Mean")) {
            for (grp in c("Pcs", "Npcs")) {
            df_to_use <- sliding_tuning_df[sliding_tuning_df$Session %in% sessions_to_use &
                                           sliding_tuning_df$Metric == metric &
                                           sliding_tuning_df$WindowP == window &
                                           sliding_tuning_df$CellFunc == cf &
                                           sliding_tuning_df$WindowFunc == cf &
                                           sliding_tuning_df$Group == grp, ]
            
            ca1pv <- .5
            ca3pv <- .5
            if (sum(is.na(df_to_use[df_to_use$Subfield == "CA1" & df_to_use$isSim == "Real",]$Slide)) == 0) {
            ca1_test <- 
              wilcox.test(df_to_use[df_to_use$Subfield == "CA1" & df_to_use$isSim == "Real",]$Slide,
                          df_to_use[df_to_use$Subfield == "CA1" & df_to_use$isSim != "Real",]$Slide,
                          signed=T,
                          alternative="less")
            df_res  <- rbind(df_res,
                             data.frame(ca1_pv=ca1_test$p.value,
                                        sbf="CA1",
                                        win=window,
                                        metric=metric,
                                        group=session_group_name,
                                        wf=wf,
                                        cf=cf,
                                        grp=grp))
            ca1pv <- ca3_test$p.value
            
            }
            
            if (sum(is.na(df_to_use[df_to_use$Subfield == "CA1" & df_to_use$isSim == "Real",]$Slide)) == 0) {
            ca3_test <- 
              wilcox.test(df_to_use[df_to_use$Subfield == "CA3" & df_to_use$isSim == "Real",]$Slide,
                          df_to_use[df_to_use$Subfield == "CA3" & df_to_use$isSim != "Real",]$Slide,
                          signed=T,
                          alternative="less")
            df_res  <- rbind(df_res,
                             data.frame(ca1_pv=ca3_test$p.value,
                                        sbf="CA3",
                                        win=window,
                                        metric=metric,
                                        group=session_group_name,
                                        wf=wf,
                                        cf=cf,
                                        grp=grp))
            ca3pv <- ca3_test$p.value
            }
            

            g <- 
            ggplot(df_to_use) +
              geom_jitter(aes(x=Subfield, y=Slide, color=interaction(isSim, Subfield, Group)), alpha=.3,
                          position=position_jitterdodge(1)) +
              geom_boxplot(aes(x=Subfield, y=Slide, fill=interaction(isSim, Subfield, Group)),
                           alpha=.3) +
              ggtitle(sprintf("%s - %s, ca1(%s), ca3(%s)", metric, window,
                              signif.num(ca1_test$p.value),
                              signif.num(ca3_test$p.value)))
              
            
          
            
            pdf(sprintf("%s//%s_Window%s_Cell%s_%s.pdf",
                        metric_path_session,
                        window,
                        wf,
                        cf,
                        grp),
                height=3,
                width=6)
            plot(g) 
            dev.off()
            }
          }
        } 
      }
    }
  }
}

figure_6_sliding_by_time <- function() {
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\sliding_by_time\\",figures_path)
  dir.create(write_path)
  
  load(sprintf("%s\\samples\\sliding_tuning\\sliding_df_all_by_time.Rda", 
               base_path), 
       verbose=T)
  
  write_path <- sprintf("%s\\figure_6\\sliding_by_time\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\big\\",write_path))
  dir.create(sprintf("%s\\big_1\\",write_path))
  dir.create(sprintf("%s\\medium_1\\",write_path))
  dir.create(sprintf("%s\\medium_2\\",write_path))
  dir.create(sprintf("%s\\medium\\",write_path))
  dir.create(sprintf("%s\\small\\",write_path))
  
  sizes <- c(big=3,
             big_1=2.75,
             medium=2.5,
             medium_1=2,
             medium_2=1.75,
             small=1.5)
  
  sliding_tuning_df <- res$df
  median_sliding_vals <- res$median
  mean_sliding_vals <- res$mean
  
  
  session_groups <- list(#all=1:16,
                         #first_2=c(1:2,9:10),
                         #first_3=c(1:3,9:11),
                         first_4=c(1:4,9:12))
  
  
  df_res <- data.frame()
  
  metrics <- unique(sliding_tuning_df$Metric)
  subfields <- unique(sliding_tuning_df$Subfield)
  windows <- unique(sliding_tuning_df$WindowP)

  ylabs_to_use <- c(peaks="Peak movement (Bins)",
        cor="Tuning correlation (sliding)",
        cosine="Tuning cosine distance (sliding)",
        eucl="Eucledian distance (sliding)",
        wasserstein="Wasserstein distance (sliding)")
  
  
  for (session_group_name in names(session_groups)) {
    sessions_to_use <- session_groups[[session_group_name]]
  
    
    for (metric in metrics) {
      for (window in windows) {
        
        median_plt_list <- list()
        mean_plt_list <- list()
        
        med_min <- c()
        med_max <- c()
        mean_min <- c()
        mean_max <- c()
        
        for (sbf in unique(sliding_tuning_df$Subfield)) {
          for (grp in c("Pcs", "Npcs")) {
            
            ind <- sliding_tuning_df$Session %in% sessions_to_use &
              sliding_tuning_df$Metric == metric &
              sliding_tuning_df$WindowP == window &
              sliding_tuning_df$Group == grp &
              sliding_tuning_df$Subfield == sbf
            
            df_to_use <- sliding_tuning_df[ind, ]
            median_sliding_vals_to_use <- median_sliding_vals[ind]
            mean_sliding_vals_to_use <- mean_sliding_vals[ind]
            
            plot_median_res_df <- data.frame()
            plot_mean_res_df <- data.frame()
            
            for (sim_group in unique(df_to_use$isSim)) {
              
              mean_res_df <- data.frame()
              median_res_df <- data.frame()
              
              real_or_sim_indices <- which(df_to_use$isSim == sim_group)
              
              for (real_or_sim_idx in real_or_sim_indices) {

            
              median_values <- 
                rep(NA, times=40)
              mean_values <- 
                rep(NA, times=40)
              
              median_values_ind <- as.numeric(cut(seq(1.025,
                                                      20,
                                                      length.out=len(median_sliding_vals_to_use[[real_or_sim_idx]])),
                                                      breaks=seq(1,20,length.out=40)))
              
              mean_values_ind <- as.numeric(cut(seq(1.025,
                                                    20,
                                                    length.out=len(mean_sliding_vals_to_use[[real_or_sim_idx]])),
                                                  breaks=seq(1,20,length.out=40)))
              
              if (len(mean_values_ind) > 1) {
                mean_values[mean_values_ind] <- mean_sliding_vals_to_use[[real_or_sim_idx]]
              }
              
              if (len(median_values_ind) > 1) {
                median_values[median_values_ind] <- median_sliding_vals_to_use[[real_or_sim_idx]]
              }
              
              mean_res_df <- rbind(mean_res_df,
                                   mean_values)
              
              median_res_df <- rbind(median_res_df,
                                   median_values)              
              }
              
              
              colnames(mean_res_df) <- seq(0.5,20,length.out=40)
              colnames(median_res_df) <- seq(0.5,20,length.out=40)
              
              
              mean_res_sd <- apply(mean_res_df, 2, sd, na.rm=T)
              mean_res_mean <- apply(mean_res_df, 2, mean, na.rm=T)
              median_res_sd <- apply(median_res_df, 2, sd, na.rm=T)
              median_res_mean <- apply(median_res_df, 2, mean, na.rm=T)
              
              non_na_indices <- !is.na(median_res_sd)
              
              plot_median_res_df <- 
                                    rbind(plot_median_res_df,
                                      data.frame(mu=median_res_mean[non_na_indices],
                                                 sd=median_res_sd[non_na_indices],
                                                 x=as.numeric(names(median_res_sd[non_na_indices])),
                                                 group=rep(sim_group, times=sum(non_na_indices))))
              
              non_na_indices <- !is.na(mean_res_sd)
              
              plot_mean_res_df <- 
                                rbind(plot_mean_res_df,
                                  data.frame(mu=mean_res_mean[non_na_indices],
                                             sd=mean_res_sd[non_na_indices],
                                             x=as.numeric(names(mean_res_sd[non_na_indices])),
                                             group=rep(sim_group, times=sum(non_na_indices))))
            }
            

            
            gmean <- ggplot(plot_mean_res_df) +
              geom_line(aes(x=x, y=mu, color=group)) +
              geom_ribbon(aes(x=x, ymin=mu-sd, ymax=mu+sd, fill=group),
                          col=NA,
                          alpha=.2) + 
              xlab("Time (minutes)") +
              ylab(ylabs_to_use[[metric]]) +
              theme_classic() +
              ggtitle(sprintf("%s - %s - %s",
                              sbf,
                              grp,
                              window)) + 
              scale_color_manual(values=c("#652D90", "#F05A28")) + 
              scale_fill_manual(values=c("#652D90","#F05A28")) + 
              base_plot_theme +
              theme(text=element_text(size=14, color="black"),
                    axis.text = element_text(color="black"),
                    plot.title = element_text(size=9))
            
            gmed <- ggplot(plot_median_res_df) +
              geom_line(aes(x=x, y=mu, color=group)) +
              geom_ribbon(aes(x=x, ymin=mu-sd, ymax=mu+sd, fill=group),
                          col=NA,
                          alpha=.2) + 
              xlab("Time (minutes)") +
              ylab(ylabs_to_use[[metric]]) +
              theme_classic() +
              ggtitle(sprintf("%s - %s - %s",
                              sbf,
                              grp,
                              window)) + 
              scale_color_manual(values=c("#652D90", "#F05A28")) + 
              scale_fill_manual(values=c("#652D90","#F05A28")) + 
              base_plot_theme + 
              theme(text=element_text(size=14, color="black"),
                    axis.text = element_text(color="black"),
                    plot.title = element_text(size=9))
            
            mean_min <- c(mean_min, min(plot_mean_res_df$mu - plot_mean_res_df$sd))
            mean_max <- c(mean_max, max(plot_mean_res_df$mu + plot_mean_res_df$sd))
            med_min <- c(med_min,min(plot_median_res_df$mu - plot_median_res_df$sd))
            med_max <- c(med_max,  max(plot_median_res_df$mu + plot_median_res_df$sd))
                                
            mean_plt_list <- append(mean_plt_list,
                                    list(gmean))
            median_plt_list <- append(median_plt_list,
                                      list(gmed))              
          }
        }
        
        median_plt_list <- 
          lapply(median_plt_list,
                 function(plt) {
                   return(plt + ylim(c(min(med_min),
                                       max(med_max))))
                 })
        
        mean_plt_list <- 
          lapply(mean_plt_list,
                 function(plt) {
                   return(plt + ylim(c(min(mean_min),
                                       max(mean_max))))
                 })
        
        median_plt_list$nrow <- 1
        mean_plt_list$nrow <- 1
        gmed_f <- do.call(plot_grid, median_plt_list)
        gmean_f <- do.call(plot_grid, mean_plt_list)
        
        gf <- plot_grid(gmean_f, gmed_f, nrow=2)
        
        for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        write_path_session <- sprintf("%s\\%s\\%s",write_path, size_name, session_group_name)
        dir.create(write_path_session)
        metric_path_session <- sprintf("%s\\%s",write_path_session, metric)
        dir.create(metric_path_session)
          
        pdf(sprintf("%s//Window_%s.pdf",
                    metric_path_session,
                    window),
            height=size *2,
            width=size * 4)
        plot(gf) 
        dev.off()
        
        }
      }
    }
  }
}


figure_6_sliding_tuning_by_window_size <- function(load=T) {
  
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\sliding_tuning_by_window\\",figures_path)
  dir.create(write_path)
  
  
  
  if (load) {
    load(file=sprintf("%s//sliding_tuning_df.Rda",write_path),verbose=T)
  } else {
    get_sliding_tunning_corr_df(file_name = "peaks_sliding_tuning_2.Rda")
  }
  
  sessions_to_use <- 1:16
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  session_groups <- list(mean_ribbon_all_peaks_final_version=1:16)
  #first_2=c(1,2,9,10),
  #first_4=c(1:4,9:13))
  
  sizes <- c(big=3,
             big_1=2.75,
             medium=2.5,
             medium_1=2,
             medium_2=1.75,
             small=1.5)
  
  group_session = T
  
  df_res <- data.frame()
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  
  metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
  metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
  colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
  
  metrics <- unique(sliding_tuning_df$Metric)
  subfields <- unique(sliding_tuning_df$Subfield)
  windows <- unique(sliding_tuning_df$WindowP)
  cell_funcs <- unique(sliding_tuning_df$CellFunc)
  window_funcs <- unique(sliding_tuning_df$WindowFunc)
  
  all_jumps <- unique(unlist(lapply(str_split(unique(sliding_tuning_df$WindowP), "_"), function(tuple) {tuple[[2]]})))
  
  all_windows <- unique(unlist(lapply(str_split(unique(sliding_tuning_df$WindowP), "_"), function(tuple) {tuple[[1]]})))
  
  all_windows <- all_windows[3:17] # Rest of the bins are irrelevant
  
  statistics_df <- data.frame()
  
  for (session_group_name in names(session_groups)) {
    sessions_to_use <- session_groups[[session_group_name]]
    
    write_path_session <- sprintf("%s\\%s",
                                  write_path, session_group_name)
    
    #dir.create(write_path_session)
    for (metric in metrics) {
      
      metric_path_session <- sprintf("%s\\%s",
                                     write_path_session, metric)
      
      dir.create(metric_path_session)
       
      for (cf in cell_funcs) {
        for (wf in window_funcs) {
          for (jmp in all_jumps) {
            
            all_plots <- list()
            minv <- c()
            maxv <- c()
            
            
            df_to_use_outer <- 
              sliding_tuning_df[sliding_tuning_df$Session %in% sessions_to_use &
                                  sliding_tuning_df$Metric == metric &
                                  sliding_tuning_df$CellFunc == cf &
                                  sliding_tuning_df$WindowFunc == wf &
                                  grepl(jmp, sliding_tuning_df$WindowP),]
            
          

            
            pval_grps <- c()
            error_plots <- list()
            max_y_all_errors <- c()
            error_pvals <- c()
            for (sbf in subfields)  {
              errors <- c()
              
              for (grp in c("Pcs", "Npcs")) {
                # minv<- min(df_to_use_outer$Slide, na.rm=T)
                # maxv <- max(df_to_use_outer$Slide, na.rm=T)
                df_grp_sbf <- data.frame()
                pval_all <- c()
                anova_all <- data.frame()
                
                for (window in all_windows) {
                  tmp_for_test <- list()
                  
                  for (isSim in c("Real", "Simulation")) {
                    
                    df_to_use_outer$Session = df_to_use_outer$Session %% 8
                    df_to_use_outer$Session[df_to_use_outer$Session == 0] <- 8
                    search_window <- sprintf("%s_%s", window, jmp)
                    df_to_use <- df_to_use_outer[df_to_use_outer$WindowP == search_window &
                                                 df_to_use_outer$Group == grp &
                                                 df_to_use_outer$Subfield == sbf &
                                                 df_to_use_outer$isSim == isSim, ]
                    
                    print(nrow(df_to_use))
                    
                    if (nrow(df_to_use) == 0) {
                      next
                    }
                      
                    
                    if (!group_session) {
                      df_grp_sbf <- rbind(df_grp_sbf,
                                          data.frame(mu=mean(df_to_use$Sliding, na.rm=T),
                                            sd=sem(df_to_use$Sliding),#, na.rm=T),
                                            isSim=isSim,
                                            window=as.numeric(str_split(window, "W")[[1]][2])))
                    } else {
                      tmp <- 
                        ddply(df_to_use, .(isSim),  
                              function(issim_df) {
                                by_ses <- 
                                  ddply(issim_df, .(Session), 
                                        function(session_df) {
                                          return(c(mean(session_df$Sliding, na.rm=T)))
                                        })
                                return(by_ses)
                              })
                      
                      colnames(tmp) <- c("isSim", "Session", "mu")
                      
                      tmp2 <- 
                        ddply(tmp, .(isSim),  
                              function(issim_df) {
                                return(c(mean(issim_df$mu, na.rm=T),
                                         sem(issim_df$mu)))
                              })
                      
                      colnames(tmp2) <- c("isSim", "mu", "sd")
                      
                      tmp2$window <- rep(as.numeric(str_split(window, "W")[[1]][2]), 
                                         nrow(tmp2))
                      
                      df_grp_sbf <- rbind(df_grp_sbf,
                                          tmp2)
                      
                      tmp_for_test <- append(tmp_for_test,
                                            list(tmp$mu))
                      
                      tmp$window <-  rep(as.numeric(str_split(window, "W")[[1]][2]), 
                                         nrow(tmp))
                      anova_all <- rbind(anova_all,
                                         tmp)
                    }
                  }
                  
                  tmp3 <- do.call(cbind, tmp_for_test)
                  pval_all <- c(pval_all,
                                wilcox.test(tmp3[,1], tmp3[,2])$p.value)
                }
                
                df_grp_sbf <- as.data.frame(df_grp_sbf)
                
                
                minv <- c(minv, min(df_grp_sbf$mu - df_grp_sbf$sd))
                maxv <- c(maxv, max(df_grp_sbf$mu + df_grp_sbf$sd))
                
                
                g <- 
                  ggplot(df_grp_sbf, aes(x=window, y=mu, group=isSim, color=isSim, fill=isSim)) +
                  geom_line() +
                  #geom_point(size=2) +
                  #geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd))+
                  geom_ribbon(aes(ymin=mu-sd, ymax=mu+sd), color=NA, alpha=.3) + 
                  theme_classic() +
                  ggtitle(sprintf("%s - %s",
                                  sbf,
                                  grp)) + 
                  scale_color_manual(values=c("#652D90", "#F05A28"), breaks=c("Real", "Simulation")) + 
                  scale_fill_manual(values=c("#652D90","#F05A28"), breaks=c("Real", "Simulation")) + 
                  base_plot_theme +
                  theme(text=element_text(size=14, color="black"),
                        axis.text = element_text(color="black"),
                        plot.title = element_text(size=9))  +
                  xlab("Window size (frames)") +
                  ylab("Peak shift (cm)")
                  #ylab(metric) # + ylim(c(minv, maxv))
                  
                
                all_plots <- append(all_plots,
                                    list(g))
                
                pval_grps <- rbind(pval_grps, pval_all)
                
                error <- df_grp_sbf[df_grp_sbf$isSim == "Real",][order(df_grp_sbf[df_grp_sbf$isSim == "Real",]$window),"mu"] - 
                  df_grp_sbf[df_grp_sbf$isSim != "Real",][order(df_grp_sbf[df_grp_sbf$isSim != "Real",]$window),"mu"] 
                errors <- rbind(errors,
                                error)
              }
              
              
              errors <- t(errors ** 2)
               colnames(errors) <- c("PC", "NPC")
              melted_errors_df <- melt(errors)
              colnames(melted_errors_df) <- c("n", "Group", "MSE")
              max_y <- max(colMeans(errors) + apply(errors, 2, sem))
              wilk_t <-  wilcox.test(errors[,1], errors[,2], alternative="less", paired=T, corrected=F)
              error_pvals <- c(error_pvals,
                               wilk_t$p.value)
              
              statistics_df <- rbind(statistics_df,
                                     data.frame(
                                       subfield=sbf,
                                       pval=wilk_t$p.value,
                                       method=wilk_t$method,
                                       statistic=wilk_t$statistic,
                                       alternative=wilk_t$alternative,
                                       signif=signif.num((wilk_t$p.value)),
                                       jump=jmp,
                                       cf=cf,
                                       wf=wf,
                                       metric=metric,
                                       mean_pc=apply(errors, 2, mean)[1],
                                       mean_npc=apply(errors, 2, mean)[2],
                                       sd_pc=apply(errors, 2, sd)[1],
                                       sd_npc=apply(errors, 2, sd)[2],
                                       sem_pc=apply(errors, 2, sem)[1],
                                       sem_npc=apply(errors, 2, sem)[2]))
              
              gerror <-
                ggplot(melted_errors_df, aes(x=Group, y=MSE)) +
                  geom_bar(stat="summary", width=.7) +
                  geom_errorbar(stat="summary", width=.3) +
                  #geom_jitter(aes(color=Group),position=position_jitterdodge(.5))
                  base_plot_theme + 
                  scale_y_continuous(expand=c(0,0)) +
                  ggtitle(sbf) + 
                  theme(text=element_text(size=14, color="black"),
                        axis.text = element_text(color="black"),
                        plot.title = element_text(size=9)) + 
                  xlab("")
              
                  max_y_all_errors <- c(max_y_all_errors, max_y + .25)
                  
                  error_plots <- append(error_plots,
                                        list(gerror))
             }
            
            all_plots <- lapply(all_plots, 
                                function(p) {p + ylim(c(min(minv), max(maxv)))})
            all_plots$nrow <- 1
            gf <- do.call(plot_grid, all_plots)
            
            max_y <- max(max_y_all_errors)
            error_plots <- lapply(1:len(error_plots),
                                 function(p_idx) {
                                      p <- error_plots[[p_idx]]
                                      p <- p + 
                                       geom_line(data=data.frame(x=c(1,1,2,2),
                                                                 y=c(max_y  + .15, max_y + .25, max_y  +.25, max_y + .15)),
                                                 aes(x=x,y=y)) + 
                                       geom_text(label=signif.num(error_pvals[p_idx]),
                                                 x=1.5,
                                                 y=max_y + .15)
                                      
                                      return(p)
                                     })
            error_plots$nrow <- 2
            gf_error <- do.call(plot_grid, error_plots)
            for (size_name in names(sizes)) {
              size = sizes[[size_name]]
              write_path_session <- sprintf("%s\\%s\\",write_path, size_name)
              dir.create(write_path_session)
              write_path_session <- sprintf("%s\\%s\\%s",write_path, size_name, session_group_name)
              dir.create(write_path_session)
              metric_path_session <- sprintf("%s\\%s",write_path_session, metric)
              dir.create(metric_path_session)
              window_function_path_session <- sprintf("%s\\window_%s", metric_path_session,wf)
              dir.create(window_function_path_session)
              
              pdf(sprintf("%s//Jump_%s_cell_%s.pdf",
                          window_function_path_session,
                          jmp,
                          cf),
                  height=size *1,
                  width=size * 3.6)
              plot(gf) 
              dev.off()
              
              pdf(sprintf("%s//Jump_%s_cell_%s_error.pdf",
                          window_function_path_session,
                          jmp,
                          cf),
                  height=size *2,
                  width=size * .65)
              plot(gf_error) 
              dev.off()
              
            }
          }          
        } 
      }
    }
  }
  
  write.csv(file=sprintf("%s//statistics.csv",write_path), statistics_df)
  save(file=sprintf("%s//sliding_tuning_df.Rda",write_path), sliding_tuning_df)
}

figure_6_prob_fire_df <- function() {
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\probability_firing_near_analysis\\",figures_path)
  dir.create(write_path)
  
  missing_paths <- c()
  statistics_df <- data.frame()
  wilcox_statistics_df <- data.frame()
  
  for (subfield_name in c("CA1")){#}, "CA3")) {
    if (subfield_name == "CA1") {
      data_paths_ind <- 1:8
    } else {
      data_paths_ind <- 9:18
    }
  sessions_to_use <- 1:16
  prob_fire_df <- data.frame()
  consec_df <- list()
  file_name="new_probability_to_fire_by_pos.Rda"
  true_data_paths <- sprintf("%s\\%s", all_data_paths[data_paths_ind], "equalized")
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  
  metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
  metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
  colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
  
  all_pos_mat <- list()
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.3,
          small=2)
  
  for (group_name in c("cyclic_pcs", "cyclic_non_pcs")) {
    for (isSim in c("Real", "Simulation")) {
      all_pos_mat[[sprintf("%s_%s", group_name, isSim)]]  <- 
        matrix(rep(0, times = 24**2), nrow=24)
    }
  }
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    lidx = 0
    for (ext in c("", "\\fit_simulation")) {
      
      #print (len(session_paths))
      # Run through all paths
      for (tpath in session_paths) {
        
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        

        
        lidx <- idx
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          #print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=F)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
        
        
        
        isSimu <- ifelse(ext == "", "Real", "Simulation")
        
        
        for (group_name in c("cyclic_pcs", "cyclic_non_pcs")) {
          # fdf_1 <- apply(final_result[[group_name]]$pos_diff[,1:2], 1, sum)
          # fdf_f <- apply(final_result[[group_name]]$pos_diff[,5:7], 1, sum)
          # fdf <- cbind(fdf_1, final_result[[group_name]]$pos_diff[,3:4])
          # fdf <- cbind(fdf, fdf_f)
          
          if (group_name == "cyclic_pcs") {
            group_ind_tu <- final_result$place_cells_metadata$ind[final_result$place_cells_metadata$cyclic_pcs]
          } else {
            group_ind_tu <- final_result$place_cells_metadata$ind[final_result$place_cells_metadata$cyclic_non_pcs]
          }
          

          fdf <- final_result[[group_name]]$pos_diff
          fdf2 <- final_result$trav_prob[group_ind_tu,]
          mp <- c(apply(final_result$trav_tuning_curve[group_ind_tu,], 1, function(r) {diff(r, lag=5)}))
          Consec <- abs(mp[!is.na(mp)])
          #if (len(mp) == 0) {next}
          df <-  c(colSums(fdf) / sum(fdf),
                   ByTrav=colSums(fdf2) / sum(fdf2),
                   ConsecRunTun=mean(Consec),
                   ConsecRunTunMed=median(Consec),
                   Session=idx,
                   Subfield=subfield,
                   Mice=mice_str,
                   isSim=isSimu,
                   Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                   Group=group_name)
          
          c_df <- list(Shifts=Consec)
          c_df$isSim <- isSimu
          c_df$Group <- group_name
          c_df$Session = idx#ifelse(idx %% 8 == 0, 8, idx %% 8)
          consec_df <- append(consec_df,
                             list(c_df))
          
          prob_fire_df <- rbind(df, prob_fire_df)
          
          all_pos_mat[[sprintf("%s_%s", group_name, isSimu)]] <- 
            all_pos_mat[[sprintf("%s_%s", group_name, isSimu)]] + 
            final_result[[group_name]]$pos_matrix
        }
      }
    } 
  }
  
  prob_df <- apply(prob_fire_df[1:17], 2, function(col) {as.numeric(unlist(col))})
  
  final_df <- cbind(prob_df, prob_fire_df[,18:22])
  
  
  probability_belt_range <- c(2,3,6,10,15,20,25)
  colnames(final_df) <- c(paste("by_", probability_belt_range, sep=""),
                          paste("trav_by_", probability_belt_range, sep=""),
                          "MeanTunPeak", "MedTunPeak",
                          "Session", "Subfield", "Mouse", "isSim", "Direction", "Group")
  

  
  
    final_df_f <- final_df#[final_df$Subfield == subfield_name,]
    final_df_f$AbsSession <- final_df_f$Session %% 8
    final_df_f$AbsSession[final_df_f$AbsSession ==0] <- 8
    
    sim_prob <- final_df_f[final_df_f$Group != "cyclic_pcs" &
                             final_df_f$isSim == "Simulation",]
    
    real_prob <- final_df_f[final_df_f$Group != "cyclic_pcs" &
                              final_df_f$isSim != "Simulation",]
    
    sim_prob <- sim_prob[order(paste(paste(sim_prob$Session, sim_prob$Mouse, sep="_"), sim_prob$Direction)),]
    real_prob <- real_prob[order(paste(paste(real_prob$Session, real_prob$Mouse, sep="_"), real_prob$Direction)),]
    
    wilc <- wilcox.test(sim_prob$trav_by_2, real_prob$trav_by_2, alternative="less", paired=T)
    
    wilcox_pvals_df <- data.frame(Df=nrow(sim_prob),
                                  `Wvalue`=wilc$statistic,
                                  `Pval`=wilc$p.value,
                                  `Group`=subfield_name,
                                  `Method`=wilc$method,
                                  `Alternative`=wilc$alternative,
                                  `Signif.code`=signif.num(wilc$p.value))
    
    wilcox_statistics_df  <- rbind(wilcox_statistics_df,
                                   wilcox_pvals_df)
    
    gprob_near <- 
      ggplot(final_df_f) +#[final_df_f$Group != "cyclic_non_pcs",]) + 
      geom_line(aes(x=AbsSession, y=MeanTunPeak , color=isSim, 
                       group=interaction(isSim, Group)),
                stat="summary") +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="top",
            text=element_text(size=15),
            plot.title=element_text(size=11)) +
      xlab("") +
      ylab("P(Next event <= 2)") +
      ggtitle(subfield_name)
    
    pcs_delta_df <- data.frame()
    delta_df <- data.frame()
    paired_vec_pcs <-  c()
    paired_vec_npcs <- c()
    
  for(ses in unique(final_df_f$AbsSession)) {
      for (mice in unique(final_df_f$Mouse)) {
        for (direc in unique(final_df_f$Direction)) {
          
          comp_df <- final_df_f[final_df_f$AbsSession == ses & 
                                  final_df_f$Direction == direc & 
                                  final_df_f$Mouse == mice & 
                                  final_df_f$Group == "cyclic_non_pcs",] 
          
          pcs_comp_df <-  final_df_f[final_df_f$Session == ses & 
                                       final_df_f$Direction == direc & 
                                       final_df_f$Mouse == mice & 
                                       final_df_f$Group == "cyclic_pcs",] 
          
          metadata_inst <- metadata_df[metadata_df$Session == ses & 
                                         metadata_df$Direction == direc & 
                                         metadata_df$Mice == mice,]
          
          
          paired_vec_npcs <- rbind(paired_vec_npcs,
                                   c(comp_df[comp_df$isSim == "Real","by_2"], 
                                     comp_df[comp_df$isSim == "Simulation","by_2"]))
          
          paired_vec_pcs <- rbind(paired_vec_pcs,
                                  c(pcs_comp_df[pcs_comp_df$isSim == "Real","by_2"], 
                                    pcs_comp_df[pcs_comp_df$isSim == "Simulation","by_2"]))
          
          pcs_delta <- pcs_comp_df[pcs_comp_df$isSim == "Real","by_2"] - 
            pcs_comp_df[pcs_comp_df$isSim == "Simulation","by_2"]
          pcs_delta_df <- rbind(pcs_delta_df,
                                data.frame(Delta=pcs_delta ** 1,
                                           Mice=mice,
                                           Sessions=ses,
                                           Direction=direc,
                                           Est=metadata_inst$Est))
          
          delta <- comp_df[comp_df$isSim == "Real","by_2"] - 
            comp_df[comp_df$isSim == "Simulation","by_2"]
          
          delta_df <- rbind(delta_df,
                            data.frame(Delta=delta ** 1,
                                       Mice=mice,
                                       Sessions=ses,
                                       Direction=direc,
                                       Est=metadata_inst$Est))
          
        }
      }
    }
    
    paired_vec_pcs <- cbind(paired_vec_pcs, pcs_delta_df$Sessions)
    paired_vec_npcs <- cbind(paired_vec_npcs, delta_df$Sessions)
    colnames(paired_vec_npcs) <- c("Real", "Data", "Session")
    colnames(paired_vec_pcs) <- c("Real", "Data", "Session")
    
    paired_vec_pcs <- as.data.frame(paired_vec_pcs)
    paired_vec_npcs <- as.data.frame(paired_vec_npcs)
    
    paired_vec_pcs$fSession <- paired_vec_pcs$Session %% 8
    paired_vec_pcs[paired_vec_pcs$fSession == 0,]$fSession <- 8
    
    paired_vec_npcs$fSession <- paired_vec_npcs$Session %% 8
    paired_vec_npcs[paired_vec_npcs$fSession == 0,]$fSession <- 8
    
    npcs_over_session <- ddply(paired_vec_npcs, .(Session), function(grp) {mean((grp[,1] - grp[,2]))})
    pcs_over_session <- ddply(paired_vec_pcs, .(Session), function(grp) {mean((grp[,1] - grp[,2]))})
    
    # delta_df$fsession <- delta_df$Sessions %% 8
    # delta_df[delta_df$fsession == 0,]$fsession <- 8
    
    all_delta_df <- rbind(pcs_delta_df, delta_df)
    g_by_est <- 
      ggplot(delta_df, aes(x=Sessions, y=Delta)) + 
      geom_line(stat="summary", size=1) +
      geom_errorbar(stat="summary", size=0.75, width=.05) +
      geom_point(stat="summary", size=2,) + 
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="NA",
            text=element_text(size=15),
            plot.title=element_text(size=11)) + 
      scale_color_manual(values=c(CA1="#CE3736", 
                                  CA3="#3F5DAB")) +
      xlab("Session") +
      ylab("P(Fire near| Real) - P(Fire near | Sim)")
    
    delta_df$Est <- factor(delta_df$Est)
    
    model.aov <- aov(data=delta_df,
                     formula=Delta ~ 
                       Est + 
                       Error(factor(Mice)))
    
    pvals <- summary(model.aov)
    
    pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
    pvals_df$Group <- subfield_name
    pvals_df$Measurement <- "Estimated"
    pvals_df$Factor <- rownames(pvals_df)
    rownames(pvals_df) <- c()
    statistics_df <- rbind(statistics_df, pvals_df)
    mprb = 0.02
    prob <- lapply(all_pos_mat, function(m) {m / sum(m)})
    #prob_ca3 <- lapply(all_pos_mat_ca3, function(m) {m / sum(m)})
    rpc <- pheatmap(prob$cyclic_pcs_Real, cluster_rows=F, cluster_cols=F, breaks=seq(0,mprb, length.out=100), border_col=NA)
    spc <- pheatmap(prob$cyclic_pcs_Simulation, cluster_rows=F, cluster_cols=F, breaks=seq(0,mprb, length.out=100), border_col=NA)
    rnpc <- pheatmap(prob$cyclic_non_pcs_Real, cluster_rows=F, cluster_cols=F, breaks=seq(0,mprb, length.out=100), border_col=NA)
    snpc <- pheatmap(prob$cyclic_non_pcs_Simulation, cluster_rows=F, cluster_cols=F, breaks=seq(0,mprb, length.out=100), border_col=NA)
    ls <- list(rpc[[4]], spc[[4]], rnpc[[4]], snpc[[4]])
    ls$nrow <- 2
    prob_plot <- do.call(plot_grid, ls)
    pf <- plot_grid(prob_plot, gprob_near, g_by_est, nrow=1)
    
    
    for (size_name in names(sizes)) {
      size = sizes[size_name]
      dir.create(sprintf("%s\\%s", 
                         write_path,
                         size_name))
      
      pdf(file=sprintf("%s\\%s\\prob_fire_analysis_%s_by_session.pdf",
                       write_path,
                       size_name,
                       subfield_name),
          height=size * 1,
          width=size * 3)
      
      
      plot(pf)
      dev.off()
      
    }
  }
  
  statistics_df$signif <- signif.num(statistics_df$`Pr(>F)`)
  write.csv(file=sprintf("%s\\wilcox_stats.csv",write_path),
            wilcox_statistics_df)
  write.csv(file=sprintf("%s\\anova_stats.csv",write_path),
            statistics_df)
}

get_sliding_tunning_corr_df <- function(true_data_paths, 
                                        ext="",
                                        verbose=F,
                                        output_name="",
                                        file_name="sliding_tuning.Rda") {
  
  
  sessions_to_use <- 1:16
  sliding_tuning_df <- data.frame()
  median_sliding_vals <- list()
  mean_sliding_vals <- list()
  
  values_list <- list()
  
  BY_TIME = F
  ALL_CELLS = F
  INFO_THRESHOLD = T
  
  if (BY_TIME) {
    cell_func <- list(Median=median)
    window_func <- list(Median=median)
  } else { 
    cell_func <- list(Median_new=median)
    window_func <- list(Median_new=median)
  }
  
  
  
  
  for (working_path in true_data_paths[1:8]) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    for (ext in c("", "\\fit_simulation")) {
      
      # Run through all paths
      for (tpath in session_paths) {
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=F)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
        
        
        place_cells_md <- final_result$place_cells_metadata
        sim_df <- final_result$similarity
        npc_ind <- place_cells_md$ind[place_cells_md$cyclic_non_pcs]
        pc_ind <- place_cells_md$ind[place_cells_md$cyclic_pcs]
        
        if (INFO_THRESHOLD) {
          
          
          filtered_ind <- order(place_cells_md$SI[[1]][place_cells_md$cyclic_non_pcs], decreasing = F)[1:ceiling(len(place_cells_md$cyclic_non_pcs) * 0.2)]
          #filtered_ind <- filtered_ind[order(place_cells_md$FR[filtered_ind], decreasing=T)[1:10]]
          npc_ind <- place_cells_md$ind[place_cells_md$cyclic_non_pcs[filtered_ind]]
          
          filtered_ind <- order(place_cells_md$SI[[1]][place_cells_md$cyclic_pcs], decreasing = F)[1:ceiling(len(place_cells_md$cyclic_pcs) * 0.2)]
          #filtered_ind <- filtered_ind[order(place_cells_md$FR[filtered_ind], decreasing=T)[1:10]]
          pc_ind <- place_cells_md$ind[place_cells_md$cyclic_pcs[filtered_ind]]
        }
        
        if (ALL_CELLS) {
          group_indices <- list(ALL=place_cells_md$ind) 
        } else {
          group_indices <- list(Pcs=pc_ind, Npcs=npc_ind)  
        }
        
        
        for (window_params in names(sim_df)) {
          for (metric in names(sim_df[[1]])) {
            for (group_name in names(group_indices)) {
              
              ind <- group_indices[[group_name]]
              
              #for (window_f in names(window_func)) {
                #for (cell_f in names(cell_func)) {
                  comp <- as.matrix(sim_df[[window_params]][[metric]][ind,])
                  
                  if (BY_TIME) {
                    
                    if (metric == "peaks") {
                      median_val <- apply(abs(comp), 2, median, na.rm=T)
                      mean_val <- apply(abs(comp), 2, mean, na.rm=T)
                    } else {
                      median_val <- apply(comp, 2, median, na.rm=T)
                      mean_val <- apply(comp, 2, mean, na.rm=T)
                    }
                    
                    df <-  data.frame(Session=idx,
                                      Subfield=subfield,
                                      Mice=mice_str,
                                      Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                                      isSim=sim,
                                      Metric=metric,
                                      WindowP=window_params,
                                      Group=group_name)
                    
                    sliding_tuning_df <- rbind(sliding_tuning_df,
                                               df) 
                    
                    mean_sliding_vals <- append(mean_sliding_vals,
                                                list(mean_val))
                    
                    median_sliding_vals <- append(median_sliding_vals,
                                                  list(median_val))
                    
                  } else {
                    
                    if (metric == "peaks") {
                      #val <- 
                      #cell_func[[cell_f]](apply(abs(comp), 1, window_func[[window_f]], na.rm=T), na.rm=T)
                      val <- median(abs(comp), na.rm=T)
                    } else {
                      val <- 
                        cell_func[[cell_f]](apply(comp, 1, window_func[[window_f]], na.rm=T), na.rm=T)
                    }
                    
                    values_list <- append(values_list, list(c(comp)))
                    
                    df <-  data.frame(Sliding=val,
                                      Session=idx,
                                      Subfield=subfield,
                                      Mice=mice_str,
                                      Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                                      isSim=sim,
                                      Metric=metric,
                                      WindowP=window_params,
                                      #CellFunc=cell_f,
                                      #WindowFunc=window_f,
                                      Group=group_name)
                    
                    sliding_tuning_df <- rbind(sliding_tuning_df,
                                               df)
                    
                  }
                #}
              #} 
            }
          }
        }
      }
    } 
  }
  
  if (BY_TIME) {
    res = list(df=sliding_tuning_df,
               median=median_sliding_vals,
               mean=mean_sliding_vals)
    
    save(file=sprintf("%s\\%s\\%s\\%s", base_path, "samples", "sliding_tuning", "sliding_df_all_by_time.Rda"), res)
    return(res)
    
  } else {
    save(file=sprintf("%s\\%s\\%s\\%s", base_path, "samples", "sliding_tuning", "liron_sliding.Rda"), sliding_tuning_df)
    return(sliding_tuning_df)
  }
}










new_heatmap_peaks_by_trav <- function(true_data_paths, 
                                        ext="",
                                        verbose=F,
                                        output_name="",
                                        file_name="peak_movement_by_runs.Rda") {
  
  
  dir.create(sprintf("%s\\figure_6\\",figures_path))
  write_path <- sprintf("%s\\figure_6\\new_heatmap_peaks_by_trav",figures_path)
  dir.create(write_path)
  
  sizes <- c(big=3,
             big_1=2.75,
             medium=2.5,
             medium_1=2,
             medium_2=1.75,
             small=1.5)
  
  sessions_to_use <- 1:16
  maxruns <- 15
  metadata_df <- data.frame()
  
  pc_matrices <- list()
  npc_matrices <- list()
  
  
  
  for (working_path in true_data_paths[1:8]) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    for (ext in c("", "\\fit_simulation")) {
      
      # Run through all paths
      for (tpath in session_paths) {
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=T)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
        
        
        place_cells_md <- final_result$place_cells_md
      
        npc_ind <- place_cells_md$ind[place_cells_md$cyclic_non_pcs]
        pc_ind <- place_cells_md$ind[place_cells_md$cyclic_pcs]
        
        
        pc_df <- do.call(rbind, final_result$sliding_df[final_result$used_neur_indices %in% pc_ind])
        npc_df <- do.call(rbind, final_result$sliding_df[final_result$used_neur_indices %in% npc_ind])
        
        pc_df <- as.data.frame(pc_df)
        npc_df <- as.data.frame(npc_df)
        
        pc_mt <- matrix(nrow=maxruns, ncol=24)
        npc_mt <- matrix(nrow=maxruns, ncol=24)
        
        for (i in 1:maxruns) { 
          for (j in 1:24) {
            pc_mt[i,j]<- nrow(pc_df[pc_df$TravDiff == i & pc_df$PeakDiff == j,])
            npc_mt[i,j]<- nrow(npc_df[npc_df$TravDiff == i & npc_df$PeakDiff == j,])
          }
        }
        
        pc_mt <- pc_mt/sum(pc_mt)
        npc_mt <- npc_mt/sum(npc_mt)
        
        pc_matrices <- append(pc_matrices, list(pc_mt))
        npc_matrices <- append(npc_matrices, list(npc_mt))
        
        df <-
               data.frame(Session=idx,
                          Subfield=subfield,
                          Mice=mice_str,
                          Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                          isSim=sim)
        
        
        metadata_df <- rbind(metadata_df, df)
        
        
      }
    } 
  }
  
  metadata_df$FixedSession = metadata_df$Session %% 8 
  metadata_df$FixedSession[metadata_df$FixedSession == 0] <- 8
  
    
    real_npcs_all <- list()
    real_pcs_all <- list()
    sim_npcs_all <- list()
    sim_pcs_all <- list()
    
    real_npcs_plot_all <- list()
    real_pcs_plot_all <- list()
    sim_npcs_plot_all <- list()
    sim_pcs_plot_all <- list()
    
    brks = seq(0, .3, length.out=100)
    rdylbu_col = rdylbu_cg(100)
    
    minmaxscale <- function(r) {(r-min(r)) / (max(r) - min(r))}
    
    
    for (ses in 1:8) { 
      
      real_mts <- which(metadata_df$FixedSession == ses & metadata_df$isSim == "Real")
      
      real_final_pc_mt <- pc_matrices[[real_mts[1]]]
      real_final_npc_mt <- npc_matrices[[real_mts[1]]]
      
      sim_mts <- which(metadata_df$Session == ses & metadata_df$isSim == "Real")
      
      sim_final_pc_mt <- pc_matrices[[sim_mts[1]]]
      sim_final_npc_mt <- npc_matrices[[sim_mts[1]]]
      
      for (idx in real_mts[-1]) {
        print(idx)
        real_final_pc_mt <- real_final_pc_mt + pc_matrices[[idx]]
        real_final_npc_mt <- real_final_npc_mt + npc_matrices[[idx]]
      }
      
      fidx <- 1:maxruns
      for (idx in sim_mts[-1]) {
        print(idx)
        sim_final_pc_mt <- sim_final_pc_mt + pc_matrices[[idx]]
        sim_final_npc_mt <- sim_final_npc_mt + npc_matrices[[idx]]
      }
      
      
      real_final_npc_mt <- real_final_npc_mt[fidx,]
      real_final_pc_mt <- real_final_pc_mt[fidx,]
      sim_final_npc_mt <- sim_final_npc_mt[fidx,]
      sim_final_pc_mt <- sim_final_pc_mt[fidx,]
      
      
      # real_final_npc_mt <- t(apply(real_final_npc_mt, 1, function(r) {r/sum(r)}))
      # real_final_pc_mt <- t(apply(real_final_pc_mt, 1, function(r) {r/sum(r)}))
      # sim_final_npc_mt <- t(apply(sim_final_npc_mt, 1, function(r) {r/sum(r)}))
      # sim_final_pc_mt <- t(apply(sim_final_pc_mt, 1, function(r) {r/sum(r)}))
      
      
      
      real_pc_plot <- ph(real_final_pc_mt, border_col=NA, legend=F, col=rdylbu_col, breaks=brks)[[4]]
      real_npc_plot <- ph(real_final_npc_mt, border_col=NA, legend=F, col=rdylbu_col, breaks=brks)[[4]]
      real_pcs_plot_all <- append(real_pcs_plot_all, list(real_pc_plot))  
      real_npcs_plot_all <- append(real_npcs_plot_all, list(real_npc_plot)) 
      
      
      sim_pc_plot <- ph(sim_final_pc_mt, border_col=NA, legend=F, col=rdylbu_col, breaks=brks)[[4]]
      sim_npc_plot <- ph(sim_final_npc_mt, border_col=NA, legend=F, col=rdylbu_col, breaks=brks)[[4]]
      sim_pcs_plot_all <- append(sim_pcs_plot_all, list(sim_pc_plot))  
      sim_npcs_plot_all <- append(sim_npcs_plot_all, list(sim_npc_plot)) 
      
      real_npcs_all <- append(real_npcs_all, list(real_final_npc_mt))
      real_pcs_all <- append(real_pcs_all, list(real_final_pc_mt))
      sim_npcs_all <- append(sim_npcs_all, list(sim_final_npc_mt))
      sim_pcs_all <- append(sim_pcs_all, list(sim_final_pc_mt))
      
    }
    
    real_npcs_plot_all$nrow <- 1
    real_pcs_plot_all$nrow <- 1
    sim_npcs_plot_all$nrow <- 1
    sim_pcs_plot_all$nrow <- 1
    
    real_pc_p <- do.call(plot_grid, real_pcs_plot_all)
    real_npc_p <- do.call(plot_grid, real_npcs_plot_all)
    sim_pc_p <- do.call(plot_grid, sim_pcs_plot_all)
    sim_npc_p <- do.call(plot_grid, sim_npcs_plot_all)
        
    pf <- plot_grid(real_pc_p,
                    real_npc_p,
                    sim_pc_p,
                    sim_npc_p,
      nrow=4)
    
    
      tmp <- lapply(1:8, function(i) {real_npcs_all[[i]] - sim_npcs_all[[i]]})
    
    diff_df <- data.frame(x=1:8, y=unlist(lapply(tmp, function(rt) {mean(rt[1:2,1:2])})))
    
    #gdiff <- 
    ggplot(data=diff_df) + 
      geom_line(aes(x=x,y=y)) +
      base_plot_theme +
      xlab("Session") +
      ylab("Difference (Sim - Real)") + 
    theme(axis.text=element_text(color="black"))
    
    for (size_name in names(sizes)) {
      size = sizes[size_name]
      dir.create(sprintf("%s\\%s", 
                         write_path,
                         size_name))
      
      pdf(file=sprintf("%s\\%s\\peak_diff_trav_diff_ALL.pdf",
                       write_path,
                       size_name),
          height=size * 2,
          width=size * 4)
      
      
      plot(pf)
      dev.off()
      
      pdf(file=sprintf("%s\\%s\\diff_quant.pdf",
                       write_path,
                       size_name),
          height=size,
          width=size)
      
      
      plot(gdiff)
      dev.off()
      
    }
  
}

get_sliding_tunning_corr_df <- function(true_data_paths, 
                                        ext="",
                                        verbose=F,
                                        output_name="",
                                        file_name="sliding_tuning.Rda") {
  
  
  sessions_to_use <- 1:16
  sliding_tuning_df <- data.frame()
  median_sliding_vals <- list()
  mean_sliding_vals <- list()
  
  values_list <- list()
  
  BY_TIME = F
  ALL_CELLS = F
  INFO_THRESHOLD = T
  
  if (BY_TIME) {
    cell_func <- list(Median=median)
    window_func <- list(Median=median)
  } else { 
    cell_func <- list(Median_new=median)
    window_func <- list(Median_new=median)
  }
  
  
  
  
  for (working_path in true_data_paths[1:8]) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    for (ext in c("", "\\fit_simulation")) {
      
      # Run through all paths
      for (tpath in session_paths) {
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=F)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
        
        
        place_cells_md <- final_result$place_cells_metadata
        sim_df <- final_result$similarity
        npc_ind <- place_cells_md$ind[place_cells_md$cyclic_non_pcs]
        pc_ind <- place_cells_md$ind[place_cells_md$cyclic_pcs]
        
        if (INFO_THRESHOLD) {
          
          
          filtered_ind <- order(place_cells_md$SI[[1]][place_cells_md$cyclic_non_pcs], decreasing = F)[1:ceiling(len(place_cells_md$cyclic_non_pcs) * 0.2)]
          #filtered_ind <- filtered_ind[order(place_cells_md$FR[filtered_ind], decreasing=T)[1:10]]
          npc_ind <- place_cells_md$ind[place_cells_md$cyclic_non_pcs[filtered_ind]]
          
          filtered_ind <- order(place_cells_md$SI[[1]][place_cells_md$cyclic_pcs], decreasing = F)[1:ceiling(len(place_cells_md$cyclic_pcs) * 0.2)]
          #filtered_ind <- filtered_ind[order(place_cells_md$FR[filtered_ind], decreasing=T)[1:10]]
          pc_ind <- place_cells_md$ind[place_cells_md$cyclic_pcs[filtered_ind]]
        }
        
        if (ALL_CELLS) {
          group_indices <- list(ALL=place_cells_md$ind) 
        } else {
          group_indices <- list(Pcs=pc_ind, Npcs=npc_ind)  
        }
        
        
        for (window_params in names(sim_df)) {
          for (metric in names(sim_df[[1]])) {
            for (group_name in names(group_indices)) {
              
              ind <- group_indices[[group_name]]
              
              #for (window_f in names(window_func)) {
              #for (cell_f in names(cell_func)) {
              comp <- as.matrix(sim_df[[window_params]][[metric]][ind,])
              
              if (BY_TIME) {
                
                if (metric == "peaks") {
                  median_val <- apply(abs(comp), 2, median, na.rm=T)
                  mean_val <- apply(abs(comp), 2, mean, na.rm=T)
                } else {
                  median_val <- apply(comp, 2, median, na.rm=T)
                  mean_val <- apply(comp, 2, mean, na.rm=T)
                }
                
                df <-  data.frame(Session=idx,
                                  Subfield=subfield,
                                  Mice=mice_str,
                                  Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                                  isSim=sim,
                                  Metric=metric,
                                  WindowP=window_params,
                                  Group=group_name)
                
                sliding_tuning_df <- rbind(sliding_tuning_df,
                                           df) 
                
                mean_sliding_vals <- append(mean_sliding_vals,
                                            list(mean_val))
                
                median_sliding_vals <- append(median_sliding_vals,
                                              list(median_val))
                
              } else {
                
                if (metric == "peaks") {
                  #val <- 
                  #cell_func[[cell_f]](apply(abs(comp), 1, window_func[[window_f]], na.rm=T), na.rm=T)
                  val <- median(abs(comp), na.rm=T)
                } else {
                  val <- 
                    cell_func[[cell_f]](apply(comp, 1, window_func[[window_f]], na.rm=T), na.rm=T)
                }
                
                values_list <- append(values_list, list(c(comp)))
                
                df <-  data.frame(Sliding=val,
                                  Session=idx,
                                  Subfield=subfield,
                                  Mice=mice_str,
                                  Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                                  isSim=sim,
                                  Metric=metric,
                                  WindowP=window_params,
                                  #CellFunc=cell_f,
                                  #WindowFunc=window_f,
                                  Group=group_name)
                
                sliding_tuning_df <- rbind(sliding_tuning_df,
                                           df)
                
              }
              #}
              #} 
            }
          }
        }
      }
    } 
  }
  
  if (BY_TIME) {
    res = list(df=sliding_tuning_df,
               median=median_sliding_vals,
               mean=mean_sliding_vals)
    
    save(file=sprintf("%s\\%s\\%s\\%s", base_path, "samples", "sliding_tuning", "sliding_df_all_by_time.Rda"), res)
    return(res)
    
  } else {
    save(file=sprintf("%s\\%s\\%s\\%s", base_path, "samples", "sliding_tuning", "peaks_sliding_df_all_3_new_med_FINAL.Rda"), sliding_tuning_df)
    return(sliding_tuning_df)
  }
}

figure_6_new_analysis_tuning_shift_by_laps <- function() {
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\probability_firing_near_analysis\\",figures_path)
  dir.create(write_path)
  
  missing_paths <- c()
  statistics_df <- data.frame()
  wilcox_statistics_df <- data.frame()
  
  for (subfield_name in c("CA1")){#}, "CA3")) {
    if (subfield_name == "CA1") {
      data_paths_ind <- 1:8
    } else {
      data_paths_ind <- 9:18
    }
    sessions_to_use <- 1:16
    prob_fire_df <- data.frame()
    consec_df <- list()
    file_name="new_probability_to_fire_by_pos.Rda"
    true_data_paths <- sprintf("%s\\%s", all_data_paths[data_paths_ind], "equalized")
    
    all_df <- get_all_df(true_data_paths, sessions_to_use)
    estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
    
    metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
    metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
    colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
    
    all_pos_mat <- list()
    
    sizes=c(big=3,
            medium=2.5,
            medium_1=2.3,
            small=2)
    
    for (group_name in c("cyclic_pcs", "cyclic_non_pcs")) {
      for (isSim in c("Real", "Simulation")) {
        all_pos_mat[[sprintf("%s_%s", group_name, isSim)]]  <- 
          matrix(rep(0, times = 24**2), nrow=24)
      }
    }
    
    for (working_path in true_data_paths) {
      
      # extract all paths
      session_paths <- 
        list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
      
      for (ext in c("", "\\fit_simulation")) {
        
        # Run through all paths
        for (tpath in session_paths) {
          
          
          idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
          if (sum(idx == sessions_to_use) == 0){
            next
          }
          
          if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
            print(tpath)
            missing_paths <- c(missing_paths, tpath)
            print("Missing")
            next
          }
          
          if (verbose) {
            print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
          }
          load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=T)
          
          
          if (grepl("CA3", working_path)) {
            mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
            mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
          } else {
            mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
            mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
          }
          
          subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
          sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
          
          
          
          isSimu <- ifelse(ext == "", "Real", "Simulation")
          
          
          for (group_name in c("cyclic_pcs", "cyclic_non_pcs")) {
            # fdf_1 <- apply(final_result[[group_name]]$pos_diff[,1:2], 1, sum)
            # fdf_f <- apply(final_result[[group_name]]$pos_diff[,5:7], 1, sum)
            # fdf <- cbind(fdf_1, final_result[[group_name]]$pos_diff[,3:4])
            # fdf <- cbind(fdf, fdf_f)
            
            if (group_name == "cyclic_pcs") {
              group_ind_tu <- final_result$place_cells_metadata$ind[final_result$place_cells_metadata$cyclic_pcs]
            } else {
              group_ind_tu <- final_result$place_cells_metadata$ind[final_result$place_cells_metadata$cyclic_non_pcs]
            }
            

            mp <- c(apply(final_result$trav_tuning_curve[group_ind_tu,], 1, function(r) {diff(r, lag=5)}))
            
            Consec <- abs(mp[!is.na(mp)])
            
            if (len(mp) == 0) {next}
            
            df <-  c(Session=idx,
                     Subfield=subfield,
                     Mice=mice_str,
                     isSim=isSimu,
                     Direction=ifelse(grepl("Left", tpath), "Left", "Right"),
                     Group=group_name)
            
            consec_df <- append(consec_df,
                                list(Consec))
            
            prob_fire_df <- rbind(prob_fire_df, 
                                  df)
          
          }
        }
      } 
    }
  
    colnames(prob_fire_df) <- c("Session",
                                "SBF",
                                "Mice",
                                "isSim",
                                "Dir",
                                "Group")
    
    
  pdf <- list()
  
  dual_sessions <- as.numeric(prob_fire_df$Session) %% 2 == 0
  dprob_fire_df <- prob_fire_df[dual_sessions,]
  dconsec_df <- consec_df[dual_sessions]
  for (isSim in unique(dprob_fire_df$isSim))  {
    #for (group in unique(dprob_fire_df$Group)) {
      
      for (ses in unique(dprob_fire_df$Session)) {
        #ind <- dprob_fire_df$Group == group & dprob_fire_df$isSim ==isSim & dprob_fire_df$Session == ses
        ind <- dprob_fire_df$isSim ==isSim & dprob_fire_df$Session == ses
      #ind <- dprob_fire_df$Group == group & dprob_fire_df$isSim ==isSim 
      
      plot_df <- data.frame(Session = rep(dprob_fire_df[ind,"Session"], times=unlist(lapply(dconsec_df[ind], function(cdf) {len(cdf)}))),
                            Group = rep(dprob_fire_df[ind,"Group"], times=unlist(lapply(dconsec_df[ind], function(cdf) {len(cdf)})))
                            Shifts = unlist(dconsec_df[ind]))
      
      plot_df$Session <- as.numeric(plot_df$Session) %% 8
      plot_df$Session[as.numeric(plot_df$Session) == 0] <- 8
      plot_df$Session <- as.character(plot_df$Session)
      
      #p <- 
        ggplot(plot_df, aes(x=Shifts)) + 
        geom_density(position = position_nudge(), aes(col=group), alpha=.25) +
        ggtitle(sprintf("%s - %s - %s", isSim, group, ses))
        pdf <- append(pdf, list(p))
      }
    #}
  }
  
  pdf$nrow <- 4
  do.call(plot_grid, pdf)
  
  
  npdf <- list()
  for (isSim in unique(prob_fire_df$isSim))  {
    for (group in unique(prob_fire_df$Group)) {
      
      fdf <- data.frame()
      for (ses in unique(prob_fire_df$Session)) {
        ind <- prob_fire_df$Group == group & prob_fire_df$isSim ==isSim & prob_fire_df$Session == ses
        
        
        df <- data.frame(Med=unlist(lapply(consec_df[ind], median)),
                          Mean=unlist(lapply(consec_df[ind], mean)))
        
        df$Sess <- as.numeric(ses)
        
        fdf <- rbind(fdf, df)
      }
      
      
      fdf$AbsSession <- fdf$Sess %% 8
      fdf$AbsSession[fdf$AbsSession == 0] <- 8
      
      p <- ggplot(fdf, aes(x=AbsSession, y=Med)) + 
        geom_point(stat="summary") + 
        geom_errorbar(stat="summary") + 
        ggtitle(sprintf("%s - %s", isSim, group))
      
      npdf <- append(npdf,
                     list(p))
    }
  }
  
  npdf$nrow <- 2
  do.call(plot_grid, npdf)
    
    for (size_name in names(sizes)) {
      size = sizes[size_name]
      dir.create(sprintf("%s\\%s", 
                         write_path,
                         size_name))
      
      pdf(file=sprintf("%s\\%s\\prob_fire_analysis_%s.pdf",
                       write_path,
                       size_name,
                       subfield_name),
          height=size * 1,
          width=size * 3)
      
      
      plot(pf)
      dev.off()
      
    }
  }
  
  statistics_df$signif <- signif.num(statistics_df$`Pr(>F)`)
  write.csv(file=sprintf("%s\\wilcox_stats.csv",write_path),
            wilcox_statistics_df)
  write.csv(file=sprintf("%s\\anova_stats.csv",write_path),
            statistics_df)
}



figure_6_sliding_tuning_by_window_size <- function(load=T) {
  
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\sliding_tuning_by_window\\",figures_path)
  dir.create(write_path)
  
  
  
  if (load) {
    load(file=sprintf("%s//sliding_tuning_df.Rda",write_path),verbose=T)
  } else {
    get_sliding_tunning_corr_df(file_name = "peaks_sliding_tuning_2.Rda")
  }
  
  sessions_to_use <- 1:16
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  session_groups <- list(mean_ribbon_all_peaks_final_version=1:16)
  #first_2=c(1,2,9,10),
  #first_4=c(1:4,9:13))
  
  sizes <- c(big=3,
             big_1=2.75,
             medium=2.5,
             medium_1=2,
             medium_2=1.75,
             small=1.5)
  
  group_session = T
  
  df_res <- data.frame()
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  
  metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
  metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
  colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
  
  metrics <- unique(sliding_tuning_df$Metric)
  subfields <- unique(sliding_tuning_df$Subfield)
  windows <- unique(sliding_tuning_df$WindowP)
  # cell_funcs <- unique(sliding_tuning_df$CellFunc)
  # window_funcs <- unique(sliding_tuning_df$WindowFunc)
  
  all_jumps <- unique(unlist(lapply(str_split(unique(sliding_tuning_df$WindowP), "_"), function(tuple) {tuple[[2]]})))
  
  all_jumps <- all_jumps[8:len(all_jumps)]
  
  all_windows <- unique(unlist(lapply(str_split(unique(sliding_tuning_df$WindowP), "_"), function(tuple) {tuple[[1]]})))
  
  #all_windows <- all_windows[3:17] # Rest of the bins are irrelevant
  
  statistics_df <- data.frame()
  
  for (session_group_name in names(session_groups)) {
    sessions_to_use <- session_groups[[session_group_name]]
    
    write_path_session <- sprintf("%s\\%s",
                                  write_path, session_group_name)
    
    #dir.create(write_path_session)
    for (metric in metrics) {
      
      metric_path_session <- sprintf("%s\\%s",
                                     write_path_session, metric)
      
      dir.create(metric_path_session)
      
      #for (cf in cell_funcs) {
      #for (wf in window_funcs) {
      for (window in all_windows) {
        
        all_plots <- list()
        minv <- c()
        maxv <- c()
        
        
        df_to_use_outer <- 
          sliding_tuning_df[sliding_tuning_df$Session %in% sessions_to_use &
                            sliding_tuning_df$Metric == metric &
                            grepl(window, sliding_tuning_df$WindowP),]
        
        
        
        
        pval_grps <- c()
        error_plots <- list()
        max_y_all_errors <- c()
        error_pvals <- c()
        
        errors <- c()
          
          for (grp in c("Pcs", "Npcs")) {
            minv<- min(df_to_use_outer$Sliding, na.rm=T)
            maxv <- max(df_to_use_outer$Sliding, na.rm=T)
            df_grp_sbf <- data.frame()
            pval_all <- c()
            anova_all <- data.frame()
            
            for (jmp in all_jumps) {         
            tmp_for_test <- list()
              
              for (isSim in c("Real", "Simulation")) {
                
                df_to_use_outer$Session = df_to_use_outer$Session %% 8
                df_to_use_outer$Session[df_to_use_outer$Session == 0] <- 8
                search_window <- sprintf("%s_%s", window, jmp)
                df_to_use <- df_to_use_outer[df_to_use_outer$WindowP == search_window &
                                               df_to_use_outer$Group == grp &
                                               df_to_use_outer$isSim == isSim, ]
                
                print(nrow(df_to_use))
                
                if (nrow(df_to_use) == 0) {
                  next
                }
                
                
                
                  tmp <- 
                    ddply(df_to_use, .(isSim),  
                          function(issim_df) {
                            by_ses <- 
                              ddply(issim_df, .(Session), 
                                    function(session_df) {
                                      return(c(mean(session_df$Sliding, na.rm=T)))
                                    })
                            return(by_ses)
                          })
                  
                  colnames(tmp) <- c("isSim", "Session", "mu")
                  
                  tmp2 <- 
                    ddply(tmp, .(isSim),  
                          function(issim_df) {
                            return(c(mean(issim_df$mu, na.rm=T),
                                     sem(issim_df$mu)))
                          })
                  
                  colnames(tmp2) <- c("isSim", "mu", "sd")
                  
                  tmp2$jmp <- rep(as.numeric(str_split(jmp, "J")[[1]][2]), 
                                     nrow(tmp2))
                  
                  df_grp_sbf <- rbind(df_grp_sbf,
                                      tmp2)
                  
                  tmp_for_test <- append(tmp_for_test,
                                         list(tmp$mu))
                  
                  tmp$jmp <-  rep(as.numeric(str_split(jmp, "J")[[1]][2]), 
                                     nrow(tmp))
                  anova_all <- rbind(anova_all,
                                     tmp)
                }
              
              
              tmp3 <- do.call(cbind, tmp_for_test)
              pval_all <- c(pval_all,
                            wilcox.test(tmp3[,1], tmp3[,2])$p.value)
            }
            
            df_grp_sbf <- as.data.frame(df_grp_sbf)
            
            
            minv <- c(minv, min(df_grp_sbf$mu - df_grp_sbf$sd))
            maxv <- c(maxv, max(df_grp_sbf$mu + df_grp_sbf$sd))
            
            
            g <- 
              ggplot(df_grp_sbf, aes(x=jmp, y=mu, group=isSim, color=isSim, fill=isSim)) +
              geom_line() +
              #geom_point(size=2) +
              #geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd))+
              geom_ribbon(aes(ymin=mu-sd, ymax=mu+sd), color=NA, alpha=.3) + 
              theme_classic() +
              ggtitle(sprintf("%s",
                              grp)) + 
              scale_color_manual(values=c("#652D90", "#F05A28"), breaks=c("Real", "Simulation")) + 
              scale_fill_manual(values=c("#652D90","#F05A28"), breaks=c("Real", "Simulation")) + 
              base_plot_theme +
              theme(text=element_text(size=14, color="black"),
                    axis.text = element_text(color="black"),
                    plot.title = element_text(size=9))  +
              xlab("Window size (frames)") +
              ylab("Peak shift (cm)")
            #ylab(metric) # + ylim(c(minv, maxv))
            
            
            all_plots <- append(all_plots,
                                list(g))
            
            pval_grps <- rbind(pval_grps, pval_all)
            
            error <- df_grp_sbf[df_grp_sbf$isSim == "Real",][order(df_grp_sbf[df_grp_sbf$isSim == "Real",]$jmp),"mu"] - 
              df_grp_sbf[df_grp_sbf$isSim != "Real",][order(df_grp_sbf[df_grp_sbf$isSim != "Real",]$jmp),"mu"] 
            errors <- rbind(errors,
                            error)
          }
          
          
          errors <- t(errors ** 2)
          colnames(errors) <- c("PC", "NPC")
          melted_errors_df <- melt(errors)
          colnames(melted_errors_df) <- c("n", "Group", "MSE")
          max_y <- max(colMeans(errors) + apply(errors, 2, sem))
          wilk_t <-  wilcox.test(errors[,1], errors[,2], alternative="less", paired=T, corrected=F)
          error_pvals <- c(error_pvals,
                           wilk_t$p.value)
          
          statistics_df <- rbind(statistics_df,
                                 data.frame(
                                   subfield="CA1",
                                   pval=wilk_t$p.value,
                                   method=wilk_t$method,
                                   statistic=wilk_t$statistic,
                                   alternative=wilk_t$alternative,
                                   signif=signif.num((wilk_t$p.value)),
                                   jump=jmp,
                                   #cf=cf,
                                   #wf=wf,
                                   metric=metric,
                                   mean_pc=apply(errors, 2, mean)[1],
                                   mean_npc=apply(errors, 2, mean)[2],
                                   sd_pc=apply(errors, 2, sd)[1],
                                   sd_npc=apply(errors, 2, sd)[2],
                                   sem_pc=apply(errors, 2, sem)[1],
                                   sem_npc=apply(errors, 2, sem)[2]))
          
          gerror <-
            ggplot(melted_errors_df, aes(x=Group, y=MSE)) +
            geom_bar(stat="summary", width=.7) +
            geom_errorbar(stat="summary", width=.3) +
            #geom_jitter(aes(color=Group),position=position_jitterdodge(.5))
            base_plot_theme + 
            scale_y_continuous(expand=c(0,0)) +
            #ggtitle(sbf) + 
            theme(text=element_text(size=14, color="black"),
                  axis.text = element_text(color="black"),
                  plot.title = element_text(size=9)) + 
            xlab("")
          
          max_y_all_errors <- c(max_y_all_errors, max_y + .25)
          
          error_plots <- append(error_plots,
                                list(gerror))
        
        
        all_plots <- lapply(all_plots, 
                            function(p) {p + ylim(c(min(minv), max(maxv)))})
        all_plots$nrow <- 1
        gf <- do.call(plot_grid, all_plots)
        
        max_y <- max(max_y_all_errors)
        error_plots <- lapply(1:len(error_plots),
                              function(p_idx) {
                                p <- error_plots[[p_idx]]
                                p <- p + 
                                  geom_line(data=data.frame(x=c(1,1,2,2),
                                                            y=c(max_y  + .15, max_y + .25, max_y  +.25, max_y + .15)),
                                            aes(x=x,y=y)) + 
                                  geom_text(label=signif.num(error_pvals[p_idx]),
                                            x=1.5,
                                            y=max_y + .15)
                                
                                return(p)
                              })
        error_plots$nrow <- 2
        gf_error <- do.call(plot_grid, error_plots)
        for (size_name in names(sizes)) {
          size = sizes[[size_name]]
          write_path_session <- sprintf("%s\\%s\\",write_path, size_name)
          dir.create(write_path_session)
          write_path_session <- sprintf("%s\\%s\\%s",write_path, size_name, session_group_name)
          dir.create(write_path_session)
          metric_path_session <- sprintf("%s\\%s",write_path_session, metric)
          dir.create(metric_path_session)
          
          
          pdf(sprintf("%s//Jump_%s_cell.pdf",
                      metric_path_session,
                      window),
              height=size *1,
              width=size * 2)
          plot(gf) 
          dev.off()
          
          pdf(sprintf("%s//Jump_%s_cell_error.pdf",
                      metric_path_session,
                      window),
              height=size * 1,
              width=size *  1)
          plot(gf_error) 
          dev.off()
          
        }
      }          

    }
  }
  
  write.csv(file=sprintf("%s//statistics.csv",write_path), statistics_df)
  save(file=sprintf("%s//sliding_tuning_df.Rda",write_path), sliding_tuning_df)
}