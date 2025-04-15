


figure_6_corlineplots_new <- function() {
  
  load(sprintf("%s\\samples\\tuning_dataframes\\KS_stability_dataframes.R", 
               base_path), 
       verbose=T)
  
  write_path_base <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path_base)
  write_path_base <- sprintf("%s\\figure_6\\corline_plots_new_all_100\\",figures_path)
  dir.create(write_path_base)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  sessions_to_use=c(1:16)
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  
  metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
  metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
  colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.3,
          small=2)
  
  
  session_groups <- list(all=1:16)
  

  
  mrs_all <- c()
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
            non_pc_error_all <- c()
            pc_error_all <- c()
            non_pcs <- c()
            
            for (mice in unique(subsample_cor_df$Mice)) {
              for (session in unique(subsample_cor_df$Session)){
                for (dir in unique(subsample_cor_df$Direction)) {
                  for (group in c("cyclic_pcs", "cyclic_non_pcs")) {
                    
                    est = 
                    metadata_df[metadata_df$Session == session & 
                                metadata_df$Mice == mice & 
                               metadata_df$Direction == dir,]$Est
                    
                    if (est != 1) {
                      #print(sprintf("%s - %s - %s estimation is not 100%% (%f)", session, mice, dir, est))
                      next
                    }
                    
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
                      pc_error = (df2[df2$isSim == "Real","Corr"]- df2[df2$isSim != "Real","Corr"]) ** 2
                      pc_error_all <- c(pc_error_all,
                                        pc_error)
                                        
                      
                      
                      pcs <- c(pcs, c(df2[,"Corr"]))
                    } else {
                      non_pcs_real <- c(non_pcs_real, df2[df2$isSim == "Real","Corr"])
                      non_pcs_sim <- c(non_pcs_sim, df2[df2$isSim != "Real","Corr"]) 
                      non_pc_error = (df2[df2$isSim == "Real","Corr"]- df2[df2$isSim != "Real","Corr"]) ** 2
                      non_pc_error_all <- c(non_pc_error_all,
                                            non_pc_error)
                      non_pcs <- c(non_pcs, c(df2[,"Corr"]))
                    }
                    
                    real <- c(real, df2[df2$isSim == "Real","Corr"])
                    sim <- c(sim, df2[df2$isSim != "Real","Corr"])
                    
                    
                  }
                }
              }
            }
            
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
            
            mrs <-  c(df_name,
                      subsample,
                      mean(pc_error_all),
                      mean(non_pc_error_all),
                      wilcox.test(pc_error_all, non_pc_error_all, paired=T)$p.value,
                      t.test(pc_error_all, non_pc_error_all, paired=T)$p.value,
                      prop)
            mrs_all <- rbind(mrs_all,
                             mrs)
            
            print(mrs)
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


figure_6_sliding_by_time_new <- function() {
  write_path <- sprintf("%s\\figure_6\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_6\\sliding_by_time\\",figures_path)
  dir.create(write_path)
  
  load(sprintf("%s\\samples\\sliding_tuning\\not_100_sliding_df_all_by_time.Rda", 
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
  
  
  # true_data_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  # sessions_to_use=c(1:16)
  # all_df <- get_all_df(true_data_paths, sessions_to_use)
  # estimated_res_cyclic <- get_estimated_df(all_df, true_data_paths, sessions_to_use, estimated_df_folder="KS_simulations_likelihood_c")
  #  
  # metadata_df <- as.data.frame(estimated_res_cyclic$session_df)
  # metadata_df <- cbind(metadata_df, unlist(estimated_res_cyclic$estimation_list))
  # colnames(metadata_df) <- c("Session", "Mice", "Direction", "Subfield", "Est")
  # metadata_df <- metadata_df[metadata_df$Est == 1,]
  # 
  # session_groups <- list(all_non_100=1:16)
  # 
  
  # new_sliding_tuning_df <- data.frame()
  # 
  # ind_all <- c()
  # 
  # for (session_to_use_idx in 1:nrow(metadata_df)) {
  #   session_to_use <- metadata_df[session_to_use_idx,]
  #   ind <- which(res$df$Session %in% session_to_use$Session &
  #                res$df$Mice %in% session_to_use$Mice &
  #                res$df$Direction %in% session_to_use$Direction)
  #   
  #   ind_all <- c(ind_all, ind)
  #   
  # }
  # 
  # res$df <- res$df[ind_all,]
  # res$median <- res$median[ind_all]
  # res$mean <- res$mean[ind_all]
  # 
  # save(file=sprintf("%s\\samples\\sliding_tuning\\not_100_sliding_df_all_by_time.Rda", 
  #                   base_path),
  #      res)
  # 
  # dim(metadata_df)
  
  
  
  
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



figure_6_non_place_cells_rasters <- function() {
  dir.create(sprintf("%s\\figure_6\\",figures_path))
  write_path <- sprintf("%s\\figure_6\\specific_rasters",figures_path)
  dir.create(write_path)
  
  npc_csv <- read.csv(sprintf("%s\\npc_simulation.csv", base_path))
  
  
  tabled <- table(npc_csv[,1:2])
  
  npc_to_sim <- list()
  for(i in 1:nrow(tabled)) {
    
    
    dataset_sessions <- tabled[i,]
    dataset_name <- as.numeric(rownames(tabled))[i]
    dataset_sessions <- as.numeric(names(dataset_sessions)[dataset_sessions != 0])
    
    for (j in dataset_sessions) {
      events <- unique(npc_csv[npc_csv$dataset == dataset_name & npc_csv$session == j,]$events)
      ind <- unique(npc_csv[npc_csv$dataset == dataset_name & npc_csv$session == j,]$cell_num)
      
      sim_list <- list(dataset=dataset_name,
                       session=j,
                       events=events,
                       cells=ind)
      npc_to_sim <- append(npc_to_sim,
                           list(sim_list))
    }
  }

  plot_events = T
  plot_cells = F
  
  
  
  for (dataset in npc_to_sim) {
    
    session_length <- 20
    p <- all_data_paths[dataset$dataset]
    a <- get_spike_train_and_stim_trace_from_path(p,
                                                  dataset$session)
    original_spike_train <- a[[1]]
    stim_trace <- a[[2]]
    
    run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * session_length), 
                         Position=stim_trace * 4) 
    run_df$Position[which(run_df$Position == 4)] <- 0
    
    tmp <- compute_tuning(original_spike_train, stim_trace)
    SI <- compute_SI(tmp[[1]], tmp[[2]], rowMeans(original_spike_train) / dt)[[1]]
    
    if (grepl("CA3", p)) {
      mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', p))
      mice_str <- substr(p, mice_str_index, mice_str_index+4) 
    } else {
      mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', p))
      mice_str <- substr(p, mice_str_index, mice_str_index+3) 
    }
    
    dir <- ifelse(grepl("Left", p), "L", "R")
    issim <- "_simulation"
    mice_path <- sprintf("%s\\mice_%s_%s_%d%s",
                         write_path,
                         mice_str,
                         dir,
                         dataset$session,
                         issim)
    dir.create(mice_path)
    

    if (plot_cells) {
      for (cell in dataset$cells) {
        
        cell_group_path <- sprintf("%s\\real_cells", mice_path, cell)
        dir.create(cell_group_path)
        dir.create(sprintf("%s\\regular_size", cell_group_path))
        dir.create(sprintf("%s\\big_size", cell_group_path))
        dir.create(sprintf("%s\\med_size", cell_group_path))
        
        #idx <- round(runif(1,1,nrow(spike_train)))
        firing_ind <- which(original_spike_train[cell,] > 0)
        
        cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
                              Times=run_df$Time[firing_ind])
        
        
        tuning_curve <- (unlist(lapply(sort(unique(stim_trace)), function (sb) {mean(original_spike_train[cell,which(stim_trace == sb)]) / dt}))) 
        barp <- barplot(tuning_curve)
        smoothed_tuning <- smooth.spline(barp, tuning_curve, all.knots = T, lambda=1e-4)
        smoothed_tuning$y[smoothed_tuning$y < 0] <- 0
        pos <- sort(unique(stim_trace)) * 4; pos[pos == 4] <- 0
        smoothed_df <- data.frame(FR=smoothed_tuning$y,
                                  Positions=pos)
        
        smoothed_df$FR <- smoothed_df$FR / sum(smoothed_df$FR)
        
        
        title = sprintf("SI (%f) E(%d)",
                        SI[[cell]],
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
        
        
        png(sprintf("%s\\regular_size\\real_neur_%d.png", cell_group_path, cell), unit="in", height=5, width=1.5, res=300); 
        plot(g); dev.off()   
        
        pdf(sprintf("%s\\regular_size\\real_neur_%d.pdf", cell_group_path, cell), height=5, width=1.5); 
        plot(g); dev.off()   
        
        pdf(sprintf("%s\\big_size\\real_neur_%d.pdf", cell_group_path, cell), height=5, width=1.5); 
        plot(g_big); dev.off()   
        
        pdf(sprintf("%s\\med_size\\real_neur_%d.pdf", cell_group_path, cell), height=5, width=1.5); 
        plot(g_med); dev.off()
      }
    }
    
    if (plot_events) {
      for (n_events in dataset$events) {
        spike_train <- original_spike_train
        fr <- rowMeans(spike_train) / dt
        processed_real <- preprocess_spike_train(spike_train, stim_trace)
        
        true_cells_spike_train <- processed_real$working_cells_spike_train
        true_firing_rate <- processed_real$working_firing_rate
        true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
        
        #params <- get_fit_params(path, estimation_path="simulations_2_new_new", likelihood_path = "simulations_3_new", pct_range=seq(0.5,1,by=0.1))
        params <- get_fit_params(sprintf("%s\\equalized\\session_%d", p, dataset$session))
        
        print(params)
        simulated_tuning_curve <-
          generate_tuning_curves_cost(n = 15000,
                                      percentage = 0.05,
                                      average_width = params$params["average_width"], 
                                      sd_width = params$params["sd_width"],
                                      fixed_fr=sample(rep(fr, 500), 15000),
                                      noise=params$params["noise"],
                                      double_peak_pct = params$params["double_peak_pct"],
                                      n_bins=24,
                                      plot=F)
        
        pois_factor = 1.5
        
        generated_spike_train <- 
          generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                     stim_trace = stim_trace,
                                     factor=pois_factor,
                                     fs=1)
        
        spike_train <- generated_spike_train[which(rowSums(generated_spike_train) == n_events),]
        simulated_tuning_curve <- simulated_tuning_curve[which(rowSums(generated_spike_train) == n_events),]
        
        
        tmp <- compute_tuning(spike_train, stim_trace)
        SI <- compute_SI(tmp[[1]], tmp[[2]], rowMeans(spike_train) / dt)[[1]]
        
        
        
        
        mice_group_path <- sprintf("%s\\%d", mice_path, n_events)
        dir.create(mice_group_path)
        dir.create(sprintf("%s\\regular_size", mice_group_path))
        dir.create(sprintf("%s\\big_size", mice_group_path))
        dir.create(sprintf("%s\\med_size", mice_group_path))
        
        for (idx in 1:nrow(spike_train)) {
          
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
